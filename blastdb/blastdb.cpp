#include <db_cxx.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

#define USAGE "\
usage: blastdb <database file> <blast file> ...\n\n"

#define MAX_LINE 1000


/*
 NCBI BLASTALL 2.2.10 -m 8 tab-delimited output
 Fields: 
 0. Query id, 
 1. Subject id, 
 2. % identity, 
 3. alignment length, 
 4. mismatches, 
 5. gap openings, 
 6. q. start, 
 7. q. end, 
 8. s. start, 
 9. s. end, 
 10. e-value, 
 11. bit score
*/


int Split(char *str, char delim, char **tokens, int max)
{
    char *token = str;
    int ntokens = 0;
    
    // walk down string
    for (char *p=str; *p; p++) {
        if (*p == delim) {
            // save token, cut string
            tokens[ntokens++] = token;
            *p = '\0';
            token = p+1;
            
            // only find max tokens
            if (ntokens >= max)
                return ntokens;
        }
    }
    
    // save last token
    tokens[ntokens++] = token;
    
    return ntokens;
}




class BlastDb
{
public:
    BlastDb(const char *filename) :
        m_db(NULL, 0),
        m_nhits(0)
    {
        try {
            // allow duplicate keys
            m_db.set_flags(DB_DUP);
        
            // open database
            m_db.open(NULL,                // Transaction pointer 
                      filename,            // Database file name 
                      NULL,                // Optional logical database name
                      DB_HASH,             // Database access method
                      DB_CREATE,           // Open flags
                      0);                  // File mode (using defaults)
        } catch(DbException &e) {
            fprintf(stderr, "bdb error\n");
        } catch(std::exception &e) {
            fprintf(stderr, "unknown error\n");
        }
    }
    
    virtual ~BlastDb()
    {
        // close database
        m_db.close(0);
    }
    
    
    inline void Add(const char *key, const char *value) {
        // do not store the trailing NULL
        Dbt dbkey((void*) key, strlen(key));
        Dbt dbvalue((void*) value, strlen(value));
        m_db.put(NULL, &dbkey, &dbvalue, 0);
    }
    
    
    bool AddBlastFile(const char *filename) {
        FILE *infile;        
        
        if ((infile = fopen(filename, "r")) == NULL) {
            fprintf(stderr, "could not open '%s'\n", filename);
            return false;
        }
        
        char line[MAX_LINE+1];
        
        while (fgets(line, MAX_LINE, infile)) {
            if (line[0] != '#')
                AddBlastHit(line);
        }
        
        fclose(infile);
    }
    
    
    void AddBlastHit(char *line)
    {
        char query[100];
        char line2[MAX_LINE+1];
        
        // get query name
        int i=0;
        for (; i<100 && line[i] != '\t'; i++)
            query[i] = line[i];
        query[i] = '\0';
        
        // remove newline from line
        int end = strlen(line)-1;
        if (line[end] == '\n')
            line[end] = '\0';
        
        // add new hit
        Add(query, line);
        
        // split line into tokens
        char *tokens[12];
        int ntokens = Split(line, '\t', tokens, 12);
        
        if (ntokens == 12) {
            // flip query & subject 
            line2[0] = '\0';
            
            int perm[] = {1,0,2,3,4,5,8,9,6,7,10,11};
            
            for (int i=0; i<11; i++) {
                strncat(line2, tokens[perm[i]], MAX_LINE);
                strncat(line2, "\t", MAX_LINE);
            }
            
            strncat(line2, tokens[11], MAX_LINE);
            
            // add flipped hit
            Add(tokens[1], line2);
        } else {
            fprintf(stderr, "bad line at hit %d\n", m_nhits);
            printf("%d\n", ntokens);
        }
        
        
        // count number of added hits
        m_nhits++;
        if (m_nhits % 1000 == 0) {
            fprintf(stderr, "added %d hits\n", m_nhits);
        }
    }
    
    Db m_db;
    int m_nhits;
};


int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, USAGE);
        exit(1);
    }

    char *dbfile = argv[1];    
    BlastDb db(dbfile);
    
    // add blast files
    for (int i=2; i<argc; i++) {
        fprintf(stderr, "processing '%s'\n", argv[i]);
        db.AddBlastFile(argv[i]);
    }
    
}
