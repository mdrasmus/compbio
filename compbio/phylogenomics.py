
# python imports
import sys, os, re
from pysqlite2 import dbapi2 as sqlite

# rasmus imports
from rasmus import tablelib, stats, treelib, util

# compbio imports
from . import genecluster, phylo, fasta, gff


def tableExists(cur, tablename):
   cur.execute("""SELECT name from sqlite_master 
                  WHERE name = '%s';""" % tablename)
   return cur.fetchone() != None



class PhyloDb:
    def __init__(self, dbfile=None, famfile=None, smapfile=None,
                 genenamefile=None, streefile=None,
                 baseDir=None,
                 treeFileExt=None,
                 fastaFileExt=None):
        self.fams = genecluster.FamilyDb(famfile)
        self.gene2species = phylo.read_gene2species(smapfile)
        self.genenames_tab = tablelib.read_table(genenamefile)
        self.gene2name = self.genenames_tab.lookup("id")
        self.stree = treelib.read_tree(streefile)
        self.baseDir = baseDir
        self.treeFileExt = treeFileExt
        self.fastaFileExt = fastaFileExt
        
        # open database
        self.con = sqlite.connect(dbfile, isolation_level="DEFERRED")
        self.cur = self.con.cursor()
    

    def close(self):
        self.con.commit()
        self.con.close()


    #==========================================================
    # Genes

    def makeGenesTable(self):
        """make genes table"""
        
        self.cur.execute("""
            CREATE TABLE Genes (
                geneid CHAR(20) PRIMARY KEY,
                common_name CHAR(20),
                species CHAR(20),
                chrom CHAR(20),
                start INTEGER,
                end INTEGER,
                strand INTEGER,        
                description VARCHAR(1000),
                famid CHAR(20)
            );""")

        self.cur.execute("""CREATE UNIQUE INDEX IndexGenes
                            ON Genes (geneid);""")
        self.cur.execute("""CREATE INDEX Index2Genes
                            ON Genes (famid);""")

    def clearGenes(self):
        if tableExists(self.cur, "Genes"):
            self.cur.execute("DROP TABLE Genes");
    
    
    def addGenes(self, species, gff_files, region_filter=lambda x: x):
        """populate genes table"""

        # clear Genes Table
        if not tableExists(self.cur, "Genes"):
            self.makeGenesTable()
        

        dups = set()
        
        util.tic("add genes")
        for sp, gff_file in zip(species, gff_files):
            for region in gff.read_gff(gff_file, regionFilter=region_filter):
                gene = region.data["ID"]
                #gene = row["name"]

                if gene in self.fams.genelookup:                
                    famid = self.fams.getFamid(gene)
                    if len(self.fams.getGenes(famid)) < 2:
                        famid = "NONE"
                else:
                    famid = "NONE"
                
                
                if gene in dups:
                    continue
                dups.add(gene)
                
                assert region.start <= region.end

                if gene in self.gene2name:
                    common = self.gene2name[gene]["name"]
                    desc = self.gene2name[gene]["description"]
                else:
                    common = ""
                    desc = ""

                cmd = """INSERT INTO Genes VALUES 
                               ("%s", "%s", "%s", "%s", %d, %d, %d, "%s", "%s");""" % \
                            (gene, 
                             common,
                             self.gene2species(gene),
                             region.seqname,
                             region.start,
                             region.end,
                             region.strand,
                             desc.replace('"', ''),
                             famid)

                self.cur.execute(cmd)
        util.toc()


    #=====================================================
    # Families

    def makeFamiliesTable(self):
        """create families table"""       

        self.cur.execute("""CREATE TABLE Families (
                            famid CHAR(20) PRIMARY KEY,
                            gene_names CHAR(20),
                            gene_rate FLOAT,
                            tree_len FLOAT,
                            median_gene_len FLOAT,
                            dup INTEGER,
                            loss INTEGER,
                            genes INTEGER,
                            description VARCHAR(1000)
                           );
                        """)

        self.cur.execute("""CREATE UNIQUE INDEX IndexFamilies
                            ON Families (famid);""")

    def clearFamilies(self):
        if tableExists(self.cur, "Families"):
            self.cur.execute("DROP TABLE Families");

    
    def addFamilies(self, eventsfile, discard=[]):
        
        if not tableExists(self.cur, "Families"):
            self.makeFamiliesTable()
    
        util.tic("add families")
        events_tab = tablelib.read_table(eventsfile)
        events_lookup = events_tab.lookup("partid")
        familyGeneNames = self.makeFamilyGeneNames()
        discard = set(discard)

        for row in events_tab:
            famid = row["partid"]
            if famid in discard:
                util.logger("discarding '%s'" % famid)
                continue
            
            tree = treelib.read_tree(self.getTreeFile(famid))
            treelen = sum(x.dist for x in tree)
            seqs = fasta.read_fasta(self.getFastaFile(famid))
            seqlen = stats.median(map(len, seqs.values()))

            self.cur.execute("""INSERT INTO Families VALUES 
                                ("%s", "%s", %f, %f, %f, %d, %d, %d,
                                 "%s");""" % 
                        (row["partid"], 
                         familyGeneNames.get(row["partid"], ("", ""))[0],
                         row["famrate"], treelen, seqlen * 3,
                         row["dup"], row["loss"], row["genes"],
                         familyGeneNames.get(row["partid"], ("", ""))[1]))
        util.toc()
    
    
    
    def getFamDescription(self, descriptions):
    
        # TODO: remove this hardcoding
        rmdesc = set([
            "",
            "Predicted ORF from Assembly 19",
            "Predicted ORF in Assemblies 19 and 20",
            "ORF Predicted by Annotation Working Group",
            "possibly spurious ORF (Annotation Working Group prediction)"])

        descs = []
        for d in descriptions:
            descs.extend(d.split("; "))
        descs = filter(lambda x: x not in rmdesc, descs)

        items = util.hist_dict(descs).items()
        items.sort(key=lambda x: x[1], reverse=True)

        desc = "; ".join(["%s[%d]" % item for item in items])
        return desc
    
    
    def makeFamilyGeneNames(self):
        """Tries to name and describe a family using its genes"""
    
        self.cur.execute("""SELECT g.famid, g.common_name, g.description
                            FROM Genes g
                         """)
        
        fams = util.groupby(lambda x: x[0], self.cur)
        
        familyGeneNames = {}
        for famid, rows in fams.iteritems():
            names = util.unique(["".join([i for i in x
                                     if not i.isdigit() and i != "-"])
                            for x in util.cget(rows, 1)
                            if x != ""])
            names.sort()
            
            description = self.getFamDescription(util.cget(rows, 2))
                        
            familyGeneNames[famid] = (",".join(names), description)
        return familyGeneNames
    
    
    #============================================
    # Phylogenetic Events
    
    def makeEventsTable(self):
        
        self.cur.execute("""
            CREATE TABLE Events (
                famid CHAR(20),
                species CHAR(20),
                genes INTEGER,                
                dup INTEGER,
                loss INTEGER,
                appear INTEGER
            );""")
    
    
        self.cur.execute("""CREATE INDEX EventsIndex
                            ON Events (famid);""")
    
    
    
    def clearEvents(self):
        if tableExists(self.cur, "Events"):
            self.cur.execute("DROP TABLE Events")
        
    
    def addEvents(self, eventsfile):
        
        if not tableExists(self.cur, "Events"):
            self.makeEventsTable()
        
        util.tic("add events")
        events_tab = tablelib.read_table(eventsfile)
        events_lookup = events_tab.lookup("partid")
        
        self.cur.execute("SELECT famid FROM Families;")
        famids = [x[0] for x in self.cur]
        
        for famid in famids:
            if famid not in events_lookup:
                continue
            row = events_lookup[famid]
            
            for sp in self.stree.nodes:
                sp = str(sp)
            
                self.cur.execute("""INSERT INTO Events VALUES 
                                ("%s", "%s", %d, %d, %d, %d);""" % 
                        (famid, sp,
                         row[sp+"-genes"],
                         row[sp+"-dup"],
                         row[sp+"-loss"],
                         row[sp+"-appear"]))
        util.toc()
        
    
    #============================================
    # GO terms
    
    def makeGoTermsTable(self):
        """add go term table"""

        self.cur.execute("""
        CREATE TABLE GoTerms (
            goid CHAR(20) PRIMARY KEY,
            name VARCHAR(100)
        );""")

        self.cur.execute("""CREATE UNIQUE INDEX IndexGoTerms
                            ON GoTerms (goid);""")
        
        # add gene2goterm table
        self.cur.execute("""CREATE TABLE GeneGoTerms (
                            geneid CHAR(20),
                            goid CHAR(20)
                            );""")

        self.cur.execute("""CREATE INDEX IndexGeneGoTerms
                            ON GeneGoTerms (goid);""")
        self.cur.execute("""CREATE INDEX Index2GeneGoTerms
                            ON GeneGoTerms (geneid);""")
                            

    def clearGoTerms(self):
        if tableExists(self.cur, "GoTerms"):
            self.cur.execute("DROP TABLE GoTerms;")
        if tableExists(self.cur, "GeneGoTerms"):
            self.cur.execute("DROP TABLE GeneGoTerms;")
        

    def addGoTerms(self, gofile):
        
        if not tableExists(self.cur, "GoTerms"):
            self.makeGoTermsTable()
        
        util.tic("add go terms")
        goterms = tablelib.read_table(gofile)
        goterms_lookup = goterms.groupby("orf")
        goterms_bygoid = goterms.groupby("goid")
        
        for goterm in goterms_bygoid:
            term = goterms_bygoid[goterm][0]

            if '"' in term["term"]:
                print term

            self.cur.execute("""INSERT INTO GoTerms VALUES ("%s", "%s")""" %
                             (term["goid"], term["term"]))

        
        for gene, terms in goterms_lookup.iteritems():
            for term in terms:
                self.cur.execute("""INSERT INTO GeneGoTerms VALUES ("%s", "%s");""" %
                                 (gene, term["goid"]))
        util.toc()


    #================================================
    # Pfams    

    def makePfamTable(self):
        """make pfam domains table"""
        
        self.cur.execute("""CREATE TABLE PfamDomains (
                            pfamid CHAR(20) PRIMARY KEY,
                            name CHAR(20),
                            description VARCHAR(100)
                            );""")

        self.cur.execute("""CREATE UNIQUE INDEX IndexPfamDomains
                            ON PfamDomains (pfamid);""")
        
        self.cur.execute("""CREATE TABLE GenePfamDomains (
                            geneid CHAR(20),
                            pfamid CHAR(20),
                            start INTEGER,
                            end INTEGER,
                            score FLOAT,
                            evalue FLOAT
                            );""")

        self.cur.execute("""CREATE INDEX IndexGenePfamDomains
                            ON GenePfamDomains (pfamid);""")
        self.cur.execute("""CREATE INDEX Index2GenePfamDomains
                            ON GenePfamDomains (geneid);""")


    def clearPfams(self):
        if tableExists(self.cur, "PfamDomains"):
            self.cur.execute("DROP TABLE PfamDomains;")
        if tableExists(self.cur, "GenePfamDomains"):
            self.cur.execute("DROP TABLE GenePfamDomains;")

    def addPfams(self, pfamfile):

        if not tableExists(self.cur, "PfamDomains"):
            self.makePfamTable()
        
        pfams = tablelib.read_table(pfamfile)
        
        for row in pfams:         
            self.cur.execute("""INSERT INTO PfamDomains VALUES ("%s", "%s", "%s");""" %
                             (row['pfamid'], 
                              row['pfam_name'], 
                              row['pfam_description']))




    def addPfamGenes(self, pfamfile):
        """add pfam domains"""

        if not tableExists(self.cur, "PfamDomains"):
            self.makePfamTable()      

        util.tic("add pfam domains")

        pfams = tablelib.read_table(pfamfile)
        
        for row in pfams:
            name = re.sub("\..*$", "", row["pfam_acc"])
           
            self.cur.execute("""INSERT INTO GenePfamDomains VALUES
                                ("%s", "%s", %d, %d, %f, %f);""" %
                             (row["locus"], name, row["start"], 
                              row["end"], row["score"], row["evalue"]))
        util.toc()


    def getTreeFile(self, famid):
        return os.path.join(self.baseDir, famid, famid + self.treeFileExt)

    def getFastaFile(self, famid):
        return os.path.join(self.baseDir, famid, famid + self.fastaFileExt)
