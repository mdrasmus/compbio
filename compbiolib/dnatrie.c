#include <stdlib.h>
#include "dnatrie.h"


#ifdef __cplusplus
extern "C" {
#endif



static char *int2dna = "ACGT";
static int dna2int_table [256] = 
{
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 9
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 19
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 29
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 39
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 49
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 59
    -1, -1, -1, -1, -1,  0, -1,  1, -1, -1,   // 69
    -1,  2, -1, -1, -1, -1, -1, -1, -1, -1,   // 79
    -1, -1, -1, -1,  3, -1, -1, -1, -1, -1,   // 89
    -1, -1, -1, -1, -1, -1, -1,  0, -1,  1,   // 99
    -1, -1, -1,  2, -1, -1, -1, -1, -1, -1,   // 109
    -1, -1, -1, -1, -1, -1,  3, -1, -1, -1,   // 119
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 129
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 139
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 149
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 159
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 169
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 179
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 189
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 199
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 209
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 219
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 229
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 239
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 249
    -1, -1, -1, -1, -1, -1                    // 255
};

#define dna2int(dna) (dna2int_table[(int) (dna)])


/* insert a dna sequence into a Trie tree 
   To initialize a new tree, pass a pointer to a NULL pointer.
*/
void dnatrie_insert(DnaTrie** tree_ptr, char* seq, unsigned char strand)
{
    while (1) {
    	if (*tree_ptr == NULL)
	    	*tree_ptr = (DnaTrie*) calloc(1, sizeof(DnaTrie));

    	if (seq[0] == '\0') {
            /* NOTE: may not be needed */
		    (*tree_ptr)->strand = (*tree_ptr)->strand | strand; 
            break;
    	} else if (dna2int(seq[0]) != -1) {
	    	tree_ptr = &((*tree_ptr)->next[dna2int(seq[0])]);
            seq++;
        }
    }
}

/* query whether a squence is in the tree */
int dnatrie_query(DnaTrie *node, char *seq)
{
    int depth = 0;
    
    while (seq[0] && dna2int(seq[0]) != -1) {
        /* follow next pointer */
        node = node->next[dna2int(seq[0])];
    
        /* quit if no pointer */
        if (!node)
            break;
        
        /* move down query */
        seq++;
        depth++;
    }
    
    /* return how much of the query matches the tree */
    return depth;
}


/* delete a Trie tree */
void dnatrie_delete(DnaTrie* tree_ptr)
{
	int i;
	if (tree_ptr == NULL)
		return;

	for (i=0; i<4; i++)
		dnatrie_delete(tree_ptr->next[i]);

	free(tree_ptr);	
}

#ifdef __cplusplus
}
#endif
