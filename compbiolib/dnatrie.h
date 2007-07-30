#ifndef DNATRIE_H
#define DNATRIE_H


#ifdef __cplusplus
extern "C" {
#endif

/* Trie datastructure */
typedef struct DnaTrie_struct
{
    unsigned char strand;
    struct DnaTrie_struct *next[4];
} DnaTrie;


void dnatrie_insert(DnaTrie** tree_ptr, char* seq, unsigned char strand);
int dnatrie_query(DnaTrie *node, char *seq);
void dnatrie_delete(DnaTrie* tree_ptr);

#ifdef __cplusplus
} // extern C
#endif

#endif /* DNA_TRIE_H */
