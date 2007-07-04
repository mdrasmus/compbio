#!/usr/bin/env python

import util
import sys
import ensembl
import genomeio

options = [
    ("d:", "desc=", "desc", "[ -d <desc>] [ --desc <desc> ]"),
    ("p", "part", "part", "[-p] [--part]")
    ]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)



exclude = {
 "receptor":1, "protein":1, "ec":1, "fragment":1,
 "family":1, "subfamily":1, "superfamily":1,
 "member":1, "putative":1, "probable": 1, 
 "precursor":1, "binding":1, "a":1, "homo":1, "sapiens":1, "factor":1, "dna":1,
 "specific":1, "subunit":1, "and":1, "of":1, "with":1, "to":1, "in":1, "like":1,
 "group":1, "associated":1, "gene":1, "homolog":1, "transcript":1, "contains":1,
 "chain":1, "domain":1, "includes":1, "containing":1, "type":1, "component":1,
 "acid":1, "amino":1, "related":1
}


desc = {}
if "desc" in param:
    for f in param["desc"]:
        desc.update(genomeio.readGeneDesc(f))



def summary(genes, desc):
    counts = {}
    for gene in genes:
        if gene in desc:
            words = {}
            for word in desc[gene].split():
                word = word.lower().replace("(", "").replace(")", "").\
                                    replace("[", "").replace("]", "").\
                                    replace(".", "").replace(",", "").\
                                    replace(";", "").replace(":", "").\
                                    replace("-", "")
                words[word] = 1
            for word in words:
                if word not in exclude and \
                   not word.isdigit() and \
                   word.find("source") == -1 and \
                   len(word) > 1:
                    counts[word] = counts.setdefault(word, 0) + 1
    return counts


wordsig = summary(desc.keys(), desc)

if "part" in param:
    lookup = {}
    
    for gene in desc:
        words = summary([gene], desc)
        keys = filter(lambda x: wordsig[x] > 50, words.keys())
        keys.sort(lambda a,b: cmp(wordsig[b], wordsig[a]))
        
        lookup.setdefault(tuple(keys), []).append(gene)
    
    for key in lookup:
        print key, lookup[key]
else:
    for gene in desc:
        words = summary([gene], desc)
        keys = filter(lambda x: wordsig[x] > 50, words.keys())
        keys.sort(lambda a,b: cmp(wordsig[b], wordsig[a]))

        if len(keys) > 0:
            print "%s\t%s" % (gene, "\t".join(keys[:5]))

