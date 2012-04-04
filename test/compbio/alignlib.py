
from compbio import fasta, alignlib


#=============================================================================
# testing

if __name__ == "__main__":
    aln = fasta.FastaDict()
    aln["a"] = "AAA-A"
    aln["b"] = "-BD-C"
    aln["c"] = "A-D--"

    alignlib.print_align(alignlib.remove_empty_columns(aln))
    alignlib.print_align(alignlib.remove_gapped_columns(aln))
    alignlib.print_align(alignlib.require_nseq(aln, 2))
