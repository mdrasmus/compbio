
from compbio import alignlib
from compbio import fasta


def test_align_subsets():
    aln = fasta.FastaDict()
    aln["a"] = "AAA-A"
    aln["b"] = "-BD-C"
    aln["c"] = "A-D--"

    aln2 = alignlib.remove_empty_columns(aln)
    assert aln2 == {'a': 'AAAA', 'c': 'A-D-', 'b': '-BDC'}

    aln2 = alignlib.remove_gapped_columns(aln)
    assert aln2 == {'a': 'A', 'c': 'D', 'b': 'D'}

    aln2 = alignlib.require_nseqs(aln, 2)
    assert aln2 == {'a': 'AAAA', 'c': 'A-D-', 'b': '-BDC'}
