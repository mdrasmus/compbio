
.PHONY: bin cq test


bin:
	./make-bin.sh

test:
	nosetests -v test test/rasmus test/compbio

cq:
	nosetests -v test/test_codequality.py
