
.PHONY: bin cq test clean pypi


bin:
	./make-bin.sh

test:
	nosetests -v test test/rasmus test/compbio

cq:
	nosetests -v test/test_codequality.py

pypi:
        python setup.py register

clean:
	rm -rf test/tmp
