
.PHONY: bin cq test clean pkg


bin:
	./make-bin.sh

pkg:
	python setup.py sdist

test:
	nosetests -v test test/rasmus test/compbio

cq:
	nosetests -v test/test_codequality.py

clean:
	python setup.py clean
	rm -rf test/tmp
