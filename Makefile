
.PHONY: bin cq test clean pkg pypi


bin:
	./make-bin.sh

pkg:
	python setup.py sdist

test:
	nosetests -v test test/rasmus test/compbio

cq:
	nosetests -v test/test_codequality.py

pypi:
	python setup.py register

clean:
	python setup.py clean
	rm -rf test/tmp
