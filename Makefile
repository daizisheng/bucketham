# bucketham -- CWEB literate program
# `make`      -> tangle + compile the program
# `make pdf`  -> weave the documentation (bucketham.pdf)
# `make check`-> run the built-in self-tests

CC     = gcc
CFLAGS = -O2 -Wall

bucketham: bucketham.c
	$(CC) $(CFLAGS) -o $@ $<

bucketham.c: bucketham.w
	ctangle bucketham.w

pdf: bucketham.pdf
bucketham.pdf: bucketham.w
	cweave bucketham.w
	pdftex bucketham.tex

check: bucketham
	./bucketham

clean:
	rm -f bucketham bucketham.c bucketham.tex bucketham.pdf \
	      bucketham.idx bucketham.scn bucketham.toc bucketham.log

.PHONY: pdf check clean
