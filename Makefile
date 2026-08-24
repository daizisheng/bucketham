# bucketham -- CWEB literate program
# `make`      -> tangle + compile the program
# `make pdf`  -> weave the documentation (bucketham.pdf)
# `make check`-> run the built-in self-tests

CC     = gcc
CFLAGS = -O2 -Wall

all: bucketham dualham

bucketham: bucketham.c
	$(CC) $(CFLAGS) -o $@ $<

dualham: dualham.c
	$(CC) $(CFLAGS) -o $@ $<

dualham.c: dualham.w
	ctangle dualham.w

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

.PHONY: all pdf check clean
