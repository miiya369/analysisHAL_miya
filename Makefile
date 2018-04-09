ROOT=.
include $(ROOT)/include/Makefile.inc

all:
	cd src;  make $@
	cd main; make $@

test:
	cd src;  make all
	cd main; make $@

miya:
	cd src;  make all
	cd main; make $@

lib:
	cd src;  make all

.PHONY: all test miya lib
