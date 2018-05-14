ROOT=.
include $(ROOT)/include/Makefile.inc

all:
	cd src;  make $@
	cd main; make $@

sample:
	cd src;  make all
	cd main; make $@

lib:
	cd src;  make all

.PHONY: all sample lib
