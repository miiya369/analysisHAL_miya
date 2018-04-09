ROOT=../../..
include $(ROOT)/include/Makefile.inc

all: $(lykwC)

$(lykwC): $(patsubst %.C, $(OBJDIR)/extern/yukawa/%.o, $(wildcard *.C))
	@if [ ! -d $(LIBDIR) ]; then mkdir -p $(LIBDIR); fi
	$(AR) rcs $@ $^

.PHONY: all
