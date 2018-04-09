ROOT=../../..
include $(ROOT)/include/Makefile.inc

all: $(lykw)

$(lykw): $(patsubst %.C, $(OBJDIR)/extern/yukawa/%.o, $(wildcard *.C))
	@if [ ! -d $(LIBDIR) ]; then mkdir -p $(LIBDIR); fi
	$(AR) rcs $@ $^

.PHONY: all
