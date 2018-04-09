ROOT=../..
include $(ROOT)/include/Makefile.inc

all: $(lpot)

$(lpot): $(patsubst %.cpp, $(OBJDIR)/%.o, $(wildcard *.cpp))
	@if [ ! -d $(LIBDIR) ]; then mkdir -p $(LIBDIR); fi
	$(AR) rcs $@ $^

.PHONY: all
