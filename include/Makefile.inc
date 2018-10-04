
### You can change here for your environment ###

CXX = g++-5
AR  = ar

OPTION_CXXFLAGS    = -fopenmp #-DDEBUG_MODE
OPTION_CXXINCLUDES = 
OPTION_CXXLIBS     = 

#===================================================================================#
#============================= DON’T CHANGE FROM HERE ==============================#
#===================================================================================#

OBJDIR = $(ROOT)/obj
LIBDIR = $(ROOT)/lib
BINDIR = $(ROOT)/bin

lana   = $(LIBDIR)/libanalysis.a
lfield = $(LIBDIR)/libfield.a
lpot   = $(LIBDIR)/libpotential.a
lykw   = $(LIBDIR)/libyukawa.a
lykwC  = $(LIBDIR)/libyukawaCORE.a

CXXFLAGS    = $(OPTION_CXXFLAGS)    -O3 -g -Wall
CXXINCLUDES = $(OPTION_CXXINCLUDES) -I$(ROOT)/include -I$(ROOT)/include/extern
CXXLIBS     = $(OPTION_CXXLIBS)     -L$(LIBDIR)

$(OBJDIR)/%.o: %.cpp
	@if [ ! -d $(OBJDIR) ]; then mkdir -p $(OBJDIR); fi
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) $(CXXLIBS) -o $@ -c $<

$(OBJDIR)/extern/yukawa/%.o: %.C
	@if [ ! -d $(OBJDIR)/extern/yukawa ]; then mkdir -p $(OBJDIR)/extern/yukawa; fi
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) $(CXXLIBS) -o $@ -c $<

Default: all

clean:
	rm -rf $(OBJDIR) $(LIBDIR) $(BINDIR)
	find . -name  "*~"       -exec rm -f {} \;
	find . -name ".*~"       -exec rm -f {} \;
	find . -name "._*"       -exec rm -f {} \;
	find . -name ".DS_Store" -exec rm -f {} \;

.PHONY: clean Default
