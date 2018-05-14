ROOT=..
include $(ROOT)/include/Makefile.inc

all:    anaHAL.x
sample: sample_effmass.x sample_pot_effSwave.x

anaHAL.x:              $(BINDIR)/anaHAL.x
sample_effmass.x:      $(BINDIR)/sample_effmass.x
sample_pot_effSwave.x: $(BINDIR)/sample_pot_effSwave.x

$(BINDIR)/anaHAL.x: ./main_anaHAL.cpp $(lpot) $(lfield) $(lana) $(lykw) $(lykwC)
	@if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR); fi
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) $(CXXLIBS) -o $@ $^

$(BINDIR)/sample_effmass.x: ./main_samples/effective_mass.cpp \
	$(lpot) $(lfield) $(lana) $(lykw) $(lykwC)
	@if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR); fi
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) $(CXXLIBS) -o $@ $^

$(BINDIR)/sample_pot_effSwave.x: ./main_samples/potential_effective_Swave.cpp \
	$(lpot) $(lfield) $(lana) $(lykw) $(lykwC)
	@if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR); fi
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) $(CXXLIBS) -o $@ $^

.PHONY: all sample anaHAL.x sample_effmass.x sample_pot_effSwave.x
