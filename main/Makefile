ROOT=..
include $(ROOT)/include/Makefile.inc

all:  anaHAL.x
test: test.x
miya: pot.x effpot.x ccp2.x ccp3.x tensor2.x makebin.x

anaHAL.x:  $(BINDIR)/anaHAL.x
test.x:    $(BINDIR)/test.x
pot.x:     $(BINDIR)/pot.x
effpot.x:  $(BINDIR)/effpot.x
ccp2.x:    $(BINDIR)/ccp2.x
ccp3.x:    $(BINDIR)/ccp3.x
tensor2.x: $(BINDIR)/tensor2.x
makebin.x: $(BINDIR)/makebin.x

$(BINDIR)/anaHAL.x: ./main_anaHAL.cpp $(lpot) $(lfield) $(lana) $(lykw) $(lykwC)
	@if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR); fi
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) $(CXXLIBS) -o $@ $^

$(BINDIR)/test.x: ./samples/test/main_test.cpp $(lpot) $(lfield) $(lana) $(lykw) $(lykwC)
	@if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR); fi
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) $(CXXLIBS) -o $@ $^

$(BINDIR)/pot.x: ./samples/potential/main_leading_order_potentials_octet.cpp \
		 $(lpot) $(lfield) $(lana) $(lykw) $(lykwC)
	@if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR); fi
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) $(CXXLIBS) -o $@ $^

$(BINDIR)/effpot.x: ./samples/potential/main_eff_Swave_potentials_octet.cpp \
		    $(lpot) $(lfield) $(lana) $(lykw) $(lykwC)
	@if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR); fi
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) $(CXXLIBS) -o $@ $^

$(BINDIR)/ccp2.x: ./samples/potential/main_eff_Swave_CCP_LxN_SxN.cpp \
		  $(lpot) $(lfield) $(lana) $(lykw) $(lykwC)
	@if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR); fi
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) $(CXXLIBS) -o $@ $^

$(BINDIR)/ccp3.x: ./samples/potential/main_eff_Swave_CCP_LxN_SxN_SxSN.cpp \
		  $(lpot) $(lfield) $(lana) $(lykw) $(lykwC)
	@if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR); fi
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) $(CXXLIBS) -o $@ $^

$(BINDIR)/tensor2.x: ./samples/potential/main_tensor_CCP_LxN_SxN.cpp \
		  $(lpot) $(lfield) $(lana) $(lykw) $(lykwC)
	@if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR); fi
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) $(CXXLIBS) -o $@ $^

$(BINDIR)/makebin.x: ./samples/wave/main_make_bin_wave_LxN_SxN_SxSN.cpp \
		  $(lfield) $(lana) $(lykw) $(lykwC)
	@if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR); fi
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) $(CXXLIBS) -o $@ $^

.PHONY: all test miya anaHAL.x test.x pot.x effpot.x ccp2.x ccp3.x tensor2.x makebin.x
