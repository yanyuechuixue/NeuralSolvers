#$Id: makefile,v 1.25 2014/04/06 10:47:26 zjcao Exp $

#include ../machines_data_base/for_compile/makefile.inc
include ../machines_data_base/for_compile/makefile_bnu

.SUFFIXES: .o .f90 .C .for

.f90.o:
	$(f90) $(f90appflags) -c $< -o $@

.C.o:
	${CXX} $(CXXAPPFLAGS) -c $< $(filein) -o $@

.for.o:
	$(f77) -c $< -o $@

# Input files

TwoPunctureFILES = TwoPunctureABE.o TwoPunctures.o 
AnsorgFILES = Ansorg.o
# file dependences

$(TwoPunctureFILES): TwoPunctures.h

# projects
TwoPunctureABE: $(TwoPunctureFILES)
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(TwoPunctureFILES) $(LDLIBS)

Ansorg: $(AnsorgFILES)
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(AnsorgFILES) $(LDLIBS)
clean:
	rm -f TwoPunctureABE $(TwoPunctureFILES)
	rm -f Ansorg $(AnsorgFILES)
