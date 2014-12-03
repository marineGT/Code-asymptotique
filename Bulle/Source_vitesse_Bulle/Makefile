# Makefile
shell = /bin/sh
FC = gfortran 
#FFLAGS = -fdefault-real-8 -O0 -g  -fbounds-check -pedantic-errors -ffpe-trap=invalid,zero,overflow,underflow -ftrapv -fimplicit-none -fno-automatic
NOM = vitesseBulle.exe
SRCS = 	calculsigma0zzB.f\
	CalculTabsB.f\
	constanteCM.f\
	constanteSn.f\
	diagonalessysteme1.f\
	DLegendre.f\
	dsigma0zB.f\
	dsigma0zxB.f\
	dsigma0zzB.f\
	du0zetaB.f\
	flegendre.f\
	fonctionb.f\
	Force0.f\
	Force1.f\
	funczeta1.f\
	gaussint1dim2.f\
	gaussint1dim3.f\
	gaussint1dim4.f\
	gaussint1dim5.f\
	gaussint1dim6.f\
	gaussint1dim7.f\
	gaussint1dim.f\
	gaussint1dimSn.f\
	gaussint1dimTest.f\
	ispline.f90\
	lecturedata.f\
	Legendre.f\
	Ntoptimal.f\
	ResolutionsystLU.f\
	secondmembreB.f\
	secondmembre.f\
	sigma0xxB.f\
	sigma0zzBCherZap.f\
	sigma0zzB.f\
	solparticuliere.f\
	solpartitest.f\
	splineA.f90\
	Sntot.f\
	surfacebulle.f\
	testcentremasse.f\
	TestCL_coordspherique.f\
	testvolume.f\
	thomas1d.f\
	u0etaB.f\
	verifsigma0zzB.f\
	vitesseBulle.f


OBJS = $(SRCS:.f=.o)

target1 : $(OBJS)
	  $(FC) $(FFLAGS) -o $(NOM) $(OBJS) -llapack -lblas
	mv $(NOM) ../

clean :
	rm $(OBJS)
