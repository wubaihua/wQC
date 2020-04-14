FC = ifort

wqc:  constant.o math.o cint.o def.o geo.o HF.o PT.o fileio.o wQC.o
	${FC} -o wqc *.f90 -mkl -lcint

constant.o: constant.f90
	${FC} -c constant.f90

constant.mod: constant.f90 constant.o
	${FC} -c constant.f90

math.o: math.f90
	${FC} -c math.f90 -mkl

math.mod: math.f90 math.o
	${FC} -c math.f90 -mkl

cint.o: cint.f90
	${FC} -c cint.f90 -lcint

cint.mod: cint.f90 cint.o
	${FC} -c cint.f90 -lcint

def.o: def.f90
	${FC} -c def.f90 

def.mod: def.f90 def.o
	${FC} -c def.f90 

geo.o: geo.f90 def.mod math.o
	${FC} -c geo.f90 -mkl #def.mod math.o -mkl

geo.mod: geo.f90 geo.o def.mod math.o
	${FC} -c geo.f90 -mkl #def.mod math.o -mkl

HF.o: HF.f90 geo.mod math.mod def.mod
	${FC} -c HF.f90 -mkl #geo.mod math.mod def.mod -mkl

PT.o: PT.f90
	${FC} -c PT.f90

fileio.o: fileio.f90 def.mod
	${FC} -c fileio.f90 #def.mod

wQC.o: wQC.f90 def.mod cint.mod
	${FC} -c wQC.f90 #def.mod cint.mod



