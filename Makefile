FC = ifort

wqc:  constant.o math.o cint.o def.o geo.o init.o HF.o PT.o fileio.o wQC.o
	${FC} -o wqc *.f90 -mkl -lcint #-g -check all -fpe0 -warn -traceback -debug extended -check bound

constant.o: constant.f90
	${FC} -c constant.f90

constant.mod: constant.f90 constant.o
	${FC} -c constant.f90

math.o: math.f90
	${FC} -c math.f90 -mkl

math.mod: math.f90 math.o
	${FC} -c math.f90 -mkl

def.o: def.f90
	${FC} -c def.f90 

def.mod: def.f90 def.o
	${FC} -c def.f90 

geo.o: geo.f90 def.mod math.o
	${FC} -c geo.f90 -mkl #def.mod math.o -mkl

geo.mod: geo.f90 geo.o def.mod math.o
	${FC} -c geo.f90 -mkl #def.mod math.o -mkl
	
init.o: init.f90 def.mod
	${FC} -c init.f90 
	
init.mod: init.f90 init.o def.mod
	${FC} -c init.f90 

cint.o: cint.f90 init.mod
	${FC} -c cint.f90 -lcint

cint.mod: cint.f90 cint.o init.mod
	${FC} -c cint.f90 -lcint

HF.o: HF.f90 geo.mod math.mod def.mod
	${FC} -c HF.f90 -mkl #geo.mod math.mod def.mod -mkl 

pop.o: pop.f90 def.mod math.o
	${FC} -c pop.f90 -mkl

PT.o: PT.f90
	${FC} -c PT.f90

fileio.o: fileio.f90 def.mod
	${FC} -c fileio.f90 #def.mod
	
order.o: order.f90
	${FC} -c order.f90

wQC.o: wQC.f90 def.mod cint.mod init.mod
	${FC} -c wQC.f90 #def.mod cint.mod



clean: 
	rm *.o *.mod wqc



