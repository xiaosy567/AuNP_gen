#All_Objects=atom.o bond.o angle.o residue.o function.o
CFlag=-O2
#all: ${All_Objects}
all: sys.o
	g++ ${CFlag} -o a.out sys.o
	rm *.o
sys.o:sys.cpp
	g++ ${CFlag} -c sys.cpp
clean :
	rm *.out
