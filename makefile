//ROOTFLAGS = `root-config --cflags`
//LIBROOT = `root-config --libs`

Potentials.o: Potentials.cpp
	g++ -c Potentials.cpp $(ROOTFLAGS) $(LIBROOT)

Schroddy.o: Schroddy.cpp
	g++ -c Schroddy.cpp $(ROOTFLAGS) $(LIBROOT)

main.o: main.cpp
	g++ -c main.cpp $(ROOTFLAGS) $(LIBROOT)

Build: Potentials.o Schroddy.o main.o
	g++ Potentials.o Schroddy.o main.o -o exec.x $(ROOTFLAGS) $(LIBROOT)

Execute: Execute exec.x
	./exec.x

