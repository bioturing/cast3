GCC_FLAG = 		-Wall -Wfatal-errors -O2 -std=c++11

OBJ_MAIN =		generate.o 		\
			genome.o

gen_read:

	g++ $(GCC_FLAG) -c -o generate.o generate.cpp
	g++ $(GCC_FLAG) -c -o genome.o genome.cpp
	g++ $(GCC_FLAG) $(OBJ_MAIN) main.cpp -o ../bin/gen_read
	rm -rf *.o