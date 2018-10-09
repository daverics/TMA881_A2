CC = gcc
OPT = -O2
OBJ = main.o
FILE = main.c
GCOV = -pg -ftest-coverage -fprofile-arcs

.PHONY	:	all


all			:	$(OBJ)
	$(CC) -std=c11 -pthread -lm $(GCOV) $(OBJ) -o main $(OPT)

main.o		:	main.c
#newton.o	:	newton.c


clean:
	rm -f file main main.gcno newton.o main.o newton_attractors_x7.ppm newton_convergence_xd.ppm
