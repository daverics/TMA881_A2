CC = gcc
OPT = -O2
OBJ = main.c
FILE = main.c
GCOV = -pg -ftest-coverage -fprofile-arcs

.PHONY	:	all


all			:	$(OBJ)
	$(CC) -std=c11 -pthread -lm $(OBJ) -o main $(OPT)

# main.o		:	main.c newton.c
# newton.o	:	newton.c


clean:
	rm -f file main main.gcno newton.o main.o newton_attractors_xd.ppm newton_attractors_xd.txt
