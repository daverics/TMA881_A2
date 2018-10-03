CC = gcc
OPT = -O2
OBJ = main.o newton.o
FILE = main.c
GCOV = -pg -ftest-coverage -fprofile-arcs

.PHONY	:	all


all			:	$(OBJ)
	$(CC) -std=c11 -lm $(OBJ) -o main $(OPT)

main.o		:	main.c newton.c
newton.o	:	newton.c


clean:
	rm -f main main.gcno newton.o main.o
