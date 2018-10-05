CC = gcc
OPT = -O2
OBJ = main.c
FILE = main
GCOV = -pg -ftest-coverage -fprofile-arcs

.PHONY	:	all


all			:	$(OBJ)
	$(CC) -std=c11 -pthread -lm $(OBJ) -o main $(OPT)

# main.o		:	main.c newton.c
# newton.o	:	newton.c


clean:
	rm -f main newton_attractors_xd.ppm

	# main.gcno newton.o main.o
