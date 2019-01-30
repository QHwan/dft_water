qdft : analysis.o math2.o ideal.o hard.o disp.o input.o program.o main.o
	gcc -o ./qdft -O3 -lm analysis.o math2.o ideal.o hard.o disp.o input.o program.o main.o 

analysis.o : analysis.c
	gcc -c -O3 analysis.c

program.o : program.c
	gcc -c -O3 program.c

input.o : input.c
	gcc -c -O3 input.c

ideal.o : ideal.c
	gcc -c -O3 ideal.c

hard.o : hard.c 
	gcc -c  -O3 hard.c 

disp.o : disp.c
	gcc -c -O3 disp.c

math2.o : math2.c
	gcc -c -O3 math2.c 

main.o : main.c
	gcc -c  -O3 main.c


