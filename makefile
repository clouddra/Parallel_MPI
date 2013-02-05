CC=mpicc


all: Training Match


Match: Match.c
	$(CC) Match.c -o Match -lrt -lm


Training:	Training.c
	$(CC) Training.c -o Training

clean:
	rm -rf Training Match
