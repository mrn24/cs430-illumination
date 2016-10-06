all: illumination.c
	gcc illumination.c -o raycast -lm

clean:
	rm -rf raycast *~
