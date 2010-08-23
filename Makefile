untitled: untitled.c physics.o
	gcc -o untitled -lm `sdl-config --cflags --libs` `pkg-config --cflags --libs cairo` untitled.c physics.o

physics.o: physics.c physics.h
	gcc -c physics.c
