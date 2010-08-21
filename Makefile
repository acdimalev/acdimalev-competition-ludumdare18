untitled: untitled.c
	gcc -o untitled -lm `sdl-config --cflags --libs` `pkg-config --cflags --libs cairo` untitled.c
