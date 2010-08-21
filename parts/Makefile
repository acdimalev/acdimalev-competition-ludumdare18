editor: editor.c path.o
	gcc -o editor `sdl-config --cflags --libs` `pkg-config cairo --cflags --libs` editor.c path.o

%.o: %.c
	gcc -c `pkg-config cairo --cflags` $<
