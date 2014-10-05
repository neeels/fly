
default: t1 t2 scop3

scop3: scop3.cpp palettes.h palettes.o ctrl_layers.h
	g++ -g -o scop3 scop3.cpp palettes.o -lGL -lGLU -lSDL -lm
  
palettes.o: palettes.h palettes.c
	gcc -g -c -o palettes.o palettes.c

%: %.cpp
	g++ -g -o $* $< -lGL -lGLU -lSDL -lm
  
# vim: noexpandtab

