
default: t1 t2 scop3 t fly

F = -g -std=c++11

scop3: scop3.cpp palettes.h palettes.o ctrl_layers.h
	g++ $F -o scop3 scop3.cpp palettes.o -lGL -lGLU -lSDL -lm
  
t: t.cpp palettes.h palettes.o ctrl_layers.h textures.o font.h draw.h foreach.h
	g++ $F -o t t.cpp palettes.o textures.o -lGL -lGLU -lSDL -lSDL_image -lm
  
fly: fly.cpp palettes.h palettes.o ctrl_layers.h textures.o font.h draw.h foreach.h
	g++ $F -g -o fly fly.cpp palettes.o textures.o -lGL -lGLU -lSDL -lSDL_image -lm
  
palettes.o: palettes.h palettes.c
	gcc $F -c -o palettes.o palettes.c

%: %.cpp
	g++ $F -o $* $< -lGL -lGLU -lSDL -lm
  
# vim: noexpandtab

