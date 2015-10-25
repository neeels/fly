
default: t1 t2 scop3 t fly

scop3: scop3.cpp palettes.h palettes.o ctrl_layers.h
	g++ -g -o scop3 scop3.cpp palettes.o -lGL -lGLU -lSDL -lm
  
t: t.cpp palettes.h palettes.o ctrl_layers.h textures.o font.h draw.h
	g++ -g -o t t.cpp palettes.o textures.o -lGL -lGLU -lSDL -lSDL_image -lm
  
fly: fly.cpp palettes.h palettes.o ctrl_layers.h textures.o font.h draw.h
	g++ -g -o fly fly.cpp palettes.o textures.o -lGL -lGLU -lSDL -lSDL_image -lm
  
palettes.o: palettes.h palettes.c
	gcc -g -c -o palettes.o palettes.c

%: %.cpp
	g++ -g -o $* $< -lGL -lGLU -lSDL -lm
  
# vim: noexpandtab

