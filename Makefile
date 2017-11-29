
default: scop3 fly

run: fly
	./run_fly.sh
	#./fly -g 1920x900

F = -g -std=c++11

scop3: scop3.cpp palettes.h palettes.o ctrl_layers.h
	g++ $F -o scop3 scop3.cpp palettes.o -lGL -lGLU -lSDL2 -lm
  
t: t.cpp palettes.h palettes.o ctrl_layers.h textures.o font.h draw.h foreach.h pt.h
	g++ $F -o t t.cpp palettes.o textures.o -lGL -lGLU -lSDL2 -lSDL2_image -lm

fly: fly.cpp palettes.h palettes.o ctrl_layers.h textures.o font.h draw.h foreach.h audio.h pt.h
	g++ $F -g -o fly fly.cpp palettes.o textures.o -lGL -lGLU -lglut -lSDL2 -lSDL2_image -lm 
  
palettes.o: palettes.h palettes.c
	gcc $F -c -o palettes.o palettes.c

%.o: %.cpp
	g++ $F -c $< 

%: %.cpp
	g++ $F -o $* $< -lGL -lGLU -lSDL2 -lm

litsphere: litsphere.c
	gcc -o litsphere litsphere.c -lGL -lGLU -lglut

clean:
	rm *.o

# vim: noexpandtab

