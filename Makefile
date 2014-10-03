
default: t1 t2 scop3

%: %.cpp
	g++ -g -o $* $< -lGL -lGLU -lSDL -lm
  
# vim: noexpandtab

