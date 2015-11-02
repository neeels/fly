#ifndef TEXTURES_H
#define TEXTURES_H

#include <GL/gl.h>
#include <SDL/SDL.h>

#ifndef GL_CLAMP_TO_EDGE
#define GL_CLAMP_TO_EDGE 0x812F
#endif


GLuint load_texture(const char * filename,bool useMipMap);
SDL_Surface * flip_surface(SDL_Surface * surface);

void drawAxis(double scale);

#endif //TEXTURES_H
