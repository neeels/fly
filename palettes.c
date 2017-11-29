#include <stdlib.h>

#include "palettes.h"

#define min(A,B) ((A) > (B)? (B) : (A))
#define max(A,B) ((A) > (B)? (A) : (B))

palette_def_t pal1 = 
    {
      11,
      {
        { 0./6, 1, 1, 1 },
        { 0.5/6, 0.1, 1, .10 },
        { 1./6, 0, .1, .5 },
        { 1.5/6, .3, .3,1 },
        { 3./6, 1, 1, 0 },
        { 3.5/6, 0, 1, 1 },
        { 4.5/6, 0, 0, 1 },
        { 4.8/6, 1, 0, 0 },
        { 5.25/6, 1, 1, 0 },
        { 5.55/6, 1, .6, 0 },
        { 5.85/6, .1,.0,0 },
      }
    };

palette_def_t pal2 =
    {
      2,
      {
        { 0, 1./255*0x3e, 1./255*0x86, 1./255*0xd7 },
        { 0.5 + 3./256, 1./255*0xd0, 1./255*0x29, 1./255*0x40 },
      }
    };

#define r1 (1./256 * 0xf8)
#define g1 (1./256 * 0xb8)
#define b1 (1./256 * 0x0c)
#define r9 (r1 *.98)
#define g9 (g1 *.98)
#define b9 (b1 *.98)
#define r0 .0
#define g0 .0
#define b0 .0
#define ri (r1 * .02)
#define gi (g1 * .02)
#define bi (b1 * .02)
#define d -.1/6

palette_def_t pal3 =
    {
      40,
      {
        {0.+0./6,   r0, gi, b0 },
     {d+ 0.+.5/6,   r0, gi, b0 },
        {0.+.5/6,   r9, g0, b1 },
     {d+ 0.+1./6,   r9, g0, b1 },
        {0.+1./6,   r0, g0, b0 },
     {d+ 0.+1.5/6,  r0, g0, b0 },
        {0.+1.5/6,  r9, g9*.5, b1 },
     {d+ 0.+3./6,   r9, g9*.6, b1 },
        {0.+3./6,   ri, g0, bi },
     {d+ 0.+3.5/6,  ri, g0, bi },
        {0.+3.5/6,  r1, gi*.1, b9 },
     {d+ 0.+4.5/6,  r1, gi*.1, b9 },
        {0.+4.5/6,  ri, gi, bi },
     {d+ 0.+4.8/6,  ri, gi, bi },
        {0.+4.8/6,  r1, g1*.6, b9 },
     {d+ 0.+5.25/6, r1, g1*.7, b9 },
        {0.+5.25/6, r0, gi, bi },
     {d+ 0.+5.7/6,  r0, gi, bi },
        {0.+5.7/6,  r9, g1*.2, b9 },
     {d+ 1.+0./6,   r9, g1*.2, b9 },
        {1.+0./6,   r0, gi, b0 },
     {d+ 1.+.5/6,   r0, gi, b0 },
        {1.+.5/6,   r9, g1*.7, b9 },
     {d+ 1.+1./6,   r9, g1*.8, b9 },
        {1.+1./6,   r0, g0, b0 },
     {d+ 1.+1.5/6,  r0, g0, b0 },
        {1.+1.5/6,  r9, g9*.3, b9 },
     {d+ 1.+3./6,   r9, g9*.3, b9 },
        {1.+3./6,   ri, g0, b0 },
     {d+ 1.+3.5/6,  ri, g0, b0 },
        {1.+3.5/6,  r1, g9*.9, b1 },
     {d+ 1.+4.5/6,  r1, g9*.8, b1 },
        {1.+4.5/6,  ri, gi, bi },
     {d+ 1.+4.8/6,  ri, gi, bi },
        {1.+4.8/6,  r1, g1*.4, b1 },
     {d+ 1.+5.25/6, r1, g1*.4, b1 },
        {1.+5.25/6, r0, gi, b0 },
     {d+ 1.+5.7/6,  r0, gi, b0 },
        {1.+5.7/6,  r9, g1, b1 },
     {d+ 2,         r0, gi, b0 },
      }
    };
#undef O
#undef g
#undef I
#undef i
#undef d

#define d .1
palette_def_t pal4 = 
    {
      4,
      {
        {0, 0, 0, 0},
        {d, 1, 1, 1},
        {0.5, 1, 1, 1},
        {0.5+d, 0, 0, 0},
      }
    };

palette_def_t pal5 =
    {
      2,
      {
        { 0, .3, .3, .3 },
        { 0.5 + 3./256, .3, .3, .3 },
      }
    };

#undef d

#define s .2
#define d .05


palette_def_t pal6 = 
    {
      6,
      {
        {s+d, 0, 0, 0},
        {s+d*2, 1, 1, 1},
        {s+d*3, 0, 0, 0},
        {s+d*4, 0, 0, 0},
        {s+d*5, 1, 1, 1},
        {s+d*6, 0, 0, 0},
      }
    };


#undef d

palette_def_t *palette_defs[n_palettes] = {
    &pal1,
    &pal2,
    &pal3,
    &pal4,
    &pal5,
    &pal6,
  };

palette_t palettes[n_palettes];


void set_color(palette_t *palette, int i, float r, float g, float b) {
  if (i >= palette->len)
    return;
  palette->colors[i].r = r;
  palette->colors[i].g = g;
  palette->colors[i].b = b;
}

void make_palette(palette_t *palette, int n_colors,
                  palette_def_t *palette_def) {
  int i;
  int n_points = palette_def->n_points;
  palette_point_t *points = palette_def->points;

  palette->colors = malloc(n_colors * sizeof(rgb_t));
  palette->len = n_colors;


  if (n_points < 1) {
    for (i = 0; i < palette->len; i++) {
      float val = (float)i / palette->len;
      set_color(palette, i, val, val, val);
    }
    return;
  }

  palette_point_t *last_p = points;
  palette_point_t *first_p = points;

  for (i = 1; i < n_points; i ++) {
    if (points[i].pos > last_p->pos)
      last_p = &points[i];
    if (points[i].pos < first_p->pos)
      first_p = &points[i];
  }
  if (last_p->pos > 1.0) {
    float norm_factor = last_p->pos * (n_points + 1) / n_points;
    for (i = 0; i < n_points; i ++)
      points[i].pos /= norm_factor;
  }
  
  // duplicate the last point to "the left", wrap back below zero.
  palette_point_t p = *last_p;
  p.pos -= 1.0;
  // ...unless another point is defined there.
  if (p.pos >= first_p->pos)
    p = *first_p;

  // also duplicate the first point to "the right".
  palette_point_t post_last = *first_p;
  post_last.pos += 1.0;

  int color_pos = 0;

  while(color_pos < n_colors) {

    // look for the next point, the one with the next largest pos after p.pos
    palette_point_t *next_p = NULL;
    for (i = 0; i < n_points; i ++) {
      float i_pos = points[i].pos;
      if ((i_pos > p.pos)
          &&
          (
           (! next_p)
           || (i_pos < next_p->pos)
          )
         )
        next_p = &points[i];
    }

    if (! next_p)
      next_p = &post_last;

    int next_color_pos = (int)(next_p->pos * n_colors) + 1;

    if (next_color_pos <= color_pos)
      next_color_pos = color_pos + 1;

    for (; color_pos < next_color_pos; color_pos ++) {
      float prevpos = p.pos;
      float nextpos = next_p->pos;
      float currentpos = ((float)color_pos) / n_colors;
      float fade;
      if ((nextpos - prevpos) < 1e-3)
        fade = 0.5;
      else
        fade = (currentpos - prevpos) / (nextpos - prevpos);
      float rfade = 1.0 - fade;
      float r = rfade * p.r  +  fade * next_p->r;
      float g = rfade * p.g  +  fade * next_p->g;
      float b = rfade * p.b  +  fade * next_p->b;

      set_color(palette, color_pos, r, g, b);
    }

    p = *next_p;
  }
}

void make_palettes() {
  int pal_i;
  for (pal_i = 0; pal_i < n_palettes; pal_i ++) {
    make_palette(&palettes[pal_i], PALETTE_LEN,
                 palette_defs[pal_i]);
  }
}

void blend_palettes(palette_t *dst, palette_t *src1, palette_t *src2, float blend) {
  assert((dst->len == src1->len) && (dst->len == src2->len));

  float ba = min(1., (1.-blend) * 2);
  float bb = min(1., blend * 2);
  int i, j;
  for (i = 0; i < dst->len; i++) {
    float *a = (float*)(&src1->colors[i]);
    float *b = (float*)(&src2->colors[i]);
    float *x = (float*)(&dst->colors[i]);
    for (j = 0; j < 3; j++) {
      x[j] = max(ba * a[j], bb * b[j]);
    }
  }
}

