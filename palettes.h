#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>

#define PALETTE_LEN_BITS 12
#define PALETTE_LEN (1 << PALETTE_LEN_BITS)
#define PALETTE_LEN_MASK (PALETTE_LEN - 1)

typedef struct{
  float r;
  float g;
  float b;
} rgb_t;

typedef struct {
  rgb_t *colors;
  unsigned int len;
} palette_t;

typedef struct {
  float pos;
  float r;
  float g;
  float b;
} palette_point_t;

typedef struct {
  unsigned int n_points;
  palette_point_t points[];
} palette_def_t;


#define n_palettes 6
extern palette_t palettes[n_palettes];
extern palette_def_t *palette_defs[n_palettes];

extern void make_palettes();
extern void blend_palettes(palette_t *dst, palette_t *src1, palette_t *src2, float blend);

/* Generate a color palette, setting palette->colors and palette->len.
 * Allocate new memory for palette->colors (is not freed or reallocd).
 * 'n_colors' defines how many colors are generated in the palette.  'points'
 * is a definition of colors at specific intervals, 'n_points' gives the number
 * of palette_point_t array elements in 'points'.
 */
extern void make_palette(palette_t *palette, int n_colors,
                  palette_def_t *palette_def);

#ifdef __cplusplus
}
#endif
