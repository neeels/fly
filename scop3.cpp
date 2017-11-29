
#include <SDL2/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <cstdlib>
#include <limits.h>
#include <time.h>
#include <unistd.h>

#include <math.h>
#include <stdio.h>

#include <vector>
using namespace std;

#include "ctrl_layers.h"
#include "palettes.h"
#include "textures.h"

#define Pf(V) printf(#V "=%f\n", (float)V)

SDL_Window *window = NULL;

void draw_scene();

bool running = true;
volatile int frames_rendered = 0;
volatile int avg_frame_period = 0;
#define AVG_SHIFTING 3
float want_fps = 25;

float frandom(void) {
  return (float)(random()) / INT_MAX;
}

class Pt {
  public:

    double x;
    double y;
    double z;

    Pt() {
      set(0,0,0);
    }

    Pt(double x, double y, double z) {
      set(x, y, z);
    }

    void set(double x, double y, double z) {
      this->x = x;
      this->y = y;
      this->z = z;
    }

    void ang2cart(double r, double ang1, double ang2) {
      set(r * cos(ang1) * cos(ang2),
          r * cos(ang1) * sin(ang2),
          r * sin(ang1));
    }

    void ang2cart() {
      ang2cart(x, y, z);
    }

    void cart2ang() {
      // dunno
      double r = sqrt(x*x + y*y + z*z);
      if (fabs(r) < 1e-5) {
        x = y = z = 0;
      }
      else {
        double ay = acos(z / r);
        double az;
        if (fabs(x) < 1e-5) {
          az = 0;
        }
        else {
          az = atan(y / x);
        }
        x = r;
        y = ay;
        z = az;
      }
    }

    void set_min(Pt &p) {
      x = min(x, p.x);
      y = min(y, p.y);
      z = min(z, p.z);
    }

    void set_max(Pt &p) {
      x = max(x, p.x);
      y = max(y, p.y);
      z = max(z, p.z);
    }

    Pt& operator+=(const Pt &p) {
      x += p.x;
      y += p.y;
      z += p.z;
      return *this;
    }

    Pt& operator-=(const Pt &p) {
      x -= p.x;
      y -= p.y;
      z -= p.z;
      return *this;
    }

    Pt& operator*=(const double f) {
      x *= f;
      y *= f;
      z *= f;
      return *this;
    }


    void glVertex3d() {
      ::glVertex3d(x, y, z);
    }

    void glTranslated() {
      ::glTranslated(x, y, z);
    }

    void glRotated() {
      if (x) {
        ::glRotated(x,1,0,0);
      }
      if (y) {
        ::glRotated(y,0,1,0);
      }
      if (z) {
        ::glRotated(z,0,0,1);
      }
    }

    void glScaled() {
      ::glScaled(x, y, z);
    }

    void print() {
      printf("x%f y%f z%f\n", x, y, z);
    }

};




class Color {
  public:
    double r;
    double g;
    double b;
    double a;

    Color() {
      set(1,1,1,1);
    }

    void set(double r, double g, double b, double a=1) {
      this->r = r;
      this->g = g;
      this->b = b;
      this->a = a;
    }

    void set(rgb_t &rgb) {
      this->r = rgb.r;
      this->g = rgb.g;
      this->b = rgb.b;
    }

    void random(double cmin, double cmax) {
      r = frandom();
      g = frandom();
      b = frandom();
      a = 1;

      double mi = min(r, min(g, b));
      double ma = max(.1, max(r, max(g, b)) - mi);
      double intens = (cmax-cmin)*frandom();
      r = cmin + intens * ((r - mi) / ma);
      g = cmin + intens * ((g - mi) / ma);
      b = cmin + intens * ((b - mi) / ma);

    }

    void glColor(double alpha=1., double greying=0) {
      alpha *= a;
      if (greying > .01)
        ::glColor4f(( r * (1.-greying) + greying * .33),
                    ( g * (1.-greying) + greying * .33),
                    ( b * (1.-greying) + greying * .33),
                    ( alpha * (1.-greying) + greying * .33)
                    );
      else
        ::glColor4f(r, g, b, alpha);
    }
};



class Point : public Pt{
  public:

    Color c;

    Point() : Pt() {
      c.random(.3, 1);
    }

    Point(double x, double y, double z) : Pt(x, y, z) {
      c.random(.3, 1);
    }

    void draw(double alpha=1., double greying=0) {
      c.glColor(alpha, greying);
      glVertex3d();
    }
};

class Particle;

class DrawAs {
  public:
    double alpha;
    
    DrawAs() {
      alpha = 1;
    }

    void draw(vector<Point> &points);
};


class AsLines : DrawAs {
  public:
    AsLines() : DrawAs() {}
    void draw(vector<Point> &points)
    {
      int l;

      glBegin(GL_LINES);

      l = points.size();
      for (int i = 1; i < l; i++) {
        points[i-1].draw(alpha);
        points[i].draw(alpha);
      }

      glEnd();
    }
};

class AsTriangles : DrawAs {
  public:
    void draw(vector<Point> &points)
    {
      int l;
      glBegin(GL_TRIANGLES);

      l = points.size();
      for (int i = 2; i < l; i++) {
        points[i-2].draw(alpha);
        points[i-1].draw(alpha);
        points[i].draw(alpha);
      }

      glEnd();
    }
};

class AsTets : DrawAs {
  public:
    void draw(vector<Point> &points)
    {
      int l;
      glBegin(GL_TRIANGLES);

      l = points.size();
      for (int i = 3; i < l; i += 2) {
        points[i-3].draw(alpha);
        points[i-2].draw(alpha);
        points[i-1].draw(alpha);

        points[i-3].draw(alpha);
        points[i-2].draw(alpha);
        points[i].draw(alpha);

        points[i-3].draw(alpha);
        points[i-1].draw(alpha);
        points[i].draw(alpha);

        points[i-2].draw(alpha);
        points[i-1].draw(alpha);
        points[i].draw(alpha);
      }

      glEnd();
    }
};


class AsPoly : DrawAs {
  public:
    void draw(vector<Point> &points)
    {
      int l;
      glBegin(GL_POLYGON);

      l = points.size();
      for (int i = 0; i < l; i ++) {
        points[i].draw(alpha);
      }

      glEnd();
    }
};

class Placed {
  public:
    Pt pos;
    Pt dir;
    Pt scale;

    void gl_placement() {
      dir.glRotated();
      scale.glScaled();
      pos.glTranslated();
    }
};


class Particle : public Placed {
  public:
    vector<Point> points;

    Particle() {}

    void draw(DrawAs &as) { 
      glPushMatrix();
      gl_placement();

      as.draw(points);

      glPopMatrix();
    }
};


class Cloud : public Placed {
  public:
    vector<Particle> particles;
    vector<DrawAs> painters;

    Cloud()
    {}

    Particle &add_particle() {
      particles.resize( particles.size() + 1 );
      return particles.back();
    }

    void draw(DrawAs &as) { 
      glPushMatrix();
      gl_placement();

      for (int i = 0; i < particles.size(); i++) {
        particles[i].draw(as);
      }

      glPopMatrix();
    }

    void draw() {
      for (int i = 0; i < painters.size(); i++) {
        draw( painters[i] );
      }
    }
};



class Motion {
  public:
    bool done;

    Motion() {
      done = false;
    }

    void step();
};

class Dragon {
  public:


    vector<Point> segments;

    vector<Point> points1;
    vector<Point> points2;

    double maxangle;

    Dragon()
    {
      maxangle = M_PI * .9;
    }

    void update(double f_segments, double r, double fold, double alpha, double angle_zero,
                palette_t &pal)
    {
      Point p;

      int n_segments = (int)ceil(f_segments);
      n_segments = max(1, n_segments);

      if (segments.size() < n_segments) {
        for (int i = segments.size(); i < n_segments; i++) {
          p.set(frandom()*.9 + 5.1,
                (frandom() - .5)*2 * maxangle + (M_PI/3),
                (frandom() - .5)*2 * maxangle + (M_PI/3));
          p.c.random(.3, 1);
          segments.push_back( p );
        }
      }


#define do_points(_points, _r, _ang) \
      { \
        double angle_accum = 0; \
        p.set(0, 0, 0); \
        p.c.set( pal.colors[0] ); \
        _points.clear(); \
        _points.push_back(p); \
        for (int i = 0; i < n_segments; i += 2) { \
          Point q = segments[i]; \
          q.y *= fold; \
          q.z *= fold; \
          q._ang += angle_accum; \
          angle_accum += angle_zero; \
          q.ang2cart(); \
          q *= _r; \
          p += q; \
          p.c.set( pal.colors[ (int)((double)i / n_segments * pal.len) ] ); \
          p.c.a = alpha; \
          _points.push_back(p); \
        } \
      }


      do_points(points1, r, y);
      do_points(points2, -r, z);
      
    }

    void draw_lines(double alpha=1.) {
      int l;

      glBegin(GL_LINES);

      l = points1.size();
      for (int i = 1; i < l; i++) {
        points1[i-1].draw(alpha);
        points1[i].draw(alpha);
      }

      l = points2.size();
      for (int i = 1; i < l; i++) {
        points2[i-1].draw(alpha);
        points2[i].draw(alpha);
      }

      glEnd();
    }

    void draw_triangles(double alpha=1) {
      int l;
      glBegin(GL_TRIANGLES);

      l = points1.size();
      for (int i = 2; i < l; i++) {
        points1[i-2].draw(alpha);
        points1[i-1].draw(alpha);
        points1[i].draw(alpha);
      }

      l = points2.size();
      for (int i = 2; i < l; i++) {
        points2[i-2].draw(alpha);
        points2[i-1].draw(alpha);
        points2[i].draw(alpha);
      }
      glEnd();
    }

    void draw_polys(double alpha=1) {
      int l;
      glBegin(GL_POLYGON);

      l = points1.size();
      for (int i = 0; i < l; i++) {
        points1[i].draw(alpha);
      }

      l = points2.size();
      for (int i = 0; i < l; i++) {
        points2[i].draw(alpha);
      }
      glEnd();
    }
};


class Param {
  public:
    double val;
    double want_val;

    double slew;
    double change;

    bool do_limit_min;
    double val_min;

    bool do_limit_max;
    double val_max;


    Param(double val = 0, double slew=0.85){
      this->val = known_want_val = want_val = val;
      this->slew = slew;
      change = 0;
      do_limit_min = false;
      val_min = -1;
      do_limit_max = false;
      val_max = 1;
    }

    void step() {
      want_val += change;
      if (do_limit_min)
        want_val = max(val_min, want_val);
      if (do_limit_max)
        want_val = min(val_max, want_val);
      val = slew * val + (1. - slew) * want_val;
    }

    bool changed() {
      return known_want_val != want_val;
    }

    void set_unchanged() {
      known_want_val = want_val;
    }

    void limit(double min_val, double max_val) {
      val_min = min_val;
      do_limit_min = true;
      val_max = max_val;
      do_limit_max = true;
    }

    void limit_min(double min_val) {
      val_min = min_val;
      do_limit_min = true;
    }

    void limit_max(double max_val) {
      val_max = max_val;
      do_limit_max = true;
    }

    Param& operator=(double v) {
      want_val = v;
      return *this;
    }

    Param& operator=(Param &v) {
      want_val = v.want_val;
      return *this;
    }

    Param& operator+=(Param &v) {
      want_val += v.want_val;
      return *this;
    }

    Param& operator+=(double v) {
      want_val += v;
      return *this;
    }

    operator double() {
      return val;
    }


  private:
    double known_want_val;
};


class Animation {
  public:
    Param rot_x;
    Param rot_z;
    Param vision;
    Param alpha;
    Param dragon_points;
    Param dragon_r;
    Param dragon_fold;
    Param rotate_shift;
    Param distance;
    Param dist_scale;
    Param angle_zero;
    Param triangles_alpha;
    Param lines_alpha;
    Param lines_scale;
    Param palette_selected;
    Param palette_blend_speed;

    palette_t palette;
    palette_t blended_palette;
    palette_t *is_palette;
    palette_t *want_palette;
    double palette_blend;
    bool stop_palette_transition;

    vector<Dragon> dragons;

    Animation(int n_dragons = 3) {
      is_palette = &palettes[0];
      want_palette = is_palette;

      make_palettes();
      make_palette(&palette, PALETTE_LEN,
                   palette_defs[0]);

      make_palette(&blended_palette, PALETTE_LEN,
                   palette_defs[0]);


      dragon_points = 52;
      dragon_points.slew = 0;
      rot_x.slew = .94;
      rot_z.slew = .94;
      vision.slew = .96;
      alpha = .9;
      dragon_r = 1;
      dragon_r.limit_min(.01);
      dragon_fold = .01;
      rotate_shift = 0;
      rotate_shift.change = .0002;
      distance = 1.;
      distance.slew = .97;
      distance.limit_min(0);
      dist_scale = -.01;
      dist_scale.slew = .97;
      dist_scale.limit(-1, 1);
      triangles_alpha = .5;
      lines_alpha = 1;
      lines_scale = 1.1;
      palette_blend = 0;
      stop_palette_transition = false;
      palette_selected = 0;
      palette_selected.slew = 0;
      palette_blend_speed = .6;

      Dragon d;
      for (int i = 0; i < n_dragons; i++) {
        dragons.push_back(d);
        dragons.back().update(dragon_points, dragon_r, dragon_fold, alpha, angle_zero, palette);
      }
    }

    void step() {
      rot_x.step();
      rot_z.step();
      vision.step();
      alpha.step();
      dragon_points.step();
      dragon_r.step();
      dragon_fold.step();
      rotate_shift.step();
      distance.step();
      dist_scale.step();
      angle_zero.step();
      triangles_alpha.step();
      lines_alpha.step();
      lines_scale.step();
      palette_selected.step();
      palette_blend_speed.step();


      int _palette_selected = palette_selected;
      if ((_palette_selected >= 0)
          && (_palette_selected < n_palettes))
      {
        palette_t *pal_selected = &palettes[_palette_selected];
        if (want_palette != is_palette) 
        {
          if (pal_selected == want_palette)
          {
            stop_palette_transition = ! stop_palette_transition;
          }
          else
          if (pal_selected == is_palette) {
            stop_palette_transition = false;
            palette_t *tmp = is_palette;
            is_palette = want_palette;
            want_palette = tmp;
            palette_blend = 1. - palette_blend;
          }
          else {
            memcpy(blended_palette.colors, palette.colors,
                   blended_palette.len * sizeof(rgb_t));
            is_palette = &blended_palette;
            want_palette = pal_selected;
            stop_palette_transition = false;
            palette_blend = 0;
          }
        }
        else
        {
          want_palette = pal_selected;
          stop_palette_transition = false;
        }
      }
      if ((want_palette != is_palette) && ! stop_palette_transition) {
        float palette_blend_was = palette_blend;
        palette_blend += palette_blend_speed / want_fps;
        blend_palettes(&palette, is_palette, want_palette, min(1., palette_blend));
        if (palette_blend >= 1.) {
          is_palette = want_palette;
          palette_blend = 0;
        }
        if ((palette_blend_was < .5) != (palette_blend < .5))
          stop_palette_transition = true;
      }

      for (int i = 0; i < dragons.size(); i++) {
        dragons[i].update(dragon_points, dragon_r, dragon_fold, alpha, angle_zero, palette);
      }

      palette_selected = -1;
    }
};

Animation animation;



void draw_scene()
{
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity( );

  double vision = animation.vision;
  gluLookAt(3,4,260 - vision*200,0,0,vision * vision * 60,0,0,1);

  double sc = 1;

  double angleZ = animation.rot_z;
  double angleX = animation.rot_x;

  glScaled(sc, sc, sc);
  glRotated(angleZ,0,1,0);
  glRotated(angleX,1,0,0);

  double rotate_shift = animation.rotate_shift;
  int n_onion_shells = 100;
  for (int i = 0; i < n_onion_shells; i++) {
    double dd = i * (animation.distance + .4 * sin(rotate_shift * 1.23));
    dd /= sc;
    double d = dd * .33;

    glPushMatrix();
      glRotated(i * (16.789 + rotate_shift),1,sin(rotate_shift),.25 * cos(rotate_shift));

      double K = animation.dist_scale;
      double sc2 = 1. + K - 2.*K * i / n_onion_shells;
      glScaled(sc2,sc2,sc2);

      dd /= sc2;
      d /= sc2;

      double z1, z2;
      if (i & (int)1) {
        z1 = -d;
        z2 = dd;
      }
      else {
        z1 = dd;
        z2 = -d;
      }

      const int n_sphere_points = 4;
      Point sphere_points[n_sphere_points] = {
          Point(dd, 0., z1),
          Point(-d, d, z1),
          Point(-d, -d, z1),
          Point(0, 0, z2)
        };

      for (int j = 0; j < n_sphere_points; j++) {
        glPushMatrix();
        sphere_points[j].glTranslated();
        if (animation.triangles_alpha > .01)
          //animation.dragons[j % animation.dragons.size()].draw_triangles(animation.triangles_alpha);
          animation.dragons[j % animation.dragons.size()].draw_polys(animation.triangles_alpha);
        if (animation.lines_alpha > .01)
        {
          double lines_scale = animation.lines_scale;
          glScaled(lines_scale, lines_scale, lines_scale);
          animation.dragons[j % animation.dragons.size()].draw_lines(animation.lines_alpha);
        }
        glPopMatrix();
      }

    glPopMatrix();
  }

  glFlush();
  SDL_GL_SwapWindow(window);
}


FILE *out_stream = NULL;
FILE *out_params = NULL;
FILE *in_params = NULL;

char *audio_path = NULL;

int W = 1024;
int H = 768;

SDL_Surface *screen = NULL;

SDL_sem *please_save;
SDL_sem *saving_done;

int save_thread(void *arg) {

  for (;;) {
    SDL_SemWait(please_save);
    if (! running)
      break;

    if (out_stream) {
      fwrite(screen->pixels, sizeof(Uint32), W * H, out_stream);
    }
    SDL_SemPost(saving_done);
  }

  return 0;
}


void on_joy_axis(ControllerState &ctrl, int axis, double axis_val) {
  switch(ctrl.selected_layer) {
    default:
      break;

    case 0:
      switch(axis)
      {
        case 0:
          animation.dist_scale.change = axis_val / 20;
          break;
        case 1:
          animation.rotate_shift.change = .00003 + .01 * axis_val;
          break;
        default:
          break;
      }
      break;

    case 1:
      switch(axis)
      {
        case 1:
          animation.rot_x.change = (axis_val*axis_val*axis_val) * 6;
          break;
        case 0:
          animation.rot_z.change = (axis_val*axis_val*axis_val) * 6;
          break;
        case 7:
          break;
        case 3:
          {
            double dragon_points = (1. + axis_val)/2;
            dragon_points *= dragon_points;
            dragon_points = 2. + dragon_points * 24;
            animation.dragon_points = dragon_points;
          }
          break;
        case 4:
          animation.lines_scale.change = axis_val * .05;
          break;
        case 5:
          animation.vision = (1. + axis_val)/2;
          break;
        case 2:
          animation.alpha = 1. - (1. + axis_val)/2;
          break;
        default:
          break;
      }
      break;

    case 2:
      switch(axis)
      {
        case 1:
          animation.angle_zero.change = axis_val / 200;
          break;
        case 0:
          break;
        case 7:
          break;
        case 3:
          animation.lines_alpha.change = axis_val / 20;
          break;
        case 4:
          animation.triangles_alpha.change = -axis_val / 20;
          break;
        case 5:
          break;
        case 2:
          break;
        default:
          break;
      }
      break;

    case 3:
      switch(axis)
      {
        case 0:
          animation.dragon_fold.change = (axis_val*axis_val*axis_val) / 30;
          break;
        case 1:
          animation.dragon_r.change = -axis_val / 10;
          break;
        case 3:
          animation.rotate_shift.change = .0001 + .01 * axis_val;
          break;
        case 4:
          animation.distance.change = -(fabs(axis_val) < .05? 0. :
                                        .1 * (axis_val * axis_val * axis_val) );
          break;
        case 7:
          break;
        case 5:
          animation.dist_scale = .05 * (axis_val);
          break;
        case 2:
          break;
        default:
          break;
      }
      break;

  }

}

void on_joy_button(ControllerState &ctrl, int button, bool down) {
  switch(ctrl.selected_layer) {

    case 3:
      switch(button) {
        case 0:
          animation.angle_zero.change = 0;
          animation.angle_zero = 0;
          break;
        case 2:
          animation.dragon_fold.change = 0;
          animation.dragon_fold = 0;
          break;
      }
      break;

    case 2:
      if (down)
        animation.palette_selected = button;
      break;

  }
}



typedef struct {
  int random_seed;
  bool start_blank;
} init_params_t;

init_params_t ip;

int main(int argc, char *argv[])
{
  bool usage = false;
  bool error = false;
  bool fullscreen = false;

  int c;

  char *out_stream_path = NULL;
  char *out_params_path = NULL;
  char *in_params_path = NULL;

  ip.random_seed = time(NULL);

  while (1) {
    c = getopt(argc, argv, "hf:g:r:p:i:o:O:F");
    if (c == -1)
      break;
   
    switch (c) {
      case 'g':
        {
          char arg[strlen(optarg) + 1];
          strcpy(arg, optarg);
          char *ch = arg;
          while ((*ch) && ((*ch) != 'x')) ch ++;
          if ((*ch) == 'x') {
            *ch = 0;
            ch ++;
            W = atoi(arg);
            H = atoi(ch);

          }
          else {
            fprintf(stderr, "Invalid -g argument: '%s'\n", optarg);
            exit(-1);
          }
        }
        break;

      case 'F':
        fullscreen = true;
        break;

      case 'f':
        want_fps = atof(optarg);
        break;

      case 'r':
        ip.random_seed = atoi(optarg);
        break;

      case 'O':
        out_stream_path = optarg;
        break;

      case 'o':
        out_params_path = optarg;
        break;

      case 'i':
        in_params_path = optarg;
        break;

      case 'p':
        audio_path = optarg;
        break;

      case '?':
        error = true;
      case 'h':
        usage = true;
        break;

    }
  }

  if (usage) {
    if (error)
      printf("\n");
    printf(
"scop3 v0.1\n"
"(c) 2014 Neels Hofmeyr <neels@hofmeyr.de>\n"
"Published under the GNU General Public License v3.\n\n"
"Scop3 produces a mesmerizing animation controlled by any game controller.\n"
"\n"
"Usage example:\n"
"  scop3 -g 320x200 -f 25\n"
"\n"
"Options:\n"
"\n"
"  -g WxH   Set window width and height in number of pixels.\n"
"           Default is '-g %dx%d'.\n"
"  -F       Start in fullscreen mode.\n"
"  -f fps   Set desired framerate to <fps> frames per second. The framerate\n"
"           may slew if your system cannot calculate fast enough.\n"
"           If zero, run as fast as possible. Default is %.1f.\n"
"  -r seed  Supply a random seed to start off with.\n"
"  -O file  Write raw video data to file (grows large quickly). Can be\n"
"           converted to a video file using e.g. ffmpeg.\n"
"  -o file  Write live control parameters to file for later playback, see -i.\n"
"  -i file  Play back previous control parameters (possibly in a different\n"
"           resolution and streaming video to file...)\n"
"  -p file  Play back audio file in sync with actual framerate.\n"
"           The file format should match your sound card output format\n"
"           exactly.\n"
, W, H, want_fps
);
    if (error)
      return 1;
    return 0;
  }

  const int maxpixels = 1e4;

  if ((W < 3) || (W > maxpixels) || (H < 3) || (H > maxpixels)) {
    fprintf(stderr, "width and/or height out of bounds: %dx%d\n", W, H);
    exit(1);
  }

  SDL_Event event;

  SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_JOYSTICK);

  const int n_joysticks = SDL_NumJoysticks();
  controller_state_init(n_joysticks);

  SDL_Joystick **joysticks = NULL;

  if (n_joysticks) {
    SDL_JoystickEventState(SDL_ENABLE);

    joysticks = (SDL_Joystick**)malloc(sizeof(SDL_Joystick*) * n_joysticks);
    
    int i;
    for (i = 0; i < n_joysticks; i++)
    {
      printf("%2d: '%s'\n", i, SDL_JoystickNameForIndex(i));

      SDL_Joystick *j = SDL_JoystickOpen(i);
      printf("    %d buttons  %d axes  %d balls %d hats\n",
             SDL_JoystickNumButtons(j),
             SDL_JoystickNumAxes(j),
             SDL_JoystickNumBalls(j),
             SDL_JoystickNumHats(j)
             );
      joysticks[i] = j;
    }
  }

  atexit(SDL_Quit);
  window = SDL_CreateWindow("SCOP3",
                            SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                            W, H,
                            SDL_WINDOW_OPENGL);
  SDL_GLContext gl_ctx = SDL_GL_CreateContext(window);

  if (fullscreen)
    SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN);

  SDL_ShowCursor(SDL_DISABLE);

  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective(70,(double)W/H,1,10000);

  glEnable(GL_DEPTH_TEST);
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  draw_scene();

  Uint32 last_time = SDL_GetTicks();
  Uint32 current_time,elapsed_time;
  Uint32 start_time;

  if (out_stream_path) {
    if (access(out_stream_path, F_OK) == 0) {
      fprintf(stderr, "file exists, will not overwrite: %s\n", out_stream_path);
      exit(1);
    }
    out_stream = fopen(out_stream_path, "w");
  }

  /*
  double dragon_r = .3;
  double dragon_r_target = .3;
  double dragon_r_change = 0;
  double dragon_fold = 1;
  double dragon_points = 1;

  double rot_x = 0;
  double rot_z = 0;
  double want_vision = 0;
  double alpha = 0;
  double want_alpha = 0;
  */

  please_save = SDL_CreateSemaphore(0);
  saving_done = SDL_CreateSemaphore(1);

  SDL_Thread *save_thread_token = NULL;
  //if (out_stream)
  //  SDL_CreateThread(save_thread, NULL);

  float want_frame_period = (want_fps > .1? 1000. / want_fps : 0);
  float last_ticks = (float)SDL_GetTicks() - want_frame_period;

  char *pixelbuf = NULL;
  if (out_stream) {
    pixelbuf = (char*)malloc(W * H * 4); // 4 = RGBA
  }

  while (running)
  {
    animation.step();

    draw_scene();
    frames_rendered ++;

    {
      static int last_ticks2 = 0;

      int t = SDL_GetTicks();
      int elapsed = t - last_ticks2;
      last_ticks2 = t;

      avg_frame_period -= avg_frame_period >>AVG_SHIFTING;
      avg_frame_period += elapsed;
    }

    if (out_stream) {
      glReadPixels(0, 0, W, H, GL_RGBA, GL_UNSIGNED_BYTE, pixelbuf);
      fwrite(pixelbuf, sizeof(Uint32), W * H, out_stream);
    }

    while (running) {
      SDL_Event event;
      while (running && SDL_PollEvent(&event)) 
      {

        switch(event.type)
        {
          default:
            if (event.type == SDL_JOYBUTTONDOWN
                && event.jbutton.button == 7) {
              running = false;
              break;
            }
            handle_joystick_events(event);
            break;

          case SDL_QUIT:
            running = false;
            break;

          case SDL_KEYDOWN:
            {
              int c = event.key.keysym.sym;

              switch(c) {
                case SDLK_ESCAPE:
                  printf("Escape key. Stop.\n");
                  running = false;
                  break;

                case 13:
                  fullscreen = !fullscreen;
                  if (fullscreen)
                    SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN);
                  else
                    SDL_SetWindowFullscreen(window, 0);
                  break;
              }
            }
            break;
        }
      } // while sdl poll event

      if (want_frame_period) {
        int elapsed = SDL_GetTicks() - last_ticks;
        if (elapsed >= want_frame_period) {
          last_ticks += want_frame_period * (int)(elapsed / want_frame_period);
          break;
        }
        SDL_Delay(min((int)want_frame_period - elapsed, 5));
      }

    } // while running, for event polling / idle waiting

  } // while running

  running = false;

  printf("\n");
  printf("%d frames rendered\n", frames_rendered);
  if (out_stream) {
    fclose(out_stream);
    out_stream = NULL;

    printf("suggestion:\n"
        "ffmpeg -vcodec rawvideo -f rawvideo -pix_fmt rgb32 -s %dx%d -i %s ",
        W, H, out_stream_path);
    if (audio_path)
      printf("-i %s -acodec ac3 ", audio_path);
    printf("-vcodec mpeg4 -q 1 %s.%d.mp4\n", out_stream_path, H);
  }

  return 0;
}

