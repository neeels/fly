
#include <SDL/SDL.h>
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

#define Pf(V) printf(#V "=%f\n", (float)V)


void draw_scene();

bool running = true;
volatile int frames_rendered = 0;
volatile int avg_frame_period = 0;
#define AVG_SHIFTING 3
float want_fps = 25;

float frandom(void) {
  return (float)(random()) / INT_MAX;
}

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

    void glColor(double greying=0) {
      if (greying > .01)
        ::glColor4f(( r * (1.-greying) + greying * .33),
                    ( g * (1.-greying) + greying * .33),
                    ( b * (1.-greying) + greying * .33),
                    ( a * (1.-greying) + greying * .33)
                    );
      else
        ::glColor4f(r, g, b, a);
    }
};


class Point {
  public:

    double x;
    double y;
    double z;

    Color c;

    Point() {
      set(0,0,0);
      c.random(.3, 1);
    }

    Point(double x, double y, double z) {
      set(x, y, z);
      c.random(.3, 1);
    }

    void set(double x, double y, double z) {
      this->x = x;
      this->y = y;
      this->z = z;
    }

    void angular(double r, double ang1, double ang2) {
      set(r * cos(ang1) * cos(ang2),
          r * cos(ang1) * sin(ang2),
          r * sin(ang1));
    }

    void angular() {
      angular(x, y, z);
    }

    void set_min(Point &p) {
      x = min(x, p.x);
      y = min(y, p.y);
      z = min(z, p.z);
    }

    void set_max(Point &p) {
      x = max(x, p.x);
      y = max(y, p.y);
      z = max(z, p.z);
    }

    Point& operator+=(const Point &p) {
      x += p.x;
      y += p.y;
      z += p.z;
      return *this;
    }

    Point& operator-=(const Point &p) {
      x -= p.x;
      y -= p.y;
      z -= p.z;
      return *this;
    }

    Point& operator*=(const double f) {
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

    void print() {
      printf("x%f y%f z%f\n", x, y, z);
    }

    void draw(bool twice=false, double greying=0) {
      c.glColor(greying);
      glVertex3d();
      if (twice)
        glVertex3d();
    }
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

    void update(double f_segments, double r, double fold, double alpha)
    {
      Point p;
      Pf(f_segments);
      Pf(r);
      Pf(fold);
      Pf(alpha);

      int n_segments = (int)ceil(f_segments);
      n_segments = max(1, n_segments);

      if (segments.size() < n_segments) {
        printf("%d --> %d\n", segments.size(), n_segments);

        for (int i = segments.size(); i < n_segments; i++) {
          p.set(frandom()*.9 + 5.1,
                (frandom() - .5)*2 * maxangle + (M_PI/3),
                (frandom() - .5)*2 * maxangle + (M_PI/3));
          p.c.random(.3, 1);
          segments.push_back( p );
        }
      }


      p.set(0, 0, 0);
      p.c = segments[0].c;

      points1.clear();
      points1.push_back(p);

      for (int i = 0; i < n_segments; i += 2) {
        Point q = segments[i];
        q.y *= fold;
        q.z *= fold;
        q.angular();
        q *= r;
        p += q;
        p.c = q.c;
        double a = ((double)i)/f_segments;
        if (alpha < .5) 
          a *= (alpha*2);
        else
          a = max(a, (alpha - .5)*2);
        p.c.a = 1. - a;
        points1.push_back(p);
        p.print();
      }

      p.set(0, 0, 0);
      p.c = segments[0].c;

      points2.clear();
      points2.push_back(p);
      
      for (int i = 1; i < n_segments; i += 2) {
        Point q = segments[i];
        q.y *= fold;
        q.z *= fold;
        q.angular();
        q *= r;
        p -= q;
        p.c = q.c;
        double a = ((double)i)/f_segments;
        if (alpha < .5) 
          a *= (alpha*2);
        else
          a = max(a, (alpha - .5)*2);
        p.c.a = 1. - a;
        points2.push_back(p);
        p.print();
      }
      
    }



    void draw() {
      int l;

      glBegin(GL_LINES);

      points1[0].draw();
      l = points1.size() - 1;
      for (int i = 1; i < l; i++) {
        points1[i].draw(true);
      }
      points1[l].draw();

      points2[0].draw();
      l = points2.size() - 1;
      for (int i = 1; i < l; i++) {
        points2[i].draw(true);
      }
      points2[l].draw();
      glEnd();

      glBegin(GL_TRIANGLES);

      l = points1.size();
      for (int i = 2; i < l; i++) {
        points1[i-2].draw();
        points1[i-1].draw();
        points1[i].draw();
      }

      l = points2.size();
      for (int i = 2; i < l; i++) {
        points2[i-2].draw();
        points2[i-1].draw();
        points2[i].draw();
      }
      glEnd();
    }
};


class Param {
  public:
    double val;
    double want_val;

    double slew;

    Param(double val = 0, double slew=0.95){
      this->val = known_want_val = want_val = val;
      this->slew = slew;
    }

    void step() {
      printf("step %f", val);
      val = slew * val + (1. - slew) * want_val;
      printf(" --> %f", val);
    }

    bool changed() {
      return known_want_val != want_val;
    }

    void change_handled() {
      known_want_val = want_val;
    }

    Param& operator=(double v) {
      want_val = v;
      printf("assign want_val=%f\n", v);
      return *this;
    }

    Param& operator=(Param &v) {
      want_val = v.want_val;
      printf("assign want_val=Param(%f)\n", v.want_val);
      return *this;
    }

    Param& operator+=(Param &v) {
      want_val += v.want_val;
      return *this;
    }

    operator double() {
      printf("returning %f\n", val);
      return val;
    }


  private:
    double known_want_val;
};


class Animation {
  public:
    Param rot_x;
    Param rot_z;
    Param fold_speed;
    Param vision;
    Param alpha;
    Param dragon_points;
    Param dragon_r_change;
    Param dragon_r;
    Param dragon_fold;
    Param angleZ;
    Param angleX;

    vector<Dragon> dragons;

    Animation(int n_dragons = 3) {
      dragon_points = 3;
      rot_x.slew = .9;
      rot_z.slew = .9;
      fold_speed.slew = .9;
      vision.slew = .96;
      alpha = 1;
      dragon_r = .3;
      dragon_fold = 1;

      Dragon d;
      for (int i = 0; i < n_dragons; i++) {
        dragons.push_back(d);
        dragons.back().update(dragon_points, dragon_r, dragon_fold, alpha);
      }
    }

    void step() {
      rot_x.step();
      rot_z.step();
      fold_speed.step();
      vision.step();
      alpha.step();
      dragon_points.step();
      dragon_r_change.step();
      dragon_r.step();
      dragon_fold.step();
      angleZ.step();
      angleX.step();

      angleZ += rot_z;
      angleX += rot_x;
      dragon_fold += fold_speed;
      dragon_r += dragon_r_change;
      Pf(dragon_r);
      Pf(dragon_fold);
      Pf(alpha);
      for (int i = 0; i < dragons.size(); i++) {
        dragons[i].update(dragon_points, dragon_r, dragon_fold, alpha);
      }
    }
};

Animation animation;



void draw_scene()
{
  static double rotate_shift = 0;
  rotate_shift += .003;

  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity( );

  double vision = animation.vision;
  gluLookAt(3,4,260 - vision*200,0,0,vision * vision * 60,0,0,1);

  double sc = 1;

  double angleZ = animation.angleZ;
  double angleX = animation.angleX;

  glScaled(sc, sc, sc);
  glRotated(angleZ,0,1,0);
  glRotated(angleX,1,0,0);

  for (int i = 0; i < 100; i++) {
    double dd = i * (1 + .4 * sin(rotate_shift * 1.23));
    dd /= sc;
    double d = dd * .33;

    glPushMatrix();
      glRotated(i * (16.789 + rotate_shift),1,sin(rotate_shift),.25 * cos(rotate_shift));
      double sc2 = 1. - .01 * i;
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
        animation.dragons[j % animation.dragons.size()].draw();
        glPopMatrix();
      }

    glPopMatrix();
  }

  glFlush();
  SDL_GL_SwapBuffers();
}


FILE *out_stream = NULL;
FILE *out_params = NULL;
FILE *in_params = NULL;

char *audio_path = NULL;

int W = 1920;
int H = 1080;

SDL_Surface *screen = NULL;

SDL_sem *please_render;
SDL_sem *please_save;
SDL_sem *rendering_done;
SDL_sem *saving_done;

int render_thread(void *arg) {

  float want_frame_period = (want_fps > .1? 1000. / want_fps : 0);
  float last_ticks = (float)SDL_GetTicks() - want_frame_period;

  for (;;) {

    SDL_SemWait(please_render);
    if (! running)
      break;

    while (want_frame_period) {
      int elapsed = SDL_GetTicks() - last_ticks;
      if (elapsed >= want_frame_period) {
        last_ticks += want_frame_period * (int)(elapsed / want_frame_period);
        break;
      }
      SDL_Delay((int)want_frame_period - elapsed);
    }

    if (out_stream) {
      SDL_SemWait(saving_done);
    }

    draw_scene();

    int t = SDL_GetTicks();

    if (out_stream) {
      SDL_SemPost(please_save);
    }

    frames_rendered ++;
    SDL_SemPost(rendering_done);

    {
      static int last_ticks2 = 0;
      int elapsed = t - last_ticks2;
      last_ticks2 = t;
      avg_frame_period -= avg_frame_period >>AVG_SHIFTING;
      avg_frame_period += elapsed;
    }

  }

  return 0;
}

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
  switch(axis)
  {
    case 1:
      animation.rot_x = (axis_val*axis_val*axis_val) * 6;
      break;
    case 0:
      animation.rot_z = (axis_val*axis_val*axis_val) * 6;
      break;
    case 4:
      animation.fold_speed = (axis_val*axis_val*axis_val) / 30;
      break;
    case 7:
      animation.dragon_r_change = axis_val / 20;
      break;
    case 3:
      {
        double dragon_points = (1. + axis_val)/2;
        dragon_points *= dragon_points;
        dragon_points = 1. + dragon_points * 240;
        animation.dragon_points = dragon_points;
      }
      break;
    case 5:
      animation.vision = (1. + axis_val)/2;
      break;
    case 2:
      animation.alpha = (1. + axis_val)/2;
      break;
    default:
      break;
  }
}

void on_joy_button(ControllerState &ctrl, int button, bool down) {
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

  int c;

  char *out_stream_path = NULL;
  char *out_params_path = NULL;
  char *in_params_path = NULL;

  ip.random_seed = time(NULL);

  while (1) {
    c = getopt(argc, argv, "hf:r:p:i:o:O:");
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
      printf("%2d: '%s'\n", i, SDL_JoystickName(i));

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
  SDL_WM_SetCaption("SCOP3", NULL);
  int w = 1920;
  int h = 1080;
  screen = SDL_SetVideoMode(w,h, 32, SDL_OPENGL);
  SDL_ShowCursor(SDL_DISABLE);

  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective(70,(double)w/h,1,1000);

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
  double fold_speed = 0;
  double want_vision = 0;
  double alpha = 0;
  double want_alpha = 0;
  */

  please_render = SDL_CreateSemaphore(0);
  please_save = SDL_CreateSemaphore(0);
  rendering_done = SDL_CreateSemaphore(0);
  saving_done = SDL_CreateSemaphore(1);

  SDL_Thread *render_thread_token = SDL_CreateThread(render_thread, NULL);
  SDL_Thread *save_thread_token = NULL;
  if (out_stream)
    SDL_CreateThread(save_thread, NULL);

  while (running)
  {
    animation.step();

    SDL_SemPost(please_render);

    while (running) {
      SDL_Event event;
      while (running && SDL_PollEvent(&event)) 
      {

        switch(event.type)
        {
          default:
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
                    SDL_WM_ToggleFullScreen(screen);
                    break;
              }
            }
            break;
        }
      } // while sdl poll event

      if (SDL_SemTryWait(rendering_done) == 0)
        break;
      else
        SDL_Delay(5);
    } // while running, for event polling / idle waiting

  } // while running

  running = false;

  // make sure threads exit.
  SDL_SemPost(please_render);
  SDL_SemPost(please_render);
  if (out_stream) {
    SDL_SemPost(please_save);
    SDL_SemPost(please_save);
  }
  printf("waiting for render thread...\n");
  SDL_WaitThread(render_thread_token, NULL);
  if (out_stream) {
    printf("waiting for save thread...\n");
    SDL_WaitThread(save_thread_token, NULL);
  }

  SDL_DestroySemaphore(please_render);
  SDL_DestroySemaphore(please_save);
  SDL_DestroySemaphore(rendering_done);
  SDL_DestroySemaphore(saving_done);

  return 0;
}

