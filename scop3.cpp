
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <cstdlib>
#include <limits.h>

#include <math.h>
#include <stdio.h>

#include <vector>
using namespace std;

#include "ctrl_layers.h"

void Dessiner();

bool running = true;
double angleZ = 0;
double angleX = 0;

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
    };

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


Dragon dragon1;
Dragon dragon2;
Dragon dragon3;

double vision = 0;

class Param {
  public:
    double val;
    double want_val;

    double slew;

    Param(double val = 0, double slew=0.5){
      this->val = known_want_val = want_val = val;
      this->slew = slew;
    }

    Param& operator=(double v) {
      if (v != want_val) {
        want_val = v;
      }
    }

    void step() {
      val = slew * val + (1. - slew) * want_val;
    }

    bool changed() {
      return known_want_val != want_val;
    }

    void change_handled() {
      known_want_val = want_val;
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

    Animation() {
    }

    void step() {
      rot_x.step();
      rot_z.step();
      fold_speed.step();
      vision.step();
      alpha.step();
      dragon_points.step();

      angleZ += rot_z;
      angleX += rot_x;
      dragon_fold += fold_speed;
      dragon_r += dragon_r_change;
    }
};

Animation animation;


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


void Dessiner()
{
  static double rotate_shift = 0;
  rotate_shift += .003;

  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity( );

  gluLookAt(3,4,260 - vision*200,0,0,vision * vision * 60,0,0,1);

#if 0
  glPushMatrix();
  //glScaled(6,6,6);

    glPushMatrix();
      glRotated(angleZ,0,.0,1);
      glRotated(angleX,1,1,0);
      dragon1.draw();
    glPopMatrix();

    glPushMatrix();
      glRotated(angleZ,0,0,1);
      glRotated(angleX,0,1,0);
      glTranslated(1, 0, 0);
      dragon2.draw();
    glPopMatrix();

    glPushMatrix();
      glRotated(angleZ,1,0,1);
      glRotated(angleX,0,1,0);
      glTranslated(-1, 0, 0);
      dragon3.draw();
    glPopMatrix();
  glPopMatrix();
#endif

  double sc = 1;

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

      glPushMatrix();
      glTranslated(dd, 0, z1);
      dragon1.draw();
      glPopMatrix();

      glPushMatrix();
      glTranslated(-d, d, z1);
      dragon2.draw();
      glPopMatrix();

      glPushMatrix();
      glTranslated(-d, -d, z1);
      dragon3.draw();
      glPopMatrix();

      glPushMatrix();
      glTranslated(0, 0, z2);
      dragon3.draw();
      glPopMatrix();
    glPopMatrix();
  }

  glFlush();
  SDL_GL_SwapBuffers();
}

int main(int argc, char *argv[])
{
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
  SDL_SetVideoMode(w,h, 32, SDL_OPENGL);

  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective(70,(double)w/h,1,1000);

  glEnable(GL_DEPTH_TEST);
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  dragon1.update(10, .3, 1, 0);
  dragon2.update(10, .3, 1, 0);
  dragon3.update(10, .3, 1, 0);

  Dessiner();

  Uint32 last_time = SDL_GetTicks();
  Uint32 current_time,elapsed_time;
  Uint32 start_time;

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

  while (running)
  {
    start_time = SDL_GetTicks();
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
            }
          }
          break;
      }
    }

    animation.step();

    current_time = SDL_GetTicks();
    elapsed_time = current_time - last_time;
    last_time = current_time;

    angleZ += rot_z;
    angleX += rot_x;
    dragon_fold += fold_speed;
    dragon_r_target += dragon_r_change;
    dragon_r = (.1 * dragon_r_target) + (.9 * dragon_r);

    dragon1.update(dragon_points, dragon_r, dragon_fold, alpha);
    dragon2.update(dragon_points, dragon_r, dragon_fold, alpha);
    dragon3.update(dragon_points, dragon_r, dragon_fold, alpha);
    Dessiner();

    elapsed_time = SDL_GetTicks() - start_time;
    if (elapsed_time < 40)
    {
        SDL_Delay(40 - elapsed_time);
    }
  }

  return 0;
}

