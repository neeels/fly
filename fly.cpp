#include <stdlib.h>
#include <time.h>

#define Pf(V) printf(#V "=%f\n", (float)V)
#define Pi(V) printf(#V "=%d\n", (int)V)

#include "ctrl_layers.h"
#include "draw.h"

static AsQuads as_quads;
static AsLines as_lines;


bool running = true;
volatile int frames_rendered = 0;
volatile int avg_frame_period = 0;
volatile double dt;
#define AVG_SHIFTING 3
float want_fps = 25;

struct Face {
  static const unsigned int N = 3;
  unsigned int point_idx[Face::N];

  void load(const unsigned int idxs[Face::N])
  {
    memcpy(point_idx, idxs, sizeof(point_idx));
  }
};


struct Orientation {
  Pt nose;
  Pt top;

  Orientation()
    : nose(0, 0, -1),
      top(0, 1, 0)
  {}

  Pt right() {
    return nose.cross(top);
  }

  void rotate(double roll_x, double roll_y, double roll_z)
  {
      nose.rot_about(top, roll_y);
      Pt r = right();
      nose.rot_about(r, roll_x);
      top.rot_about(r, roll_x);
      top.rot_about(nose, roll_z);
  }

  void rotate(Pt r3)
  {
    rotate(r3.x, r3.y, r3.z);
  }

  void glRotated() const {
    rot3().glRotated();
  }

  void check_normals() {
      nose = nose.unit();
      top = top.unit();

      double is_ortho = nose.dot(top);
      if (fabs(is_ortho) > 1.e-6)
        printf("ORTHO! %g\n", is_ortho);
  }


  Pt rot3() const {
    double l = nose.len();
    if (fabs(l) < 1.e-6) {
      return Pt();
    }

    /* "first roll left/right about the nose axis by angle az (around z axis),
        then turn the nose up/down by angle ax (around the x axis),
        then turn left/right by angle ay (around the vertical y axis)."

        ^ y
        |
        | az
     ax +------> x
       / ay
      v
      z

      This vector is the result of those rotations, the top vector points
      "upwards", and we're trying to find out ax, ay, and az.
     */

    Pt u = nose.unit();

    /* ax is the angle between the y=0 plane and this vector.
                                              
                  .                             
                 /|                                
               /  |                               
      y ^    /    |
        |  /      |
        |/        |
        +__ ax    |. . . . .  --> x
       /   --__   |
      v        --_|
      z

      */
    double ax = asin(u.y);


    /* ay is the angle between the z axis and the projection of this vector
       onto the y=0 plane. 
                  .                             
                 /                                 
               /                                  
      y ^    /     
        |  /       
        |/         
        +__. . . . . . . . .  --> x
       /   --__    
      /  ay    --_ 
     /-------------
    v
    z

    */
    double ay;
    bool pos = u.z > 1e-6;
    bool neg = u.z < -1e-6;
    if (!(pos || neg))
      ay = (u.x > 0? 1. : -1.) * M_PI/2;
    else
      ay = atan(u.x / fabs(u.z));

    if (neg)
      ay = M_PI - ay;

    /* Construct the top vector that would have az == 0, according to ax and
       ay: it is this vector but with an additional ax turn of pi/2 (90°).

       So it is a unit vector, starting from directly on y axis (pointing upward),
       tilted "backward" (towards negative z axis) by ax,
       and then turned around the y axis by ay.

       Turns out the zero-top's new y coordinate is the same as the length of
       this vector's projection onto the y=0 plane: ly = sqrt(x²+z²)

                    | y
             .___   |         . 
              \  ---+ly      /|                                
               \    ^      /  |                               
                \   |    /    |
                 \ax|  /      |
                  \ |/        |
          . . . . . +__ax     |. . . . .  --> x
                   /   --__   |
                  /        --_|
                 /            > ly
                v             
                z

       Also, the new distance from the z axis == this vector's y coordinate:

                    | 
          "y".___   |         .y
              \  -->+    /   /|                                
               \    |   /  /  |                               
                \   |  / /    |
             -_<-\ay+>//      |
                -.\ |/        |
          . . . . . +__       |. . . . .  --> x
                   /   --__   |
                  /<-ay--> --_V
                 /            -
                v               
                z

      So use sin and cos to get x and z coords.
     */

    Pt top_zero(-u.y * sin(ay),
                sqrt(u.x*u.x + u.z*u.z),
                -u.y * cos(ay));

    /* Find the angle between the "zero" top and the user supplied top vector.
     * Both top and top_zero must be in the plane perpendicular to this vector.
     * Remove any component from top that's pointing in this vector's dir. */
    Pt topu = top.unit();
    if (fabs(u.dot(topu)) > 1e-6) {
      printf("ORTHO! %g " __FILE__ " line %d\n", fabs(u.dot(topu)), __LINE__);
    }
    //topu -= u * topu.project(u);

    /* Angle of the "zero" top to the plane-ized and unit-ized user supplied
     * top, angle sign relative to this vector. */
    double az = top_zero.angle(topu, u);

    /* I measured the angle from the z axis rightwards. That's counter the
     * canonical angle direction OpenGL uses, so let's correct for that... */
    return Pt(M_PI - ax, ay, - az);
  }

  Matrix33 rot_matrix()
  {
    return Matrix33::from_rot3(rot3());
  }

  bool operator!=(const Orientation &other) const {
    return (nose != other.nose) && (top != other.top);
  }

  bool operator==(const Orientation &other) const {
    return !(operator!=(other));
  }
};

class Visible {
  private:
    double _radius;
    vector<Point> rotated_points;
    Orientation rotated_ori;
    bool shape_changed;
  public:
    vector<Point> points;
    vector<Face> faces;
    Pt scale;
    Pt pos;
    Orientation ori;

    Visible()
      :scale(1, 1, 1),
       shape_changed(false),
       _radius(0)
    {
      scale.set(1, 1, 1);
    }

    Point &add_point()
    {
      points.resize(points.size() + 1);
      shape_changed = true;
      return points.back();
    }

    Face &add_face()
    {
      faces.resize(faces.size() + 1);
      return faces.back();
    }

    void draw()
    {
      glPushMatrix();

      pos.glTranslated();
      ori.glRotated();
      scale.glScaled();

      int l;
      Draw::begin(GL_TRIANGLES);

      l = faces.size();
      for (int i = 0; i < l; i++) {
        Face &face = faces[i];
        for (int j = 0; j < Face::N; j++) {
          Draw::color_point(points[face.point_idx[j]]);
        }
      }

      Draw::end();

      glPopMatrix();
    }

    void load_points(const double points[][3], int n_points)
    {
      for (int i = 0; i < n_points; i++)
        add_point().load(points[i]);
    }

    void load_faces(const unsigned int faces[][3], int n_faces)
    {
      for (int i = 0; i < n_faces; i++)
        add_face().load(faces[i]);
    }

    void load_colors(const double colors[][4], int n_colors, bool cycle=true)
    {
      int p = 0;
      while (p < points.size()) {
        for (int i = 0; i < n_colors; i++, p++) {
          points[p].c.load(colors[i]);
        }
        if (! cycle)
          break;
      }
    }

    double radius()
    {
      if (shape_changed) {
        shape_changed = false;
        _radius = 0;
        foreach(p, points) {
          _radius = max(_radius, p->scaled(scale).len());
        }
      }
      return _radius;
    }

    bool radius_overlaps(Visible &v)
    {
      return (pos - v.pos).len() < (radius() + v.radius());
    }

  private:
    void update_rotated_points()
    {
      if (rotated_ori != ori) {
        rotated_ori = ori;
        Pt rot3 = rotated_ori.rot3();

        rotated_points.resize(points.size());
        foreach(p, rotated_points) {
          *p = points[(long unsigned int)_i];

          //rotated_points[_i] = p->rot_about
        }

      }
    }
};

#define ARRAY_SIZE(X) (sizeof(X)/sizeof(*X))

void make_block(Visible &v) {
  static double points[][3] = {
    {-.5, -.5, -.5},
    {-.5, -.5,  .5},
    {-.5,  .5,  .5},
    {-.5,  .5, -.5},
    { .5, -.5, -.5},
    { .5, -.5,  .5},
    { .5,  .5,  .5},
    { .5,  .5, -.5},
  };

  // side: a rectangle from two triangles
#define side(A, B, C, D) \
  {A, B, C}, \
  {C, D, A} 

  static unsigned int faces[][3] = {
    side(0, 1, 2, 3),
    side(0, 1, 5, 4),
    side(1, 2, 6, 5),
    side(4, 5, 6, 7),
    side(2, 3, 7, 6),
    side(3, 0, 4, 7)
  };
#undef side

  v.load_points(points, ARRAY_SIZE(points));
  v.load_faces(faces, ARRAY_SIZE(faces));
}

void color_scheme(Visible &b, const Pt &base_rgb) {
  for (int i = 0; i < b.points.size(); i ++) {
    b.points[i].c = (base_rgb + Pt::random(-.1, .1)).limit(0, 1);
  }
}

class Ufo : public Visible {
  public:

    Mass mass;

    /* accumulated for each step */
    Pt a;

    bool moving()
    {
      return mass.v.len() > 1.e-6;
    }

    virtual void step() {
      pos += mass.v * dt;
    };

    void apply_force(Pt at, Pt f)
    {
      a += f / mass.m;
    }

    void step_accelerate(double dt)
    {
      mass.accelerate(a, dt);
      a = 0;
    }

    void collide(Ufo &o)
    {
      if (!radius_overlaps(o))
        return;

      Pt dir = o.pos - pos;
      Pt n = dir.unit();

      /*
      printf("\n");
      Pt i1 = mass.impulse() + o.mass.impulse();
      printf("imp\t");mass.impulse().print();
      printf("o.imp\t");o.mass.impulse().print();
      printf("i1\t"); i1.print();
      printf("m %f  o.m %f\n", mass.m, o.mass.m);
      printf("v\t");mass.v.print();
      printf("o.v\t"); o.mass.v.print();
      printf("n\t"); n.print();
      printf("vdiff\t"); (o.mass.v - mass.v).print();
      printf("vdiffp\t"); (o.mass.v - mass.v).project(n).print();
      printf("vdiffl\t%f\n", (o.mass.v - mass.v).project(n).len());
      */
      Pt J = n * (2.* (mass.v - o.mass.v).project(n).len() / (1./mass.m + 1./o.mass.m));

      mass.v -= J / mass.m;
      o.mass.v += J / o.mass.m;
      /*
      printf("J\t"); J.print();
      Pt i2 = mass.impulse() + o.mass.impulse();
      printf("v\t");mass.v.print();
      printf("o.v\t"); o.mass.v.print();
      printf("imp\t");mass.impulse().print();
      printf("o.imp\t");o.mass.impulse().print();
      printf("i2\t"); i2.print();
      printf("i1 %6.2f  i2 %6.2f\n",
             i1.len(), i2.len());
             */

      Pt at = ((pos + n * radius()) + (o.pos - n * o.radius())) / 2;
      /*
      printf("pos\t");pos.print();
      printf("at\t");at.print();
      printf("o.pos\t");o.pos.print();
*/
      pos = at - n * (1.01 * radius());
      o.pos = at + n * (1.01 * o.radius());
    }
};

class FlyingBlock : public Ufo {
  public:
    FlyingBlock() {
      make_block(*this);
      color_scheme(*this, Pt::random(0.1, 1));
    }
};

class Fly : public Ufo {
  public:
  Pt top_lag;

  Param top_angle;
  Param roll_x;
  Param roll_y;
  Param roll_z;

  Param propulsion_forward;
  Param propulsion_break;

  Param wings;

  bool do_cruise;
  Param cruise_v;

  double engines_strength;

	Fly() : Ufo() {
    ori.nose.set(0, 0, -1);
    ori.top.set(0, 1, 0);
    top_lag = ori.top;
    roll_x.slew = .94;
    roll_y.slew = .94;
    roll_z.slew = .94;
    propulsion_forward.slew = .1;
    propulsion_break.slew = .1;
    cruise_v.slew = 0;
    wings.limit(0, 1);
    wings = .1;

    do_cruise = true;
    cruise_v = .05;
    mass.v.set(0, 0, -.05);
    mass.m = 5000;
    engines_strength = mass.m * 100;

    make_block(*this);
    color_scheme(*this, Pt(.2, .6, .2));

    int l = points.size();
    for (int i = 0; i < l; i++) {
      Point &p = points[i];
      if (p.z < 0) {
        p.x *= .1;
      }
      p.y = max(min(p.y, .1), -.1);
    }
	}


  virtual void step() {
      roll_x.step();
      roll_y.step();
      roll_z.step();
      top_angle.step();
      propulsion_forward.step();
      propulsion_break.step();
      wings.step();
      cruise_v.step();

      ori.rotate(roll_x, roll_y, roll_z);
      ori.check_normals();

      const int div = 10;
      top_lag = (top_lag * (div-1) + ori.top) / div;


      Pt wings_force;

      if (dt < 1.e-5)
        return;

      if (wings > 1.e-5) {
        Pt v_want = ori.nose * mass.v.len();
        wings_force = (v_want - mass.v) * mass.m * wings / dt;
      }

      double want_propulsion_forward = propulsion_forward;
      if (do_cruise) {
        if ((want_propulsion_forward > 1e-6) || (propulsion_break > 1e-6)) {
          /* user is changing speed. Take current speed as new desired speed. */
          cruise_v = mass.v.len();
        }
        else {
          double diff = cruise_v - mass.v.len();
          if (diff > 1e-6) {
            /* we're slowing down below cruising speed. Hit the accelerator a bit. */
            want_propulsion_forward = (6. * mass.m * diff) / engines_strength; // get percentile
            if (want_propulsion_forward < propulsion_forward)
              want_propulsion_forward = propulsion_forward;
          }
        }
      }

      // 1. is pedal to the metal
      want_propulsion_forward = min(1., want_propulsion_forward);

      Pt forward_want = ori.nose * (want_propulsion_forward * engines_strength);
      Pt break_want = mass.v * mass.m * (-propulsion_break);

      Pt propulsion = forward_want + break_want;
      if (propulsion.len() > engines_strength)
        propulsion = propulsion.unit() * engines_strength;

      mass.accelerate(wings_force + propulsion, dt);

      Ufo::step();

      static int skip = 0;
      if ((skip ++) > 50) {
        skip = 0;
        printf("v=%5.2f p=%5.2f g=%5.2f\n",
               mass.v.len(),
               propulsion.len(),
               wings_force.len()
              );
      }
  }

};


struct Camera {
  Pt at;
  Pt from;
  Pt top;

  void look_at(const Pt &look_at, const Pt &from_rel, const Pt &top) {
    at = look_at;
    from = at + from_rel;
    this->top = top;
  }

  void gluLookAt() {
    ::gluLookAt(from.x, from.y, from.z,
                at.x, at.y, at.z,
                top.x, top.y, top.z);
  }
};

const int world_r = 1000;
const int max_block_size = 30;

class World {
  public:
    Fly fly;
    Camera cam;
    FlyingBlock debris[10 * world_r / max_block_size];
    int inanimate;

    vector<Ufo*> ufos;

    World() {
      for (int i = 0; i < ARRAY_SIZE(debris); i++) {
        FlyingBlock &b = debris[i];
        b.ori.rotate(Pt::random() * 2 * M_PI);

        if (i == 0) {
          b.pos.set(0, 0, -5);
          b.scale = 1.0;
        }
        else
        if (i == 1) {
          b.pos.set(3, 3, -5);
          b.mass.v.set(-.7, -.7, 0);
          b.scale = 1.001;
        }
        else {
          b.pos = Pt::random(-world_r, world_r);
          b.scale = Pt::random(1, max_block_size);
        }

        b.mass.m = b.scale.x * b.scale.y * b.scale.z;

        add(b);
      }



      add(fly);
    }

    void add(Ufo &u) {
      ufos.resize(ufos.size() + 1);
      ufos.back() = &u;
    }

    void load() {
    }

    void step() {
      inanimate = 0;
      foreach(u, ufos) {
        (*u)->step();
        if (!(*u)->moving())
          inanimate ++;
      }

      static int skip = 0;
      if (skip++ > 50) {
        skip = 0;
        if (! inanimate) {
          printf("YOU WIN!");
        }
        printf("%d / %d\n", ufos.size() - inanimate, ufos.size());
      }

      wrap(fly.pos + fly.ori.nose * (world_r/2), world_r);

      collide();

      /*
      foreach(u, ufos) {
        (*u)->step_accelerate(dt);
      }
      */

      Pt dir = fly.mass.v.unit();
      if (! dir.zero())
        dir = (dir + fly.ori.nose.unit()) / 2;
      else
        dir = fly.ori.nose;

      cam.look_at(fly.pos + fly.top_lag,
                  dir.unit() * (-3),
                  fly.top_lag);
    }


    void draw() {
      for (int i = 0; i < ufos.size(); i++) {
        ufos[i]->draw();
      }
    }

    void wrap(const Pt &center, double dist) {
      for (int i = 0; i < ufos.size(); i++) {
        ufos[i]->pos.wrap(center, dist);
      }
    }

    void collide() {
      for (int i = 0; i < ufos.size(); i++) {
        for (int j = i+1; j < ufos.size(); j++) {
          if (ufos[i]->moving() || ufos[j]->moving())
            ufos[i]->collide(*ufos[j]);
        }
      }
    }
};


class Osd {
  public:
    DrawBank bank;

    FlyingBlock speed;
    bool draw_want_speed;
    FlyingBlock want_speed;

    int was_inanimate;
    double show_got_one;
    FlyingBlock got_one;
    FlyingBlock all;
    FlyingBlock inanimates;

    Osd()
    {
      speed.pos.set(-1, -.7, -1);
      want_speed.pos.set(-1, -.7, -1.002);

      all.pos.set(0, -.71, -1.004);
      all.scale.set(.0001 + 2., .03, .001);

      inanimates.pos.set(0, -.71, -1.002);
      
      got_one.pos.set(0, -.71, -1);
      show_got_one = 0;

      double c[][4] = {{.1, 0.7, .1, .4},};
      all.load_colors(c, ARRAY_SIZE(c));
      double d[][4] = {{.7, 0.1, .1, .4},};
      inanimates.load_colors(d, ARRAY_SIZE(c));
    }

    void update(const World &world) {
      double v = sqrt(world.fly.mass.v.len()) / 10;
      speed.scale.set(.001 + v/10, .02, .001);
      double c[][4] = {{1, min(.8*v, 1.), .3, .7}, {.95, min(.8*v * 1.1, 1.), .28, .7}};
      speed.load_colors(c, ARRAY_SIZE(c));

      if (world.fly.do_cruise) {
        v = sqrt(world.fly.cruise_v) / 10;
        want_speed.scale.set(.001 + v/10, .02, .001);
        double d[][4] = {{1, min(v, 1.), .0, 1}, {.95, min(v * 1.1, 1.), .0, 1}};
        want_speed.load_colors(d, ARRAY_SIZE(d));
        draw_want_speed = true;
      }
      else
        draw_want_speed = false;

      if (show_got_one > 0) {
        show_got_one -= dt;
        double e[][4] = {{.5, 1, .5, min(show_got_one, 1.)}, {.4, 1, .4, min(show_got_one, 1.)}};
        got_one.load_colors(e, ARRAY_SIZE(c));
      };
      if (was_inanimate != world.inanimate) {
        printf("CHANGE %d %d\n", was_inanimate, world.inanimate);
        show_got_one = 1.1;
        was_inanimate = world.inanimate;
      }
      inanimates.scale.set(.0001 + 2. * world.inanimate / world.ufos.size(), .03, .001);
      got_one.scale.set(.0001 + 2. * world.inanimate / world.ufos.size(), .03, .001);
    }

    void draw() {
      speed.draw();
      if (draw_want_speed)
        want_speed.draw();
      all.draw();
      inanimates.draw();
      if (show_got_one > 0) {
        got_one.draw();
      }
    }
};

void draw_scene(World &world, Osd &osd)
{
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity( );

  world.cam.gluLookAt();
  world.draw();

  glLoadIdentity( );
  osd.draw();
  glFlush();
  SDL_GL_SwapBuffers();
}


FILE *out_stream = NULL;
FILE *out_params = NULL;
FILE *in_params = NULL;

char *audio_path = NULL;

int W = 1920;
int H = 900;//1200;

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

World *gl_world;

void on_joy_axis(ControllerState &ctrl, int axis, double axis_val) {
  double v;
  switch(axis)
  {
  default:
    printf("%d %f\n", axis, (float)axis_val);
    break;
  case 3:
    gl_world->fly.roll_y = -(axis_val*axis_val*axis_val) / 10;
    break;
  case 1:
  case 4:
    gl_world->fly.roll_x = (axis_val*axis_val*axis_val) / 10;
    break;
  case 0:
    gl_world->fly.roll_z = (axis_val*axis_val*axis_val) / 10;
    break;
  case 5:
    // analog trigger ... -1 == not pressed, 0 = half, 1 = full
    // accel
#if 0
    world.fly.velocity.change = ((axis_val + 1) / 2) / 200;
#else
    v = (axis_val + 1) / 2;
    gl_world->fly.propulsion_forward = .2 * v;
#endif
    break;
  case 2:
    // analog trigger ... -1 == not pressed, 0 = half, 1 = full
    // break
#if 0
    world.fly.velocity.change = -((axis_val + 1) / 2) / 200;
#else
    gl_world->fly.propulsion_break = (axis_val + 1) / 2;
#endif
    break;
  }
}

void on_joy_button(ControllerState &ctrl, int button, bool down) {

  switch (button) {
  case 3:
    if (! down) {
      if (gl_world->fly.wings < 1.e-5) {
        gl_world->fly.wings = .1;
        gl_world->fly.do_cruise = true;
      }
      else
      if (gl_world->fly.wings < .9) {
        gl_world->fly.wings = 1;
        gl_world->fly.do_cruise = true;
      }
      else {
        gl_world->fly.wings = 0;
        gl_world->fly.do_cruise = false;
      }
    }
    break;

  case 0:
    gl_world->fly.propulsion_forward = down? 1 : 0;
    break;

  case 1:
    gl_world->fly.propulsion_break = down? 1 : 0;
    break;

  default:
    printf("button %d\n", button);
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

  int c;

  char *out_stream_path = NULL;
  char *out_params_path = NULL;
  char *in_params_path = NULL;

  ip.random_seed = time(NULL);

  while (1) {
    c = getopt(argc, argv, "hf:g:r:p:i:o:O:");
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
  screen = SDL_SetVideoMode(W,H, 32, SDL_OPENGL | SDL_FULLSCREEN);
  SDL_ShowCursor(SDL_DISABLE);

  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective(80,(double)W/H,.5,world_r);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


  printf("seed %d\n", (int)ip.random_seed);
  srandom(ip.random_seed);

  World world;
  gl_world = &world;

  Osd osd;
  osd.was_inanimate = world.ufos.size() - 2;

  printf("%d\n", (int)(frandom() * 1000));

  world.load();

  draw_scene(world, osd);

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

    world.step();
    osd.update(world);

    draw_scene(world, osd);
    frames_rendered ++;

    {
      static int last_ticks2 = 0;

      int t = SDL_GetTicks();
      int elapsed = t - last_ticks2;

      last_ticks2 = t;

      avg_frame_period -= avg_frame_period >>AVG_SHIFTING;
      avg_frame_period += elapsed;

      dt = 1.e-3 * elapsed;
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

                  case 'f':
                  case 13:
                    SDL_WM_ToggleFullScreen(screen);
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

