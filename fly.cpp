#include <stdlib.h>
#include <time.h>

#define DRAW_SPHERES 0
#if DRAW_SPHERES
#include <GL/glut.h>
#endif

#define Pf(V) printf(#V "=%f\n", (float)V)
#define Pi(V) printf(#V "=%d\n", (int)V)

#include "draw.h"
#include "audio.h"

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
  Pt n;

  void load(const unsigned int idxs[Face::N], const Pt &normal)
  {
    memcpy(point_idx, idxs, sizeof(point_idx));
    n = normal.unit();
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

  void rotate_e(Pt e)
  {
    nose.rot_e(e);
    top.rot_e(e);
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

  Matrix33 rot_matrix() const
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

#if DRAW_SPHERES
      glutSolidSphere(.5, 20, 16);
#else
      int l;
      Draw::begin(GL_TRIANGLES);

      l = faces.size();
      for (int i = 0; i < l; i++) {
        Face &face = faces[i];
        face.n.glNormal();
        for (int j = 0; j < Face::N; j++) {
          Draw::color_point(points[face.point_idx[j]]);
        }
      }

      Draw::end();
#endif
      glPopMatrix();
    }

    void load_points(const double points[][3], int n_points)
    {
      for (int i = 0; i < n_points; i++)
        add_point().load(points[i]);
    }

    void load_faces(const double points[][3], int n_points,
                    const unsigned int faces[][3], int n_faces)
    {
      load_points(points, n_points);
      for (int i = 0; i < n_faces; i++) {
        Pt a = this->points[faces[i][0]] - this->points[faces[i][1]];
        Pt b = this->points[faces[i][1]] - this->points[faces[i][2]];
        add_face().load(faces[i], b.cross(a).unit());
      }
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
        for (int i = 0; i < rotated_points.size(); i++) {
          rotated_points[i] = points[i];
        }
      }
    }
};

#define ARRAY_SIZE(X) (sizeof(X)/sizeof(*X))

void make_block(Visible &v) {
  static double points[][3] = {
    {-.5, -.5,  .5},
    {-.5, -.5, -.5},
    {-.5,  .5, -.5},
    {-.5,  .5,  .5},
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
    side(1, 4, 7, 2),
    side(4, 5, 6, 7),
    side(5, 0, 3, 6),
    side(3, 2, 7, 6),
    side(5, 4, 1, 0)
  };
#undef side

  v.load_faces(points, ARRAY_SIZE(points),
               faces, ARRAY_SIZE(faces));
}

void color_scheme(Visible &b, const Pt &base_rgb)
{
  for (int i = 0; i < b.points.size(); i ++) {
    b.points[i].c = (base_rgb + Pt::random(-.1, .1)).limit(0, 1);
  }
}

void color_grey(Visible &b, double intens)
{
  Pt g = intens;
  for (int i = 0; i < b.points.size(); i ++) {
    b.points[i].c = g * frandom(0.9, 1.1);
  }
}

class Ufo : public Visible {
  public:

    Mass mass;

    bool moving(double min_len=1e-6)
    {
      return mass.v.len() > min_len;
    }

    bool rotating(double min_len=1e-6)
    {
      return mass.v_ang.len() > min_len;
    }

    bool animate(double min_len=1e-6)
    {
      return moving(min_len) || rotating(min_len);
    }

    virtual void step() {
      pos += mass.v * dt;
      ori.rotate_e(mass.v_ang * dt);
    };

    bool collide(Ufo &o)
    {
      if (!radius_overlaps(o))
        return false;

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

      mass.v *= .59;
      o.mass.v *= .59;
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

      return true;
    }

    void play_bump(double vol=1) const
    {
      double v = max(.01, scale.volume());
      double dens = max(.05, min(10., sqrt(mass.m / v))) * frandom(.9, 1.1);
      double l = scale.len();

      Mix *m = new Mix();
      m->add(new Sine(dens*100./max(.01,scale.x/l), vol/(12), frandom()*2.*M_PI));
      m->add(new Sine(dens*100./max(.01,scale.y/l), vol/(12), frandom()*2.*M_PI));
      m->add(new Sine(dens*100./max(.01,scale.z/l), vol/(12), frandom()*2.*M_PI));

      Audio.play(new Envelope(max(.01,.02/dens), max(.01,.03/dens), max(.01,1./dens), m));
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
    engines_strength = mass.m * 20;

    make_block(*this);
    color_scheme(*this, Pt(.2, .6, .2));

    int l = points.size();
    for (int i = 0; i < l; i++) {
      Point &p = points[i];
      if (p.z < 0) {
        p.x *= .1;
        p.y = max(min(p.y, .01), -.01);
      }
      else {
        p.y = max(min(p.y, .03), -.03);
      }
    }
	}


  virtual void step() {
    mass.v_ang *= .5;
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

class Light {
  public:
    Pt anchor;
    Pt ofs;
    GLuint id;
    bool do_wrap;

    void on() {
      glEnable(id);
    }

    void off() {
      glDisable(id);
    }

    void init(const Pt &at, bool do_wrap=false, const Pt &anchor=Pt())
    {
      this->anchor = anchor;
      this->ofs = at - anchor;
      this->do_wrap = do_wrap;
    }

    void step()
    {
      Pt pos = anchor + ofs;
      glLight(GL_POSITION, pos.x, pos.y, pos.z);
    }

    void glLight(GLuint do_id, float r, float g, float b)
    {
      GLfloat v[] = { r, g, b, 1.0 };
      glLightfv(id, do_id, v);
    }

    void ambient(double v)
    {
      ambient(v, v, v);
    }

    void diffuse(double v)
    {
      diffuse(v, v, v);
    }

    void specular(double v)
    {
      specular(v, v, v);
    }


    void ambient(float r, float g, float b)
    {
      glLight(GL_AMBIENT, r, g, b);
    }

    void diffuse(float r, float g, float b)
    {
      glLight(GL_DIFFUSE, r, g, b);
    }

    void specular(float r, float g, float b)
    {
      glLight(GL_SPECULAR, r, g, b);
    }

    void atten_lin(double f)
    {
      glLightf(id, GL_LINEAR_ATTENUATION, f);
    }

    void atten_quadr(double f)
    {
      glLightf(id, GL_QUADRATIC_ATTENUATION, f);
    }

};


class Lights {
  public:
    static vector<Light> lights;
    static const int GL_LIGHT_ID_LIST[6];

    Lights()
    {
      /* a default light. */
      clear();
      Light &l = add();
      l.init(Pt(5e8, 5e8, -1e8));
      l.ambient(.5);
      l.diffuse(1);
      l.specular(1);
    }


    void step()
    {
      foreach(l, lights) {
        l->step();
      }
    }

    Light &add()
    {
      int lights_n = lights.size();
      if (lights_n >= ARRAY_SIZE(GL_LIGHT_ID_LIST)) {
        printf("Too many lights: %d\n", lights_n + 1);
        exit(1);
      }
      lights.resize(lights_n + 1);
      Light &l = lights.back();
      l.id = GL_LIGHT_ID_LIST[lights_n];
      l.on();
      return l;
    }

    void clear()
    {
      foreach(l, lights) {
        l->off();
      }
      lights.resize(0);
    }

    void wrap(const Pt &center, double dist) {
      foreach(l, lights) {
        if (l->do_wrap)
          l->anchor.wrap_cube(center, dist);
        else
          l->anchor -= center;
      }
    }
};

const int Lights::GL_LIGHT_ID_LIST[6] = {
      GL_LIGHT0,
      GL_LIGHT1,
      GL_LIGHT2,
      GL_LIGHT3,
      GL_LIGHT4,
      GL_LIGHT5
    };
vector<Light> Lights::lights;

class World {
  public:
    Lights lights;
    vector<Ufo*> ufos;
    Pt wrap_ofs;

    World() {
    }

    void add(Ufo &u) {
      ufos.resize(ufos.size() + 1);
      ufos.back() = &u;
    }

    void step() {
      foreach(u, ufos) {
        (*u)->step();
      }
      lights.step();
    }

    void draw() {
      foreach(u, ufos) {
        (*u)->draw();
      }
    }

    void wrap(const Pt &center, double dist) {
      wrap_ofs += center;
      lights.wrap(center, dist);
      foreach(u, ufos) {
        (*u)->pos.wrap_cube(center, dist);
      }
    }

};

class Game {
  public:
    World &world;
    Camera cam;
    bool done;
    int won;
    double redraw_dist;

    Game (World &w) :world(w), done(false), won(0), redraw_dist(0)
    {
    }

    void win()
    {
      won = 1;
      done = true;
    }

    void lose()
    {
      won = -1;
      done = true;
    }

    void quit()
    {
      won = 0;
      done = true;
    }

    virtual void start()
    {
      won = 0;
      done = false;
    }

    virtual void step() = 0;
    virtual void osd_draw() {};

    void run() {
      world.step();
      step();
      draw();
    }

    virtual void on_joy_axis(int axis, double axis_val) {};
    virtual void on_joy_button(int button, bool down) {};

    virtual void on_collision(Ufo &u, Ufo &v) {
      u.play_bump(min(1., 1./((u.pos - cam.from).len()))/3);
      v.play_bump(min(1., 1./((v.pos - cam.from).len()))/3);
    };

    void draw_scene()
    {
      glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );


      glLoadIdentity( );
      
      cam.gluLookAt();

      world.draw();

      if (redraw_dist > .1) {
        Pt dir = cam.at - cam.from;
        Pt p[] = {
          Pt(-redraw_dist, -redraw_dist, -redraw_dist),
          Pt(           0, -redraw_dist, -redraw_dist),
          Pt( redraw_dist, -redraw_dist, -redraw_dist),
          Pt(-redraw_dist,            0, -redraw_dist),
          Pt(           0,            0, -redraw_dist),
          Pt( redraw_dist,            0, -redraw_dist),
          Pt(-redraw_dist,  redraw_dist, -redraw_dist),
          Pt(           0,  redraw_dist, -redraw_dist),
          Pt( redraw_dist,  redraw_dist, -redraw_dist),

          Pt(-redraw_dist, -redraw_dist,            0),
          Pt(           0, -redraw_dist,            0),
          Pt( redraw_dist, -redraw_dist,            0),
          Pt(-redraw_dist,            0,            0),

          Pt( redraw_dist,            0,            0),
          Pt(-redraw_dist,  redraw_dist,            0),
          Pt(           0,  redraw_dist,            0),
          Pt( redraw_dist,  redraw_dist,            0),

          Pt(-redraw_dist, -redraw_dist,  redraw_dist),
          Pt(           0, -redraw_dist,  redraw_dist),
          Pt( redraw_dist, -redraw_dist,  redraw_dist),
          Pt(-redraw_dist,            0,  redraw_dist),
          Pt(           0,            0,  redraw_dist),
          Pt( redraw_dist,            0,  redraw_dist),
          Pt(-redraw_dist,  redraw_dist,  redraw_dist),
          Pt(           0,  redraw_dist,  redraw_dist),
          Pt( redraw_dist,  redraw_dist,  redraw_dist),
        };

        for (int i = 0; i < ARRAY_SIZE(p); i++) {
          if (! p[i].project(dir).zero()) {
            glPushMatrix();
            p[i].glTranslated();
            world.draw();
            glPopMatrix();
          }
        }
      }
    }

    void draw()
    {
      glEnable(GL_LIGHTING);
      draw_scene();
      glDisable(GL_LIGHTING);
      osd_draw();
      glFlush();
      SDL_GL_SwapBuffers();
    }

    void collide() {
      for (int i = 0; i < world.ufos.size(); i++) {
        for (int j = i+1; j < world.ufos.size(); j++) {
          Ufo *a = world.ufos[i];
          Ufo *b = world.ufos[j];
          if (a->moving() || b->moving()) {
            if (a->collide(*b))
              on_collision(*a, *b);
          }
        }
      }
    }
};

class BlockSpace : public Game {
  public:
    double world_r;
    double max_block_size;
    double blocks_count;

    Fly fly;
    vector<FlyingBlock> debris;

    BlockSpace(World &w)
      : Game(w),
        world_r(15),
        max_block_size(2),
        blocks_count(23)
    {
      //redraw_dist = world_r * 2.;//10.0;
    }

    virtual void start()
    {
      Game::start();

      debris.resize(blocks_count);
      foreach(b, debris) {
        b->ori.rotate(Pt::random() * 2 * M_PI);
        b->mass.v_ang = Pt::random() / 2;
        b->pos = Pt::random(-world_r, world_r);
        b->scale = Pt::random(max_block_size / 10, max_block_size);

        double dens = frandom(.5, 1.5);
        b->mass.m = dens * b->scale.volume();

        world.add(*b);
      }

      world.add(fly);

      osd_init();
    }

    virtual void step()
    {
      collide();

      world.wrap(fly.pos + fly.ori.nose * (world_r/2), world_r);

      Pt dir = fly.mass.v.unit();
      if (! dir.zero())
        dir = (dir + fly.ori.nose.unit()) / 2;
      else
        dir = fly.ori.nose;

      cam.look_at(fly.pos + fly.top_lag,
                  dir.unit() * (-3),
                  fly.top_lag);

      osd_update();

      static int skip = 0;
      if ((skip ++) > 20) {
        skip = 0;
        printf("v=%5.2f p=",
               fly.mass.v.len()
              );
        fly.pos.print();
        (world.wrap_ofs + fly.pos).print();
      }
    }

    /* OSD */
    FlyingBlock speed;
    bool draw_want_speed;
    FlyingBlock want_speed;

    void osd_init()
    {
      speed.pos.set(0, -.68, -1);
      want_speed.pos = speed.pos - Pt(0, 0, .002);
    }

    void osd_update() {
      double v = sqrt(fly.mass.v.len()) / 10;
      speed.scale.set(.001 + v/10, .02, .001);
      double c[][4] = {{1, min(.8*v, 1.), .3, .7}, {.95, min(.8*v * 1.1, 1.), .28, .7}};
      speed.load_colors(c, ARRAY_SIZE(c));

      if (fly.do_cruise) {
        v = sqrt(fly.cruise_v) / 10;
        want_speed.scale.set(.001 + v/10, .02, .001);
        double d[][4] = {{1, min(v, 1.), .0, 1}, {.95, min(v * 1.1, 1.), .0, 1}};
        want_speed.load_colors(d, ARRAY_SIZE(d));
        draw_want_speed = true;
      }
      else
        draw_want_speed = false;
    }

    virtual void osd_draw()
    {
      glLoadIdentity();

      speed.draw();
      if (draw_want_speed)
        want_speed.draw();
    }


    /* User input */

    virtual void on_joy_axis(int axis, double axis_val)
    {
      double v;

      switch(axis)
      {
      default:
        printf("%d %f\n", axis, (float)axis_val);
        break;

      case 0:
        fly.roll_y = -(axis_val*axis_val*axis_val) / 10;
        fly.roll_z = (axis_val*axis_val*axis_val) / 50;
        break;

      case 1:
      case 4:
        fly.roll_x = (axis_val*axis_val*axis_val) / 10;
        break;

      case 3:
        fly.roll_z = (axis_val*axis_val*axis_val) / 10;
        break;

      case 5:
        // analog trigger ... -1 == not pressed, 0 = half, 1 = full
        v = (axis_val + 1) / 2;
        fly.propulsion_forward = v;
        break;
      case 2:
        // analog trigger ... -1 == not pressed, 0 = half, 1 = full
        fly.propulsion_break = (axis_val + 1) / 2;
        break;
      }
    }

    void on_joy_button(int button, bool down)
    {
      switch (button) {
      case 3:
        if (! down) {
          if (fly.wings < 1.e-5) {
            fly.wings = .1;
            fly.do_cruise = true;
          }
          else
          if (fly.wings < .9) {
            fly.wings = 1;
            fly.do_cruise = true;
          }
          else {
            fly.wings = 0;
            fly.do_cruise = false;
          }
        }
        break;

      case 1:
        fly.propulsion_forward = down? 1 : 0;
        break;

      case 0:
        fly.propulsion_break = down? 1 : 0;
        break;

      case 2:
        if (down)
          fly.play_bump(1.);
        break;

      default:
        printf("button %d\n", button);
        break;
      }
    }
};


class MoveAllBlocks : public BlockSpace {
  public:
    int inanimate;

    MoveAllBlocks(World &w) :BlockSpace(w)
    {}

    virtual void start()
    {
      BlockSpace::start();

      for (int i = 0; i < debris.size(); i++) {
        FlyingBlock &b = debris[i];
        b.mass.v_ang = 0;

        if (i == 0) {
          b.pos.set(0, 0, -5);
          //b.scale = 1.0;
          //b.mass.m = 1;
        }
        else
        if (i == 1) {
          b.pos.set(3, 3, -5);
          b.mass.v.set(-2, -2, 0);
          //b.scale = 1.001;
          //b.mass.m = 1.5;
        }

        color_grey(b, .1);
      }

      move_osd_init();
    }

    void count_inanimate()
    {
      inanimate = 0;
      foreach(_u, world.ufos) {
        Ufo &u = **_u;
        if (! u.animate(1e-2)) {
          if (u.animate(0)) {
            u.mass.v = 0;
            u.mass.v_ang = 0;
            color_grey(u, .1);
          }
          inanimate ++;
        }
      }
    }

    virtual void step()
    {
      count_inanimate();
      if (inanimate == 0)
        win();

      BlockSpace::step();

      move_osd_update();

      foreach(u, world.ufos) {
        (*u)->mass.v *= .998;
        (*u)->mass.v_ang *= .998;
      }
    }

    virtual void on_collision(Ufo &u, Ufo &v) {
      BlockSpace::on_collision(u, v);
      u.mass.v_ang = Pt::random();
      v.mass.v_ang = Pt::random();
      color_scheme(u, Pt::random(0.1, 1));
      color_scheme(v, Pt::random(0.1, 1));
    };


    // OSD

    int was_inanimate;
    double show_got_one;
    FlyingBlock got_one;
    FlyingBlock bar_all;
    FlyingBlock bar_inanimate;

    void move_osd_init()
    {
      bar_all.pos.set(0, -.71, -1.004);
      bar_all.scale.set(.0001 + 2., .03, .001);

      bar_inanimate.pos.set(0, -.71, -1.002);

      got_one.pos.set(0, -.71, -1);
      show_got_one = 0;

      double c[][4] = {{.1, 0.7, .1, .4},};
      bar_all.load_colors(c, ARRAY_SIZE(c));
      double d[][4] = {{.7, 0.1, .1, .4},};
      bar_inanimate.load_colors(d, ARRAY_SIZE(c));

      count_inanimate();
      was_inanimate = inanimate;
    }

    void move_osd_update()
    {
      if (show_got_one > 0) {
        show_got_one -= dt;
        double c[][4] = {{.5, 1, .5, min(show_got_one, 1.)}, {.4, 1, .4, min(show_got_one, 1.)}};
        got_one.load_colors(c, ARRAY_SIZE(c));
      };
      if (was_inanimate != inanimate) {
        printf("got %d of %d\n", inanimate, world.ufos.size());
        show_got_one = 1.1;
        was_inanimate = inanimate;
      }
      bar_inanimate.scale.set(.0001 + 2. * inanimate / world.ufos.size(), .03, .001);
      got_one.scale.set(.0001 + 2. * inanimate / world.ufos.size(), .03, .001);
    }

    virtual void osd_draw()
    {
      BlockSpace::osd_draw();

      bar_all.draw();
      bar_inanimate.draw();
      if (show_got_one > 0) {
        got_one.draw();
      }
    }
};



class FindTheLight : public BlockSpace {
  public:
    Pt light_pos;

    FindTheLight(World &w) :BlockSpace(w)
    {
    }

    virtual void start()
    {
      light_pos = Pt::randoml(500, 600);

      world.lights.clear();
      {
        Light &l = world.lights.add();
        l.init(light_pos);
        l.ambient(.8);
        l.diffuse(1);
        l.specular(1);
        l.atten_quadr(.001);
      }

      {
        Light &l = world.lights.add();
        l.init(light_pos);
        l.ambient(.0);
        l.diffuse(.1);
        l.specular(0);
      }

      BlockSpace::start();
    }

};


char *audio_path = NULL;

int W = 1920;
int H = 900;//1200;

SDL_Surface *screen = NULL;


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

  ip.random_seed = time(NULL);

  while (1) {
    c = getopt(argc, argv, "hf:g:r:");
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
, W, H, want_fps
);
    if (error)
      return 1;
    return 0;
  }

  const int maxpixels = 1e4;

#if DRAW_SPHERES
  int _argc = 1;
   glutInit(&_argc, argv);
#endif

  if ((W < 3) || (W > maxpixels) || (H < 3) || (H > maxpixels)) {
    fprintf(stderr, "width and/or height out of bounds: %dx%d\n", W, H);
    exit(1);
  }

  SDL_Event event;

  SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_JOYSTICK);

  const int n_joysticks = SDL_NumJoysticks();
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
  SDL_WM_SetCaption("fly", NULL);
  screen = SDL_SetVideoMode(W,H, 32, SDL_OPENGL | SDL_FULLSCREEN);
  SDL_ShowCursor(SDL_DISABLE);

  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective(80, (double)W/H, .5, 1000);

  glClearColor (0.0, 0.0, 0.0, 0.0);
  glShadeModel (GL_SMOOTH);

  glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE ) ;
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_NORMALIZE);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glMatrixMode( GL_MODELVIEW );


  printf("seed %d\n", (int)ip.random_seed);
  srandom(ip.random_seed);

  Audio.start();

  World world;

  MoveAllBlocks game(world);

  game.start();
  game.draw();

  Uint32 last_time = SDL_GetTicks();
  Uint32 current_time,elapsed_time;
  Uint32 start_time;

  float want_frame_period = (want_fps > .1? 1000. / want_fps : 0);
  float last_ticks = (float)SDL_GetTicks() - want_frame_period;

  while (running)
  {
    game.run();
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

    while (running) {
      SDL_Event event;
      while (running && SDL_PollEvent(&event)) 
      {
        switch(event.type)
        {

        case SDL_JOYAXISMOTION:
          game.on_joy_axis(event.jaxis.axis,
                           ((double)event.jaxis.value) / 32768.);
          break;


        case SDL_JOYBUTTONDOWN:
        case SDL_JOYBUTTONUP:
          game.on_joy_button(event.jbutton.button,
                             event.type == SDL_JOYBUTTONDOWN);
          break;

        case SDL_JOYBALLMOTION:  /* Handle Joyball Motion */
          printf("%2d: ball %d += %d, %d\n",
                 event.jball.which, event.jball.ball,
                 event.jball.xrel, event.jball.yrel);
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

  Audio.stop();

  running = false;

  printf("\n");
  printf("%d frames rendered\n", frames_rendered);

  return 0;
}

