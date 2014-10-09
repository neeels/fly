
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
#include "palettes.h"

#define Pf(V) printf(#V "=%f\n", (float)V)
#define Pi(V) printf(#V "=%d\n", (int)V)


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

    void random() {
      x = -1. + 2.*frandom();
      y = -1. + 2.*frandom();
      z = -1. + 2.*frandom();
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

    void scale3(Pt &p) {
      x *= p.x;
      y *= p.y;
      z *= p.z;
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

    Pt& operator+=(const double f) {
      x += f;
      y += f;
      z += f;
      return *this;
    }

    Pt& operator-=(const double f) {
      x -= f;
      y -= f;
      z -= f;
      return *this;
    }

    Pt& operator*=(const double f) {
      x *= f;
      y *= f;
      z *= f;
      return *this;
    }

    Pt& operator/=(const double f) {
      x /= f;
      y /= f;
      z /= f;
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

    void blend(Color &c, double fade) {
      double rfade = 1.0 - fade;
      r = rfade * r  +  fade * c.r;
      g = rfade * g  +  fade * c.g;
      b = rfade * b  +  fade * c.b;
      a = rfade * a  +  fade * c.a;
    }

    void random(double cmin=.3, double cmax=1.) {
      r = frandom();
      g = frandom();
      b = frandom();
      a = 1;
      double is_intens = max(.01, (r + g + b) / 3.);
      double intens_change = (cmin + (cmax-cmin)*frandom()) / is_intens;

      r = intens_change * r;
      g = intens_change * g;
      b = intens_change * b;
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

    void print() {
      printf("r%f g%f b%f a%f\n", r,g,b,a);
    }

};



class Point : public Pt{
  public:

    Color c;

    Point() : Pt() {
      c.random();
    }

    Point(double x, double y, double z) : Pt(x, y, z) {
      c.random();
    }

    void draw(double alpha=1., double greying=0) {
      c.glColor(alpha, greying);
      glVertex3d();
    }

    void random() {
      Pt::random();
      c.random();
    }
};


class FilterContext {
  public:
    int point_i;
    int color_i;
    int face_i;
    int particle_i;
    int cloud_i;
    int frame_i;

    GLenum what;

    FilterContext() {
      clear();
      frame_i = 0;
    }

    void clear() {
      point_i = 0;
      color_i = 0;
      face_i = 0;
      particle_i = 0;
      cloud_i = 0;
      what = -1;
    }
};

class Filter {
  public:
    FilterContext *fc;

    Filter() {
      fc = NULL;
    }

    virtual void ign() {};
};

typedef enum {
  level_unknown = 0,
  level_camera,
  level_cloud,
  level_particle
} level_e;


class PointFilter : public Filter {
  public:
    virtual void point(Pt &p) = 0;
};

class ColorFilter : public Filter {
  public:
    virtual void color(Color &c) = 0;
};

class TranslateFilter : public Filter {
  public:
    virtual void translate(Pt &p, level_e l) = 0;
};

class RotateFilter : public Filter {
  public:
    virtual void rotate(Pt &p, level_e l) = 0;
};

class ScaleFilter : public Filter {
  public:
    virtual void scale(Pt &p, level_e l) = 0;
};

template <typename V, typename T>
void vector_rm(V &v, T &item) {
  typename V::iterator i = v.begin();
  for (; i != v.end(); i++) {
    if (*i == item) {
      v.erase(i);
      return;
    }
  }
}

class DrawBank {
  public:
    FilterContext fc;

    vector<PointFilter*> point_filters;
    vector<ColorFilter*> color_filters;
    vector<TranslateFilter*> translate_filters;
    vector<RotateFilter*> rotate_filters;
    vector<ScaleFilter*> scale_filters;

    void start() {
      fc.clear();
    }

    void add(Filter &f) {
      add(&f);
    }

    void add(Filter *f) {
      if (dynamic_cast<PointFilter*>(f)) {
        f->fc = &fc;
        point_filters.push_back((PointFilter*)f);
      }
      else if (dynamic_cast<ColorFilter*>(f)) {
        f->fc = &fc;
        color_filters.push_back((ColorFilter*)f);
      }
      else if (dynamic_cast<TranslateFilter*>(f)) {
        f->fc = &fc;
        translate_filters.push_back((TranslateFilter*)f);
      }
      else if (dynamic_cast<RotateFilter*>(f)) {
        f->fc = &fc;
        rotate_filters.push_back((RotateFilter*)f);
      }
      else if (dynamic_cast<ScaleFilter*>(f)) {
        f->fc = &fc;
        scale_filters.push_back((ScaleFilter*)f);
      }
    }

    void remove(Filter &f) {
      remove(&f);
    }

    void remove(Filter *f) {
      if (dynamic_cast<PointFilter*>(f)) {
        f->fc = &fc;
        PointFilter *ff = (PointFilter*)f;
        vector_rm<vector<PointFilter*>,PointFilter*>(point_filters, ff);
      }
      else if (dynamic_cast<ColorFilter*>(f)) {
        f->fc = &fc;
        ColorFilter *ff = (ColorFilter*)f;
        vector_rm<vector<ColorFilter*>,ColorFilter*>(color_filters, ff);
      }
      else if (dynamic_cast<TranslateFilter*>(f)) {
        f->fc = &fc;
        TranslateFilter *ff = (TranslateFilter*)f;
        vector_rm<vector<TranslateFilter*>,TranslateFilter*>(translate_filters, ff);
      }
      else if (dynamic_cast<RotateFilter*>(f)) {
        f->fc = &fc;
        RotateFilter *ff = (RotateFilter*)f;
        vector_rm<vector<RotateFilter*>,RotateFilter*>(rotate_filters, ff);
      }
      else if (dynamic_cast<ScaleFilter*>(f)) {
        f->fc = &fc;
        ScaleFilter *ff = (ScaleFilter*)f;
        vector_rm<vector<ScaleFilter*>,ScaleFilter*>(scale_filters, ff);
      }
    }

    void point(Pt &p) {
      if (point_filters.size()) {
        Pt q = p;
        for (int i = 0; i < point_filters.size(); i++)
          point_filters[i]->point(q);
        q.glVertex3d();
      }
      else
        p.glVertex3d();
      fc.point_i ++;
    }

    void color(Color &c) {
      if (color_filters.size()) {
        Color d = c;
        for (int i = 0; i < color_filters.size(); i++)
          color_filters[i]->color(d);
        d.glColor();
      }
      else
        c.glColor();
      fc.color_i ++;
    }

    void color_point(Point &p) {
      color(p.c);
      point(p);
    }

    void translate(Pt &p, level_e l) {
      if (translate_filters.size()) {
        Pt q = p;
        for (int i = 0; i < translate_filters.size(); i++)
          translate_filters[i]->translate(q, l);
        q.glTranslated();
      }
      else
        p.glTranslated();
    }

    void rotate(Pt &p, level_e l) {
      if (rotate_filters.size()) {
        Pt q = p;
        for (int i = 0; i < rotate_filters.size(); i++)
          rotate_filters[i]->rotate(q, l);
        q.glRotated();
      }
      else
        p.glRotated();
    }

    void scale(Pt &p, level_e l) {
      if (scale_filters.size()) {
        Pt q = p;
        for (int i = 0; i < scale_filters.size(); i++)
          scale_filters[i]->scale(q, l);
        q.glScaled();
      }
      else
        p.glScaled();
    }

    void begin(GLenum what) {
      fc.what = what;
      glBegin(what);
    }

    void end() {
      glEnd();
      fc.what = -1;
    }

    void end_of_face() {
      fc.face_i ++;
    }

    void end_of_particle() {
      fc.particle_i ++;
    }

    void end_of_cloud() {
      fc.face_i ++;
    }

    void end_of_frame() {
      fc.frame_i ++;
      fc.clear();
    }
};



class DrawAs {
  public:
    double alpha;
    
    DrawAs() {
      alpha = 1;
    }

    virtual void draw(vector<Point> &points, DrawBank &d) = 0;
};


class AsLines : public DrawAs {
  public:
    AsLines() : DrawAs() {}
    virtual void draw(vector<Point> &points, DrawBank &d)
    {
      int l;

      d.begin(GL_LINES);

      l = points.size();
      for (int i = 1; i < l; i++) {
        d.color_point(points[i-1]);
        d.color_point(points[i]);
        d.end_of_face();
      }
      d.color_point(points.back());
      d.color_point(points.front());
      d.end_of_face();

      d.end();
    }
};

class AsTriangles : public DrawAs {
  public:
    virtual void draw(vector<Point> &points, DrawBank &d)
    {
      int l;
      d.begin(GL_TRIANGLES);

      l = points.size();
      for (int i = 2; i < l; i++) {
        d.color_point(points[i-2]);
        d.color_point(points[i-1]);
        d.color_point(points[i]);
        d.end_of_face();
      }

      d.end();
    }
};

class AsTets : public DrawAs {
  public:
    virtual void draw(vector<Point> &points, DrawBank &d)
    {
      int l;
#if 1
      d.begin(GL_TRIANGLES);

      l = points.size();
      for (int i = 3; i < l; i += 2) {
        d.color_point(points[i-3]);
        d.color_point(points[i-2]);
        d.color_point(points[i-1]);
        d.end_of_face();

        d.color_point(points[i-3]);
        d.color_point(points[i-2]);
        d.color_point(points[i]);
        d.end_of_face();

        d.color_point(points[i-3]);
        d.color_point(points[i-1]);
        d.color_point(points[i]);
        d.end_of_face();

        d.color_point(points[i-2]);
        d.color_point(points[i-1]);
        d.color_point(points[i]);
        d.end_of_face();
      }
#else

      glBegin(GL_TRIANGLE_STRIP);

      l = points.size();
      for (int i = 3; i < l; i += 2) {
        points[i-3].draw(alpha);
        points[i-2].draw(alpha);
        points[i-1].draw(alpha);

        points[i].draw(alpha);

        points[i-3].draw(alpha);

        points[i-2].draw(alpha);

      }
#endif

      glEnd();
    }
};


class AsPoly : public DrawAs {
  public:
    AsPoly() : DrawAs() {}
    virtual void draw(vector<Point> &points, DrawBank &d)
    {
      int l;
      d.begin(GL_POLYGON);

      l = points.size();
      for (int i = 0; i < l; i ++) {
        d.color_point(points[i]);
      }
      d.end_of_face();

      d.end();
    }
};

class Placed {
  public:
    Pt pos;
    Pt dir;
    Pt scale;
    level_e level;

    Placed() {
      pos.set(0, 0, 0);
      dir.set(0, 0, 0);
      scale.set(1, 1, 1);
      level = level_unknown;
    }

    void placement(DrawBank &d) {
      d.translate(pos, level);
      d.rotate(dir, level);
      d.scale(scale, level);
    }

    void draw(DrawAs &as, DrawBank &d) { 
      glPushMatrix();
      placement(d);

      _draw(as, d);

      glPopMatrix();
    }

    virtual void _draw(DrawAs &as, DrawBank &d) = 0;

};


class Particle : public Placed {
  public:
    vector<Point> points;

    Particle() {
      level = level_particle;
    }

    virtual void _draw(DrawAs &as, DrawBank &d) {
      as.draw(points, d);
      d.end_of_particle();
    }

    Point &add_point() {
      points.resize(points.size() + 1);
      return points.back();
    }
};


class Cloud : public Placed {
  public:
    vector<Particle> particles;

    Cloud()
    {
      level = level_cloud;
    }

    Particle &add_particle() {
      particles.resize( particles.size() + 1 );
      return particles.back();
    }

    virtual void _draw(DrawAs &as, DrawBank &d) { 
      for (int i = 0; i < particles.size(); i++) {
        particles[i].draw(as, d);
      }
      d.end_of_cloud();
    }

    void step() {
    }

};

class ParticleGenesis {
  public:
    virtual void generate(Cloud &in) = 0;
};

class PointGenesis {
  public:
    virtual void generate(Particle &in) = 0;
};

class RandomPoints : public PointGenesis {
  public:
    int n;
    RandomPoints() {
      n = 4;
    }

    virtual void generate(Particle &in) {
      for (int i = 0; i < n; i++) {
        Point &p = in.add_point();
        p.random();
      }
    }
};

class RandomParticles : public ParticleGenesis {
  public:
    int n;
    double pos_range;
    double scale_min, scale_max;
    Pt dir_range;
    PointGenesis *point_genesis;

    RandomParticles() {
      defaults();
      point_genesis = NULL;
    }

    RandomParticles(PointGenesis *point_genesis) {
      defaults();
      this->point_genesis = point_genesis;
    }

    void defaults() {
      n = 50;
      pos_range = 10.;
      scale_min = -1.;
      scale_max = 1.;
      dir_range.set(360., 360., 360.);
    }

    virtual void generate(Cloud &in) {
      for (int i = 0; i < n; i++) {
        Particle &p = in.add_particle();
        p.pos.random();
        p.pos *= pos_range;

        p.scale.random();
        p.scale /= scale_max - scale_min;
        p.scale += scale_min;

        p.dir.random();
        p.dir.scale3(dir_range);

        if (point_genesis) {
          point_genesis->generate(p);
        }

      }
    }
};


class Quake : public TranslateFilter {
  public:
    double strength;
    double decay;
    level_e on_level;
    int seen_frame;

    Quake() {
      strength = .5;
      decay = 0.1;
      on_level = level_particle;
      seen_frame = 0;
    }

    virtual void translate(Pt &p, level_e l) {
      if (seen_frame != fc->frame_i) {
        seen_frame = fc->frame_i;
        strength *= 1. - min(1., max(0., decay));
      }

      if ((strength > 1e-4) && (l == on_level)) {
        Pt q;
        q.random();
        q *= strength;

        p += q;
      }
    };
};

class Scale : public ScaleFilter {
  public:
    Pt factor;
    level_e on_level;

    Scale() {
      on_level = level_particle;
      factor.set(1, 1, 1);
    }

    virtual void scale(Pt &p, level_e l){
      if (l == on_level) {
        p.scale3(factor);
      }
    }
};

class Explode : public PointFilter {
  public:
    double speed;
    double decay;
    int spread;
    int _faces;
    bool _reverse;

    vector<Pt> sum;
    vector<double> speeds;

    Explode() {
      speed = .1;
      decay = .97;
      spread = 10;
      _faces = 0;
      _reverse = false;
    }

    virtual void point(Pt &p) {
      int face_i = fc->face_i;

      if (! face_i)
        _faces += spread;

      if (face_i > _faces)
        return;

      if (_reverse) {
        if (face_i >= sum.size())
          return;
        sum[face_i] /= 1. + speed;
      }
      else {
        int have = sum.size();
        if (sum.size() < (face_i+1)) {
          sum.resize(face_i + 1);
          speeds.resize(face_i + 1);
          for (int i = have; i < sum.size(); i++) {
            sum[i].random();
            speeds[i] = speed;
          }
        }

        sum[face_i] *= 1. + speeds[face_i];
        speeds[face_i] *= decay;
      }
      p += sum[face_i];
    }

    void reverse() {
      for (int i = 0; i < speeds.size(); i++) {
        speeds[i] = - speed;
      }
      _reverse = true;
    }

    void clear() {
      sum.clear();
      speeds.clear();
      _faces = 0;
      _reverse = false;
    }
};

class ChangeColor : public ColorFilter {
  public:
    int n_cycle;
    palette_t *pal;

    ChangeColor() {
      n_cycle = 50;
      pal = NULL;
    }

    virtual void color(Color &c) {
      double col = (double)fc->color_i / n_cycle;
      c.set( pal->colors[(int)(col * pal->len) & PALETTE_LEN_MASK] );
    }
};

class StrobeCloud : public ColorFilter {
  public:
    Color dark;
    Color bright;
    int n_cycle;
    int n_flash;
    int n_blend;
    int which_mask;
    bool trigger;

    StrobeCloud() {
      dark.set(0, 0.1, 1, 0.01);
      bright.set(1, 1, 1, 1);
      which_mask = 0xff;
      n_cycle = 6;
      n_flash = 2;
      n_blend = 13;
    }

    int seen_frame;
    int which;
    int cycle;
    double fcycle;

    virtual void color(Color &c) {
      if (seen_frame != fc->frame_i) {
        seen_frame = fc->frame_i;

        if (cycle >= n_cycle) {
          cycle = 0;
          which = random() & which_mask;
        }
        else {
          cycle ++;
        }
        fcycle = max(0., 1. - ((double)cycle / n_flash));
      }

      double a = .1;
      //c = dark;
      //a = dark.a;

      int d = abs((which) - (fc->particle_i & which_mask));
      
      if (d <= n_blend) {
        double blend = 1. - (double)d / n_blend;
        blend *= blend;
        a += (1. - a) * fcycle * blend;
        //c.blend(bright, fcycle * blend);
      }

      c.a *= a;
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

    vector<Cloud> clouds;

    AsTriangles as_triangles;
    AsLines as_lines;
    AsPoly as_poly;
    AsTets as_tets;

    DrawBank bank;
    Quake quake;
    Explode explode;
    Scale particle_scale;
    ChangeColor change_color;
    StrobeCloud strobe_cloud;


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

      clouds.resize(1);
      {
        Cloud &c = clouds[0];

        RandomPoints rpo;
        RandomParticles rpa(&rpo);
        rpa.n = 100;
        rpa.scale_min = 10;
        rpa.scale_max = 20;

        rpa.generate(c);
        c.pos.set(-10, 0, 0);

        clouds.push_back(c);
        clouds[1].pos.set(10, 0, 0);
        clouds[1].scale.set(-1, 1, 1);
      }

      bank.add(quake);
      bank.add(particle_scale);

      change_color.pal = &palette;
      bank.add(change_color);
      bank.add(strobe_cloud);
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

      for (int i = 0; i < clouds.size(); i++) {
        clouds[i].step();
      }

      palette_selected = -1;

      clouds[0].dir.y -= .3;
      clouds[1].dir.y += .3;
    }


    void draw() {
      for (int i = 0; i < clouds.size(); i++) {
        clouds[i].draw(as_tets, bank);
      }
      bank.end_of_frame();
    }
};

Animation animation;




void draw_scene()
{
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity( );

  double vision = animation.vision;
  gluLookAt(0,3,20 + vision*200,0,0,vision * vision * 60,0,0,1);

  double sc = animation.lines_scale;

  double angleZ = animation.rot_z;
  double angleX = animation.rot_x;

  glScaled(sc, sc, sc);
  glRotated(angleZ,0,1,0);
  glRotated(angleX,1,0,0);


  animation.draw();

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
          animation.particle_scale.factor.y = 1.1 + axis_val;
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

    case 0:
      if (down) {
        switch(button) {
          case 2:
            animation.quake.strength = 1.;
            break;
          case 0:
            animation.explode.clear();
            animation.bank.remove(animation.explode);
            animation.bank.add(animation.explode);
            break;
          case 1:
            animation.explode.reverse();
            break;
          case 3:
            animation.bank.remove(animation.explode);
            break;
        }
      }
      break;

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

