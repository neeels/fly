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

#include "palettes.h"
#include "textures.h"

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

    bool zero() const {
      return (!x) && (!y) && (!z);
    }

    void ang2cart_xy(double ax, double ay, double r) {
      set(r * cos(ax) * sin(ay),
          r * sin(ax),
          r * cos(ax) * sin(ay)
          );
    }

    Pt ang2cart_xy() const {
      Pt p;
      p.ang2cart_xy(x, y, z);
      return p;
    }

    Pt cart2ang(const Pt &top) const {
      double l = len();
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

      Pt u = unit();

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
      topu -= u * topu.project(u);

      /* Angle of the "zero" top to the plane-ized and unit-ized user supplied
       * top, angle sign relative to this vector. */
      double az = top_zero.angle(topu, u);

      /* I measures the angle from the z axis rightwards. That's counter the
       * canonical angle direction OpenGL uses, so let's correct for that... */
      return Pt(M_PI - ax, ay, - az);
    }


    void ang2cart_yz(double r, double ay, double az) {
      set(r * cos(ay) * cos(az),
          r * cos(ay) * sin(az),
          r * sin(ay));
    }

    Pt ang2cart_yz() const {
      Pt p;
      p.ang2cart_yz(x, y, z);
      return p;
    }

    Pt cart2ang_yz() const {
      // dunno
      double r = len();
      if (fabs(r) < 1e-5) {
        return Pt();
      }

      double ay = acos(z / r);
      double az;
      if (fabs(x) < 1e-5) {
        az = 0;
      }
      else {
        az = atan(y / x);
      }

      return Pt(r, ay, az);
    }

    double len() const {
      return sqrt(x*x + y*y + z*z);
    }

    Pt unit() const {
      double l = len();
      if (l < 1.e-6)
        return Pt();
      return (*this) / l;
    }

    void rot_about(const Pt &axis, double rad) {
      double xx, yy, zz;
      double l = axis.len();
      double u = axis.x / l;
      double v = axis.y / l;
      double w = axis.z / l;
      double cos_rad = cos(rad);
      double sin_rad = sin(rad);
      double f1 = (u*x + v*y + w*z)*(-cos_rad + 1);
      xx = u*f1 + x*cos_rad + (- w*y + v*z)*sin_rad;
      yy = v*f1 + y*cos_rad + (  w*x - u*z)*sin_rad;
      zz = w*f1 + z*cos_rad + (- v*x + u*y)*sin_rad;
      x = xx;
      y = yy;
      z = zz;
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

    double dot(const Pt &p) const
    {
      return
        x * p.x + y * p.y + z * p.z;
    }

    Pt cross(const Pt &p) const
    {
      return Pt(y * p.z - z * p.y,
                z * p.x - x * p.z,
                x * p.y - y * p.x);
    }

    double min_angle(const Pt &p) const
    {
      Pt a = unit();
      Pt b = p.unit();
      return acos(a.dot(b));
    }

    double angle(const Pt &p, const Pt &n) const
    {
      Pt a = unit();
      Pt b = p.unit();
      Pt nu = n.unit();
      double ma = acos(a.dot(b));
      if (nu.dot( a.cross(b) ) < 0)
        return - ma;
      return ma;
    }

    double project(const Pt &p) const
    {
      double d = dot(p);
      double sign = (d >= 0? 1. : -1.);
      return sign * sqrt(fabs(d));
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

    Pt operator+(const Pt &p) const {
      Pt x(*this);
      x += p;
      return x;
    }

    Pt operator-(const Pt &p) const {
      Pt x(*this);
      x -= p;
      return x;
    }

    Pt operator*(const double f) const {
      Pt x(*this);
      x *= f;
      return x;
    }

    Pt operator/(const double f) const {
      Pt x(*this);
      x /= f;
      return x;
    }


    void glVertex3d() {
      ::glVertex3d(x, y, z);
    }

    void glTranslated() {
      ::glTranslated(x, y, z);
    }

    void glRotated() {
      if (y) {
        ::glRotated(y, 0, 1, 0); //cos(x)*M_PI/180, sin(x*M_PI/180));
      }
      if (x) {
        ::glRotated(x,1,0,0);
      }
      if (z) {
        ::glRotated(z,0,0,1);
      }
    }

    void glScaled() {
      ::glScaled(x, y, z);
    }

    void print() {
      printf("x%5.2f y%5.2f z%5.2f\n", x, y, z);
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


class Texture {
  public:
    GLuint id;

    Texture() {}

    void load(const char *path) {
      id = load_texture(path, false);
      printf("Loaded texture %x: %s\n", (int)id, path);
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
      else {
        p.glVertex3d();
      }
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

    void texture_coord(double cx, double cy) {
      ::glTexCoord2d(cx, cy);
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

      glLineWidth(2);
      glEnable( GL_LINE_SMOOTH );
      //glEnable( GL_POLYGON_SMOOTH );
      glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
      //glHint( GL_POLYGON_SMOOTH_HINT, GL_NICEST );

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

class AsPoints : public DrawAs {
  public:
    AsPoints() : DrawAs() {}
    virtual void draw(vector<Point> &points, DrawBank &d)
    {
      int l;

      d.begin(GL_POINTS);

      l = points.size();
      for (int i = 0; i < l; i++) {
        d.color_point(points[i]);
        d.end_of_face();
      }

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
    int pitch;

    AsTets() {
      pitch = 2;
    }

    virtual void draw(vector<Point> &points, DrawBank &d)
    {
      int l;
#if 1
      d.begin(GL_TRIANGLES);

      l = points.size();
      for (int i = 3; i < l; i += pitch) {
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

class AsQuads : public DrawAs {
  public:
    int pitch;

    AsQuads(){
      pitch = 4;
    }

    virtual void draw(vector<Point> &points, DrawBank &d)
    {
      int l;
      d.begin(GL_QUADS);

      l = points.size();
      for (int i = 3; i < l; i += pitch) {
        d.color_point(points[i-3]);
        d.color_point(points[i-2]);
        d.color_point(points[i-1]);
        d.color_point(points[i]);
        d.end_of_face();
      }

      glEnd();
    }
};

class AsTexturePlanes : public DrawAs {
  public:
    int pitch;
    Texture *texture;

    AsTexturePlanes() {
      pitch = 4;
      texture = NULL;
    }

    virtual void draw(vector<Point> &points, DrawBank &d)
    {
      int l;
      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D, texture->id);
      d.begin(GL_QUADS);

      l = points.size();
      for (int i = 3; i < l; i += pitch) {
        d.texture_coord(0, 0);
        d.color_point(points[i-3]);
        d.texture_coord(1, 0);
        d.color_point(points[i-2]);
        d.texture_coord(1, 1);
        d.color_point(points[i-1]);
        d.texture_coord(0, 1);
        d.color_point(points[i]);
        d.end_of_face();
      }

      glEnd();
      glDisable(GL_TEXTURE_2D);
    }
};



class Placed {
  public:
    Pt pos;
    Pt rot3;
    Pt scale;
    level_e level;

    Placed() {
      pos.set(0, 0, 0);
      rot3.set(0, 0, 0);
      scale.set(1, 1, 1);
      level = level_unknown;
    }

    void placement(DrawBank &d) {
      d.translate(pos, level);
      d.rotate(rot3, level);
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
        p /= 2;
      }
    }
};

class Block : public PointGenesis {
  public:
    Block() {
    }

    virtual void generate(Particle &in) {
#define P(x,y,z) {\
          in.add_point().set(x, y, z); \
        }

#define A  P(-.5, -.5, -.5);
#define B  P(-.5, -.5,  .5);
#define C  P(-.5,  .5,  .5);
#define D  P(-.5,  .5, -.5);
#define E  P( .5, -.5, -.5);
#define F  P( .5, -.5,  .5);
#define G  P( .5,  .5,  .5);
#define H  P( .5,  .5, -.5);

      A B C D
      A B F E
      B C G F
      E F G H
      C D H G
      D A E H

#undef A
#undef B
#undef C
#undef D
#undef E
#undef F
#undef G
#undef H
#undef P
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
        p.scale *= scale_max - scale_min;
        p.scale += scale_min;

        p.rot3.random();
        p.rot3.scale3(dir_range);

        if (point_genesis) {
          point_genesis->generate(p);
        }

      }
    }
};

#include "font.h"

class WriteInBlocks : public ParticleGenesis {
  public:
    const char *text;
    Pt block_pitch;
    Pt block_rel_size;
    PointGenesis *point_genesis;

    WriteInBlocks() {
      defaults();
    }

    WriteInBlocks(const char *text) {
      defaults();
      this->text = text;
    }

    void defaults() {
      text = "written";
      block_pitch.set(.5, .5, .5);
      block_rel_size.set(.95, .95, .95);
    }

    virtual void generate(Cloud &in) {
      size_t l = strlen(text);

      Pt letter_size(5, 8, 1);
      letter_size.scale3(block_pitch);

      Pt letter_mid = letter_size;
      letter_mid /= 2;
      letter_mid.y = - letter_mid.y;

      Pt letter_direction(letter_size.x + block_pitch.x/5, 0, 0);

      Pt letter_pos = letter_direction;
      letter_pos /= -2;
      letter_pos *= l - 1;

      Pt block_pos;

      Pt block_size = block_pitch;
      block_size.scale3(block_rel_size);

      for (size_t i = 0; i < l; i++) {

        block_pos = letter_pos;
        block_pos -= letter_mid;

        const uint8_t *bits = font[ text[i] & 0x7f ];


        Pt y0_pos = block_pos;
        for (int x = 0; x < 5; x++) {
          block_pos = y0_pos;
          uint8_t b = bits[x];
          for (int y = 0; y < 8; y++) {
            if (b & 1) {
              Particle &p = in.add_particle();
              point_genesis->generate(p);

              p.pos = block_pos;
              p.scale = block_size;
            }
            b >>= 1;
            block_pos.y -= block_pitch.y;
          }
          y0_pos.x += block_pitch.x;
        }

        letter_pos += letter_direction;
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

class PointQuake : public PointFilter {
  public:
    double strength;
    double decay;
    level_e on_level;
    int seen_frame;

    PointQuake() {
      strength = .5;
      decay = 0.1;
      on_level = level_cloud;
      seen_frame = 0;
    }

    virtual void point(Pt &p) {
      if (seen_frame != fc->frame_i) {
        seen_frame = fc->frame_i;
        strength *= 1. - min(1., max(0., decay));
      }

      if (strength > 1e-4) {
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
        speeds[face_i] /= decay;
        sum[face_i] /= 1. + speeds[face_i];
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
        speeds[i] = speed * 5;
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
      n_cycle = 3000;
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
    int _n_cycle;

    StrobeCloud() {
      dark.set(0, 0.1, 1, 0.01);
      bright.set(1, 1, 1, 1);
      which_mask = 0xff;
      _n_cycle = n_cycle = 6;
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

        if (cycle >= _n_cycle) {
          cycle = 0;
          _n_cycle = 1 + frandom() * n_cycle;
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


class Revolve : public RotateFilter {
  public:
    Revolve(double max_speed=10) {
      this->max_speed = max_speed;
    }

    double max_speed;
    vector<Pt> state;
    vector<Pt> sum;

    virtual void rotate(Pt &p, level_e l){
      if (l == level_particle) {

        int have = state.size();
        int particle_i = fc->particle_i;
        if (state.size() < (particle_i+1)) {
          state.resize(particle_i + 1);
          sum.resize(particle_i + 1);
          for (int i = have; i < state.size(); i++) {
            state[i].random();
            state[i] *= max_speed * frandom();
            sum[i].set(0, 0, 0);
          }
        }

        sum[particle_i] += state[particle_i];
        p += sum[particle_i];
      }
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

class Mass {
  public:
    double m;
    Pt v;

    Mass() 
      :m(1), v(0, 0, 0)
    {}

    void accelerate(const Pt &F, double dt)
    {
      v += F * (dt / m);
    }
};


