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

double frandom(void) {
  return (double)(random()) / INT_MAX;
}

double frandom(double _min, double _max) {
  return _min + (_max - _min) * frandom();
}

/* vector<Moo> goo;
 * foreach(g, goo) {
 *   if (g->pritzical)
 *     g->frobnicate();
 * }
 */
#define foreach(I,VECTOR) \
  for (decltype(&*VECTOR.begin()) _i=0, I=&(VECTOR)[0]; \
       ((unsigned long int)_i) < (VECTOR).size(); \
       (*((unsigned long int*)&_i))++,I=&(VECTOR)[(unsigned long int)_i])

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

    void load(const double coords[3]) {
      x = coords[0];
      y = coords[1];
      z = coords[2];
    }

    static Pt random(double min_l=-1, double max_l=1)
    {
      return Pt(min_l + (max_l - min_l) * frandom(),
                min_l + (max_l - min_l) * frandom(),
                min_l + (max_l - min_l) * frandom());
    }

    bool zero() const {
      return (!x) && (!y) && (!z);
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

    Pt &limit(double c_min, double c_max) {
      x = min(c_max, max(c_min, x));
      y = min(c_max, max(c_min, y));
      z = min(c_max, max(c_min, z));
      return *this;
    }

    int wrap(const Pt &center, double dist) {
      Pt d = (center - (*this));
      if (d.len() > dist) {
        *this += d.unit() * (dist*2);
        return 1;
      }
      return 0;
    }

    double dot(const Pt &p) const
    {
      return
        x * p.x + y * p.y + z * p.z;
    }

    double udot(const Pt &p) const
    {
      Pt u = unit();
      Pt pu = p.unit();
      return
        min(1.0, max(-1.0, u.x * pu.x + u.y * pu.y + u.z * pu.z));
    }

    Pt cross(const Pt &p) const
    {
      return Pt(y * p.z - z * p.y,
                z * p.x - x * p.z,
                x * p.y - y * p.x);
    }

    double min_angle(const Pt &p) const
    {
      return acos(this->udot(p));
    }

    double angle(const Pt &p, const Pt &n) const
    {
      double ma = acos(this->udot(p));
      if (n.udot( this->cross(p) ) < 0)
        return - ma;
      return ma;
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

    Pt operator+(const double f) const {
      Pt x(*this);
      x += f;
      return x;
    }

    Pt operator-(const double f) const {
      Pt x(*this);
      x -= f;
      return x;
    }

    bool operator==(const Pt& other) const {
      return (x == other.x) && (y == other.y) && (z == other.z);
    }

    bool operator!=(const Pt& other) const {
      return (x != other.x) || (y != other.y) || (z != other.z);
    }

    void glVertex3d() const
    {
      ::glVertex3d(x, y, z);
    }

    void glTranslated() const
    {
      ::glTranslated(x, y, z);
    }

    void glRotated() const
    {
      // last, turn left/right
      if (y) {
        ::glRotated(y * (180./M_PI), 0,1,0);
      }
      // second, turn nose up or down
      if (x) {
        ::glRotated(x * (180./M_PI), 1,0,0);
      }
      // first, roll left or right
      if (z) {
        ::glRotated(z * (180./M_PI), 0,0,1);
      }
    }

    void glScaled() const
    {
      ::glScaled(x, y, z);
    }

    void print() const
    {
      printf("x%5.2f y%5.2f z%5.2f\n", x, y, z);
    }

    Pt &operator=(double v)
    {
      x = y = z = v;
      return *this;
    }

    Pt scaled(const Pt &factors) const
    {
      return Pt(x * factors.x, y * factors.y, z * factors.z);
    }

    Pt without(const Pt &axis) const
    {
      return (*this) - project(axis);
    }

    Pt project(const Pt &axis, bool neg_means_zero=true) const
    {
      double d = dot(axis);
      if (neg_means_zero && (d < 0))
        return Pt();
      return axis * d;
    }
};

struct Matrix33 {
    double a, b, c,
           d, e, f,
           g, h, i;

    static Matrix33 from_rot3(Pt rot3) {
      double cx = cos(rot3.x);
      double cy = cos(rot3.y);
      double cz = cos(rot3.z);
      double sx = sin(rot3.x);
      double sy = sin(rot3.y);
      double sz = sin(rot3.z);
      Matrix33 mx = {
        1, 0, 0,
        0, cx, sx,
        0, -sx, cx
      };

      Matrix33 my {
        cy, 0, -sy,
        0, 1, 0,
        sy, 0, cy
      };

      Matrix33 mz {
        cz, sz, 0,
        -sz, cz, 0,
        0, 0, 1
      };

      return my * mx * mz;
    }

    Matrix33 operator*(const Matrix33 &o) const
    {
      return Matrix33{
        a*o.a + b*o.d + c*o.g, a*o.b + b*o.e + c*o.h, a*o.c + b*o.f + c*o.i,
        d*o.a + e*o.d + f*o.g, d*o.b + e*o.e + f*o.h, d*o.c + e*o.f + f*o.i,
        g*o.a + h*o.d + i*o.g, g*o.b + h*o.e + i*o.h, g*o.c + h*o.f + i*o.i
      };
    }

    Pt operator*(const Pt &p) const
    {
      /*     a b c   x    a*x + b*y + c*z
             d e f * y =  d*x + e*y + f*z
             g h i   z    g*x + h*y + i*z */
      return Pt( a*p.x + b*p.y + c*p.z,
                 d*p.x + e*p.y + f*p.z,
                 g*p.x + h*p.y + h*p.z);
    }

    double det() const
    {
      return a*e*i + b*f*g + c*d*h - g*e*c - h*f*a - i*d*b;
    }

    Matrix33 inverse() const
    {
      double dd = det();
      if (fabs(dd) < 1e-6)
        return Matrix33();
      return Matrix33{
        e*i - f*h, c*h - b*i, b*f - c*e,
        f*g - d*i, a*i - c*g, c*d - a*f, 
        d*h - e*g, b*g - a*h, a*e - b*d
      } / dd;
    }

    Matrix33 operator/(double x) const
    {
      return Matrix33{
          a/x, b/x, c/x,
          d/x, e/x, f/x,
          g/x, h/x, i/x
      };
    }

};

class Plane {
  public:
    Pt pos;
    Pt n;

    Plane(const Pt &pos, const Pt &n)
    {
      this->pos = pos;
      this->n = n;
    }

    Plane(const Pt &p0, const Pt &p1, const Pt &p2)
    {
      pos = p0;
      n = (p1 - p0).cross(p2 - p0);
    }
};

class Line {
  public:
    Pt pos;
    Pt len;

    Line(const Pt &pos, const Pt &len)
    {
      this->pos = pos;
      this->len = len;
    }
};

bool intersect(const Plane &plane, const Line &line,
               double *at=NULL, Pt *intersection=NULL)
{
  double D = plane.n.dot(line.len);
  if (fabs(D) < 1e-6)
    return false;

  double N = plane.n.dot(plane.pos - line.pos);
  double _at = N / D;

  if (at)
    *at = _at;

  if ((_at < -1e-6) || (_at > (1. + 1e-6)))
    return false;

  if (intersection)
    *intersection = line.pos + line.len * _at;

  return true;
}

bool intersect(const Pt &p0, const Pt &p1, const Pt &p2, Line &line,
               double *at=NULL, Pt *intersection=NULL,
               bool parallelogram=false)
{
  double _at;
  Pt _intersection;

  Pt u = p1 - p0;
  Pt v = p2 - p0;

  if (! intersect(Plane(p0, u.cross(v)), line, &_at, &_intersection))
    return false;

  if (at)
    *at = _at;

  if (intersection)
    *intersection = _intersection; 

  Pt w = _intersection - p0;

  double udotv = u.dot(v);
  double wdotv = w.dot(v);
  double vdotv = v.dot(v);
  double wdotu = w.dot(u);
  double udotu = u.dot(u);

  double D = udotv * udotv - udotu * vdotv;

  double uu = ( udotv * wdotv - vdotv * wdotu ) / D;
  if ((uu < -1e-6) || (uu > (1. + 1e-6)))
    return false;

  double vv = ( udotv * wdotu - udotu * wdotv ) / D;
  if ((vv < -1e-6) || (vv > (1. + 1e-6)))
    return false;

  return parallelogram || ((vv + uu) < (1. + 1e-6));
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

    void set(rgb_t &rgb) {
      this->r = rgb.r;
      this->g = rgb.g;
      this->b = rgb.b;
    }

    void load(const double color[4]) {
      this->r = color[0];
      this->g = color[1];
      this->b = color[2];
      this->a = color[3];
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

    Color& operator=(const Pt &p) {
      r = p.x;
      g = p.y;
      b = p.z;
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

    void begin() {
      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D, id);
    }

    void end() {
      glDisable(GL_TEXTURE_2D);
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

/* Laugh if you might, yet this gives a nice overview of drawing primitives. */
class Draw {
  public:
    static void point(Pt &p) {
      p.glVertex3d();
    }

    static void color(Color &c) {
      c.glColor();
    }

    static void color_point(Point &p) {
      color(p.c);
      point(p);
    }

    static void texture_coord(double cx, double cy) {
      ::glTexCoord2d(cx, cy);
    }

    static void translate(Pt &p) {
      p.glTranslated();
    }

    static void rotate(Pt &p) {
      p.glRotated();
    }

    static void scale(Pt &p) {
      p.glScaled();
    }

    static void begin(GLenum what) {
      glBegin(what);
    }

    static void end() {
      glEnd();
    }
};

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
      texture->begin();
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
      texture->end();
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
    Revolve(double max_speed=(M_PI/5)) {
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

    operator double() const {
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

    Pt impulse() const
    {
      return v * m;
    }
};


/**
adapted from http://www.euclideanspace.com/physics/dynamics/collision/threed/index.htm

This function calulates the velocities after a 3D collision vaf, vbf, waf and wbf from information about the colliding bodies
@param double e coefficient of restitution which depends on the nature of the two colliding materials
@param double ma total mass of body a
@param double mb total mass of body b
@param matrix Ia inertia tensor for body a in absolute coordinates (if this is known in local body coordinates it must
                 be converted before this is called).
@param matrix Ib inertia tensor for body b in absolute coordinates (if this is known in local body coordinates it must
                 be converted before this is called).
@param vector ra position of collision point relative to centre of mass of body a in absolute coordinates (if this is
                 known in local body coordinates it must be converted before this is called).
@param vector rb position of collision point relative to centre of mass of body b in absolute coordinates (if this is
                 known in local body coordinates it must be converted before this is called).
@param vector n normal to collision point, the line along which the impulse acts.
@param vector vai initial velocity of centre of mass on object a
@param vector vbi initial velocity of centre of mass on object b
@param vector wai initial angular velocity of object a
@param vector wbi initial angular velocity of object b
@param vector vaf final velocity of centre of mass on object a
@param vector vbf final velocity of centre of mass on object a
@param vector waf final angular velocity of object a
@param vector wbf final angular velocity of object b
*/
void collide(double restitution,
             double ma,
             double mb,
             const Matrix33 &Ia,
             const Matrix33 &Ib,
             const Pt &ra,
             const Pt &rb,
             const Pt &n,
             const Pt &vai,
             const Pt &vbi,
             const Pt &wai,
             const Pt &wbi,
             Pt &vaf,
             Pt &vbf,
             Pt &waf,
             Pt &wbf)
{
  Matrix33 inv_Ia = Ia.inverse();
  Matrix33 inv_Ib = Ib.inverse();

  Pt anga_a = inv_Ia * n.cross(ra);
  Pt angv_a = anga_a.cross(ra);  // calculate the linear velocity of collision point on a due to rotation of a

  Pt anga_b = inv_Ib * n.cross(rb);
  Pt angv_b = anga_b.cross(rb);  // calculate the linear velocity of collision point on b due to rotation of b

  double scalar = 1/ma + angv_a.dot(n) + 1/mb + angv_b.dot(n);

  double Jmod = (restitution+1) * (vai-vbi).len() / scalar;
  Pt J = n * Jmod;
  vaf = vai - J * (1/ma);
  vbf = vbi - J * (1/mb);
  waf = wai - anga_a;
  wbf = wbi - anga_b;
}
