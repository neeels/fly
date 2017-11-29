
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
#include "draw.h"

#define Pf(V) printf(#V "=%f\n", (float)V)
#define Pi(V) printf(#V "=%d\n", (int)V)

SDL_Window *window = NULL;

void draw_scene();

bool running = true;
volatile int frames_rendered = 0;
volatile int avg_frame_period = 0;
#define AVG_SHIFTING 3
float want_fps = 25;

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

		Textures textures;
    Texture *dot;

    vector<Cloud> clouds;

    AsTriangles as_triangles;
    AsLines as_lines;
    AsPoly as_poly;
    AsTets as_tets;
    AsQuads as_quads;
    AsTexturePlanes as_texture_planes;
    AsPoints as_points;

    DrawBank bank;
    Quake quake;
    PointQuake point_quake;
    Explode explode;
    Scale particle_scale;
    ChangeColor change_color;
    StrobeCloud strobe_cloud;
    Revolve revolve;


    Animation() : textures(1) {
      is_palette = &palettes[0];
      want_palette = is_palette;

      make_palettes();
      make_palette(&palette, PALETTE_LEN,
                   palette_defs[0]);

      make_palette(&blended_palette, PALETTE_LEN,
                   palette_defs[0]);

      dot = textures.load("images/dot_white_200.png");
      //dot.load("images/metal091.jpg");
      as_texture_planes.texture = dot;

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

        if (0) {
          RandomPoints rpo;
          RandomParticles rpa(&rpo);
          rpa.n = 100;
          rpa.scale_min = 10;
          rpa.scale_max = 20;

          rpa.generate(c);
        }
        if (0)
        {
          Particle p;
          p.points.resize(4);
          p.points[0].set(-1, -1, 0);
          p.points[3].set( 1, -1, 0);
          p.points[2].set( 1,  1, 0);
          p.points[1].set(-1,  1, 0);
          p.pos.set(0, 0, -0.5);

          c.particles.push_back(p);

          p.pos.set(0, 0, 0.5);
          c.particles.push_back(p);

          int i;
          for (i = 0; i < c.particles.size(); i++) {
            printf("%p ", &(c.particles[i]));
            c.particles[i].pos.print();
            p.points[3].print();
            printf(" %p\n", &(c.particles[i].pos));
          }
          c.particles[1].pos.set(0, 0, .6);
          for (i = 0; i < c.particles.size(); i++) {
            printf("%p ", &(c.particles[i]));
            c.particles[i].pos.print();
            p.points[3].print();
            printf(" %p\n", &(c.particles[i].pos));
          }
        }
        if (0)
        {
          Particle &p = c.add_particle();
          p.points.resize(4);
          p.points[0].set(-1, -1, 0);
          p.points[1].set( 1, -1, 0);
          p.points[2].set( 1,  1, 0);
          p.points[3].set(-1,  1, 0);
          p.pos.set(0, 0, -0.5);
        }
        
        
        //RandomPoints b;
        //b.n = 16;
        Block b;
        WriteInBlocks say;
        say.block_pitch.set(1, 1, 1);
        say.point_genesis = &b;
        say.generate(c);
         
        c.pos.set(-0, 0, 0);

        /*
        clouds.push_back(c);
        clouds[1].pos.set(0, 0, 0);
        clouds[1].scale.set( 1, 1, 1);

        clouds.resize(3);
        clouds[2].pos.set(0, 0, 0);
        clouds[2].scale.set(1,  1, 1);
        say.generate(clouds[2]);
        */
      }

      bank.add(quake);
      bank.add(point_quake);
      bank.add(particle_scale);
      //particle_scale.factor.set(.1, .1, .1);

      change_color.pal = &palette;
      bank.add(change_color);
    }

    void load() {
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


      as_quads.pitch = 4;
      as_tets.pitch = 2;
    }


    void draw() {
      for (int i = 0; i < clouds.size(); i++) {
        clouds[i].draw(as_quads, bank);
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
  SDL_GL_SwapWindow(window);
}


FILE *out_stream = NULL;
FILE *out_params = NULL;
FILE *in_params = NULL;

char *audio_path = NULL;

int W = 1400;//1920;
int H = 900;//1080;

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
          animation.rot_x.change = (axis_val*axis_val*axis_val) / 6;
          break;
        case 0:
          animation.rot_z.change = (axis_val*axis_val*axis_val) / 6;
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
          animation.particle_scale.factor += axis_val / 20;
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
          case 0:
            animation.quake.strength = .3;
            break;
          case 1:
            animation.point_quake.strength = .1;
            break;
          case 2:
            animation.bank.remove(animation.strobe_cloud);
            animation.bank.add(animation.strobe_cloud);
            break;
          case 3:
            animation.bank.remove(animation.strobe_cloud);
            break;
        }
      }
      break;

    case 1:
      if (down) {
        switch(button) {
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
      if (down) {
        switch(button) {
          case 0:
            animation.bank.remove(animation.revolve);
            animation.bank.add(animation.revolve);
            break;
          case 1:
            animation.bank.remove(animation.revolve);
            break;
        }
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

  window = SDL_CreateWindow("T",
                            SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                            W, H,
                            SDL_WINDOW_OPENGL);
  SDL_GLContext gl_ctx = SDL_GL_CreateContext(window);

  if (fullscreen)
    SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN);

  SDL_ShowCursor(SDL_DISABLE);

  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective(70,(double)W/H,1,100);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


  animation.load();

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

              case 'f':
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

