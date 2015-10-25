#define Pf(V) printf(#V "=%f\n", (float)V)
#define Pi(V) printf(#V "=%d\n", (int)V)

#include "ctrl_layers.h"
#include "draw.h"


void draw_scene();

bool running = true;
volatile int frames_rendered = 0;
volatile int avg_frame_period = 0;
#define AVG_SHIFTING 3
float want_fps = 25;

class World {
  public:
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

    Texture metal;

    World() {
      is_palette = &palettes[0];
      want_palette = is_palette;

      make_palettes();
      make_palette(&palette, PALETTE_LEN,
                   palette_defs[0]);

      make_palette(&blended_palette, PALETTE_LEN,
                   palette_defs[0]);

      as_texture_planes.texture = &metal;

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
         
        c.pos.set(0, 0, -20);

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

    }

    void load() {
      metal.load("images/metal091.jpg");
    }

    void step() {
      for (int i = 0; i < clouds.size(); i++) {
        clouds[i].step();
      }

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

struct Fly {
	Param x, y, z;
	Param rot_x, rot_z;
  Param velocity;
	Param vision;

	Fly() {
    x = 0;
    y = 0;
    z = 0;
    rot_x.slew = .94;
    rot_z.slew = .94;
    vision.slew = .96;
	}

  void step(void) {
      x.step();
      y.step();
      z.step();
      rot_x.step();
      rot_z.step();
      velocity.step();
      vision.step();

      //printf("%f\n", (float)velocity);
      if (fabs(velocity) > 1.e-5) {
        Pt dir;
        dir.ang2cart(velocity, rot_z + M_PI/2, rot_x);
        printf("%f  rx %f rz %f  ", (float)velocity, (float)rot_x, (float)rot_z); dir.print();
        x += dir.x;
        y += dir.y;
        z += dir.z;
      }
  }
};


World world;
Fly camera;


void draw_scene()
{
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity( );

  double vision = camera.vision;
  gluLookAt(0, 0, 0,
            0, 0, -2,
            0, 1, 0);

  double x = camera.x;
  double y = camera.y;
  double z = camera.z;
  double angleZ = (180./M_PI) * camera.rot_z;
  double angleX = (180./M_PI) * camera.rot_x;

  glRotated(angleZ,0,1,0);
  glRotated(angleX,1,0,0);
  glTranslated(x, y, z);

  world.draw();

  glFlush();
  SDL_GL_SwapBuffers();
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
  switch(axis)
  {
  case 1:
    camera.rot_x.change = (axis_val*axis_val*axis_val) / 10;
    break;
  case 0:
    camera.rot_z.change = (axis_val*axis_val*axis_val) / 10;
    break;
  case 3:
    //camera.x.change = - axis_val;
    break;
  case 4:
    camera.velocity = -axis_val / 10;
    break;
  case 5:
    camera.vision = (1. + axis_val)/2;
    break;
  case 2:
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
  gluPerspective(70,(double)W/H,1,100);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


  world.load();

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


    world.step();
    camera.step();
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

