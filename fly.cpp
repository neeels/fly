#define Pf(V) printf(#V "=%f\n", (float)V)
#define Pi(V) printf(#V "=%d\n", (int)V)

#include "ctrl_layers.h"
#include "draw.h"


void draw_scene();

bool running = true;
volatile int frames_rendered = 0;
volatile int avg_frame_period = 0;
volatile double dt;
#define AVG_SHIFTING 3
float want_fps = 25;

const double maneuver_engines_strength = 0.5;
const double forward_engines_strength = 2;

class Fly : public Particle {
  public:
  Pt nose;
  Pt top;
  Pt top_lag;
  Pt right;

  Param top_angle;
  Param roll_x;
  Param roll_y;
  Param roll_z;

  Mass mass;
  Param propulsion_forward;
  Param propulsion_break;

  Param wings;

  bool do_cruise;
  Param cruise_v;

	Fly() : Particle() {
    nose.set(0, 0, -1);
    top.set(0, 1, 0);
    top_lag = top;
    right.set(1, 0, 0);
    update_normals();
    roll_x.slew = .94;
    roll_y.slew = .94;
    roll_z.slew = .94;
    wings.limit(0, 1);
    wings = .1;

    do_cruise = true;
    cruise_v = .05;

    Block b;
    b.generate(*this);

    int l = points.size();
    for (int i = 0; i < l; i++) {
      Point &p = points[i];
      if (p.z < 0) {
        p.x *= .1;
      }
      p.y = max(min(p.y, .1), -.1);
    }
	}

  void update_normals() {
      nose = nose.unit();

      //top.set(nose.y, nose.z, nose.x);
      //top.rot_about(nose, top_angle);
      top = top.unit();

      //right.set(nose.z, nose.x, nose.y);
      //right.rot_about(nose, top_angle);
      right = right.unit();

      double is_ortho_top = nose.x * top.x + nose.y * top.y + nose.z * top.z;
      double is_ortho_right = nose.x * right.x + nose.y * right.y + nose.z * right.z;
      if ((fabs(is_ortho_top) > 1.e-6)
          || (fabs(is_ortho_right) > 1.e-6)) {
        printf("ORTHO! %g %g\n", is_ortho_top, is_ortho_right);
      }

      rot3 = nose.cart2ang(top) * (180./M_PI);
  }


  void step(void) {
      roll_x.step();
      roll_y.step();
      roll_z.step();
      top_angle.step();
      propulsion_forward.step();
      propulsion_break.step();
      wings.step();
      cruise_v.step();

      nose.rot_about(top, roll_y);
      right.rot_about(top, roll_y);

      nose.rot_about(right, roll_x);
      top.rot_about(right, roll_x);

      top.rot_about(nose, roll_z);
      right.rot_about(nose, roll_z);

      update_normals();

      const int div = 10;
      top_lag = (top_lag * (div-1) + top) / div;


      Pt wings_force;

      if (dt < 1.e-5)
        return;

      if (wings > 1.e-5) {
        Pt v_want = nose * mass.v.len();
        wings_force = (v_want - mass.v) * wings / dt;
      }

      Pt break_want = mass.v * (-propulsion_break / dt);

      /* main engines component of break force */
      double forward_drive_break_amount = break_want.project(nose);
      forward_drive_break_amount = min(forward_drive_break_amount, forward_engines_strength);
      forward_drive_break_amount = max(forward_drive_break_amount, 0.);

      /* remove main engines part from break force */
      break_want += nose * forward_drive_break_amount;

      /* remaining components with weaker engines */
      double break_drive = break_want.len();
      if (break_drive > maneuver_engines_strength) {
        break_want *= maneuver_engines_strength / break_drive;
        break_drive = maneuver_engines_strength;
      }
      Pt maneuver_engines = break_want;

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
            want_propulsion_forward = (diff / forward_engines_strength) / dt;
            if (want_propulsion_forward < propulsion_forward)
              want_propulsion_forward = propulsion_forward;
          }
        }
      }

      double forward_drive_amount = max(forward_drive_break_amount, want_propulsion_forward);
      Pt main_engine_force = nose * (forward_drive_amount * forward_engines_strength);

      Pt propulsion = main_engine_force + maneuver_engines;
      mass.accelerate(wings_force + propulsion, dt);

      pos += mass.v * dt;

      static int skip = 0;
      if ((skip ++) > 10) {
        skip = 0;
        printf("v=%5.2f p=%5.2f g=%5.2f break_fwd=%5.2f break_maneuver=%5.2f\n",
               mass.v.len(),
               propulsion.len(),
               wings_force.len(),
               forward_drive_break_amount,
               maneuver_engines.len());
      }
  }


  void draw(DrawBank &bank)
  {
    static AsQuads as_quads;
    Particle::draw(as_quads, bank);
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


class World {
  public:
    palette_t palette;

    vector<Cloud> clouds;

    AsTriangles as_triangles;
    AsLines as_lines;
    AsPoly as_poly;
    AsTets as_tets;
    AsQuads as_quads;
    AsTexturePlanes as_texture_planes;
    AsPoints as_points;

    DrawBank cloud_bank;
    DrawBank fixed_bank;
    Quake quake;
    PointQuake point_quake;
    Explode explode;
    Scale particle_scale;
    ChangeColor change_color;
    ChangeColor change_color2;
    StrobeCloud strobe_cloud;
    Revolve revolve;

    Texture metal;

    Fly fly;
    Camera cam;

    World() : revolve(1) {

      make_palettes();
      make_palette(&palette, PALETTE_LEN,
                   palette_defs[0]);

      as_texture_planes.texture = &metal;

      clouds.resize(2);
      {
        Cloud &c = clouds[0];

        if (1) {
          RandomPoints rpo;
          Block blockpoints;
          RandomParticles rpa(&blockpoints);
          rpa.n = 1000;
          rpa.pos_range = 500;
          rpa.scale_min = 1;
          rpa.scale_max = 30;

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
        WriteInBlocks say("front");
        say.block_pitch.set(1, 1, 1);
        say.point_genesis = &b;
        say.generate(c);
         
        c.pos.set(0, 0, -40);

        Cloud &c2 = clouds[1];
        c2.pos.set(0, 0, -45);
        WriteInBlocks say2("back");
        say2.block_pitch.set(1, 1, 1);
        say2.point_genesis = &b;
        say2.generate(c2);

        /*
        clouds.resize(3);
        clouds[2].pos.set(0, 0, 0);
        clouds[2].scale.set(1,  1, 1);
        say.generate(clouds[2]);
        */
      }

      change_color.pal = &palette;
      cloud_bank.add(change_color);
      cloud_bank.add(revolve);

      change_color2.pal = &palette;
      fixed_bank.add(change_color2);
    }

    void load() {
      metal.load("images/metal091.jpg");
    }

    void step() {

      fly.step();
      Pt dir = fly.mass.v.unit();
      if (! dir.zero())
        dir = (dir + fly.nose.unit()) / 2;
      else
        dir = fly.nose;

      cam.look_at(fly.pos + fly.top_lag,
                  dir.unit() * (-3),
                  fly.top_lag);

      for (int i = 0; i < clouds.size(); i++) {
        clouds[i].step();
      }
    }


    void draw() {
      for (int i = 0; i < clouds.size(); i++) {
        clouds[i].draw(as_quads, cloud_bank);
      }
      fly.draw(fixed_bank);
      fixed_bank.end_of_frame();
      cloud_bank.end_of_frame();
    }
};

World world;


void draw_scene()
{
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity( );

  world.cam.gluLookAt();
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
  default:
    printf("%d %f\n", axis, (float)axis_val);
    break;
  case 0:
    world.fly.roll_y = -(axis_val*axis_val*axis_val) / 10;
    break;
  case 4:
    world.fly.roll_x = (axis_val*axis_val*axis_val) / 10;
    break;
  case 3:
    world.fly.roll_z = (axis_val*axis_val*axis_val) / 10;
    break;
  case 5:
    // analog trigger ... -1 == not pressed, 0 = half, 1 = full
    // accel
#if 0
    world.fly.velocity.change = ((axis_val + 1) / 2) / 200;
#else
    world.fly.propulsion_forward = (axis_val + 1) / 2;
#endif
    break;
  case 2:
    // analog trigger ... -1 == not pressed, 0 = half, 1 = full
    // break
#if 0
    world.fly.velocity.change = -((axis_val + 1) / 2) / 200;
#else
    world.fly.propulsion_break = (axis_val + 1) / 2;
#endif
    break;
  }
}

void on_joy_button(ControllerState &ctrl, int button, bool down) {
  if (! down)
    return;

  switch (button) {
  case 0:
    if (world.fly.wings < 1.e-5) {
      world.fly.wings = .1;
      world.fly.do_cruise = true;
    }
    else
    if (world.fly.wings < .9) {
      world.fly.wings = 1;
      world.fly.do_cruise = true;
    }
    else {
      world.fly.wings = 0;
      world.fly.do_cruise = false;
    }
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
  gluPerspective(80,(double)W/H,.5,500);

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

    draw_scene();
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

