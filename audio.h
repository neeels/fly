#include <list>

#include "foreach.h"

class Sound {
  public:
    virtual ~Sound() {}

    /* return true when done, and this sound will be dropped. */
    virtual bool add_next(double *dest) = 0;
};

int Audio_render_thread(void *arg);
void Audio_play_callback(void *userdata, Uint8 *stream, int len);

struct _Audio {
  volatile bool is_on = false;

  const int channels = 2;
  const int frame_size = channels * sizeof(Sint16);
  const int sample_size = sizeof(Sint16);
  int buf_len = 1024;
  int buf_samples = buf_len / sample_size;
  int buf_frames = buf_len / frame_size;
  int freq = 44100;
  double *buf_play = 0;
  double *buf_render = 0;
  double *buf_zeros = 0;
  SDL_Thread *render_thread_token = 0;

  SDL_sem *please_render;

  SDL_AudioSpec audio_spec;

  list<Sound*> playing;

  void play(Sound *sound)
  {
    printf("Adding sound\n");
    playing.push_front(sound);
  }

  bool start()
  {
    please_render = SDL_CreateSemaphore(0);

    SDL_AudioSpec audio_want_spec;
    audio_want_spec.freq = freq;
    audio_want_spec.format = AUDIO_S16;
    audio_want_spec.channels = channels;
    audio_want_spec.samples = buf_frames;
    audio_want_spec.callback = Audio_play_callback;
    audio_want_spec.userdata = static_cast<void*>(this);

    if (SDL_OpenAudio(&audio_want_spec, &audio_spec) < 0) {
      fprintf(stderr, "Cannot open audio: %s\n", SDL_GetError());
      SDL_DestroySemaphore(please_render);
      return false;
    }

    printf("Opened audio:\n"
           "  %d ch %d Hz %d buf",
           audio_spec.channels,
           audio_spec.freq,
           audio_spec.samples);
    if (audio_spec.format == AUDIO_S16) {
      printf(" 16bit WAV");
    }
    printf("\n");

    if ((audio_spec.channels != channels)
        || (audio_spec.format != audio_want_spec.format))
    {
      printf("Audio: can't work with this format, sorry.\n");
      stop();
      return false;
    }

    is_on = true;

    buf_frames = audio_spec.samples;
    buf_samples = buf_frames * channels;
    buf_len = buf_samples * sample_size;

    buf_zeros = new double[buf_samples];
    for (int i = 0; i < buf_samples; i++)
      buf_zeros[i] = 0.;

    buf_play = new double[buf_samples];
    buf_render = new double[buf_samples];
    memcpy(buf_play, buf_zeros, buf_len);
    memcpy(buf_render, buf_zeros, buf_len);

    freq = audio_spec.freq;

    next();
    SDL_Thread *render_thread_token = SDL_CreateThread(Audio_render_thread, NULL);
    SDL_PauseAudio(0);
    return true;
  }

  void stop()
  {
    is_on = false;
    if (render_thread_token) {
      SDL_WaitThread(render_thread_token, NULL);
      render_thread_token = NULL;
    }
    SDL_CloseAudio();
    SDL_DestroySemaphore(please_render);
    delete[] buf_play;
    delete[] buf_render;
    delete[] buf_zeros;
    buf_play = 0;
    buf_render = 0;
    buf_zeros = 0;
  }

  void next()
  {
    double m;
    double *dest = buf_render;
    for (int i = 0; i < buf_samples; i++)
      dest[i] = 0.;

    m = 0;
    for (int i = 0; i < buf_samples; i++)
      m = max(m, dest[i]);

    int n = 0;
    foreach(p, playing) {
      n ++;
      if ((*p)->add_next(dest)) {
        Sound *s = *p;
        p = playing.erase(p);
        delete s;
        printf("drop\n");
      }
    }
    printf("%d sounds\n", n);
  }

  void render_thread()
  {
    double *tmp;
    while(is_on) {
      tmp = buf_play;
      buf_play = buf_render;
      buf_render = tmp;

      next();
      SDL_SemWait(please_render);
    }
  }

};

_Audio Audio;

void Audio_play_callback(void *userdata, Uint8 *stream, int len)
{
  double *p = (double*)Audio.buf_play;
  int n = len / Audio.sample_size;
  Sint16 *s = (Sint16*)stream;
  for (int i = 0; i < n; i++) {
    s[i] = p[i] * 32786.;
  }
  SDL_SemPost(Audio.please_render);
}

int Audio_render_thread(void *arg)
{
  printf("Audio render thread launched.\n");
  Audio.render_thread();
  printf("Audio render thread concluded.\n");
  return 0;
}

class Sine : public Sound {
  public:

    double amp;
    double freq_rad;
    double phase;

    Sine(double freq=440, double amp=.5, double phase = 0)
    {
      freq_rad = 2. * M_PI * freq;
      this->amp = amp;
      this->phase = phase;
    }

    virtual bool add_next(double *dest)
    {
      double f = (double)Audio.freq * Audio.channels;
      double t;
      double v;
      for (int i = 0; i < Audio.buf_samples; i+=2) {
        t = ((double)i) / f;
        v = amp * sin(phase + t * freq_rad);
        dest[i] = v;
        dest[i+1] = v;
      }
      t = ((double)Audio.buf_samples)/f;
      phase += t * freq_rad;
      phase -= trunc(phase / (2.*M_PI)) * (2.*M_PI);
      return false;
    }
};

class Mix : public Sound {
  public:

    list<Sound*> in;

    Mix()
    {}

    Mix *add(Sound *in)
    {
      this->in.push_front(in);
      return this;
    }

    virtual ~Mix()
    {
      foreach(s, in) {
        delete (*s);
        s = in.erase(s);
      }
    }

    virtual bool add_next(double *dest)
    {
      bool done = true;

      foreach(s, in) {
        done = (*s)->add_next(dest) && done;
      }

      return done;
    }

};

class Envelope : public Sound {
  public:

    Sound *in;

    int t;
    int a, d;

    Envelope(double a, double d, Sound *in)
    {
      this->in = in;
      this->a = a * Audio.freq;
      this->d = d * Audio.freq;
      t = 0;
    }

    virtual ~Envelope()
    {
      delete in;
      in = 0;
    }

    virtual bool add_next(double *dest)
    {
      double buf[Audio.buf_samples] = { 0. };
      bool done = in->add_next(buf);

      int i = 0;
      double v;
      for (;(t < a) && (i < Audio.buf_samples); i+=2, t++) {
        v = ((double)t)/a;
        v *= v;
        dest[i] = buf[i] * v;
        dest[i+1] = buf[i+1] * v;
      }
      for (;(t < d) && (i < Audio.buf_samples); i+=2, t++) {
        v = 1. - ((double)(t-a))/(double)d;
        v *= v;
        dest[i] = buf[i] * v;
        dest[i+1] = buf[i+1] * v;
      }
      return done || (t > (a + d));
    }

};
