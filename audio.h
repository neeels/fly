#include <list>

#include "foreach.h"

class Sound {
    bool is_dropped = false;
  public:
    virtual ~Sound() {}

    void drop();

    inline bool done() const
    {
      return is_dropped;
    }

    /* return true when done, and this sound will be done. */
    virtual void add_next(double *dest) = 0;
    virtual void gc() {};
};

class Mix : public Sound {
  public:

    list<Sound*> in;

    Mix()
    {}

    Mix *add(Sound *new_in)
    {
      in.push_back(new_in);
      return this;
    }

    virtual ~Mix()
    {
      foreach(s, in) {
        delete (*s);
        s = in.erase(s);
      }
    }

    virtual void add_next(double *dest)
    {
      bool all_done = true;

      foreach(s, in) {
        (*s)->add_next(dest);
        all_done = all_done && (*s)->done();
      }

      if (all_done) {
        drop();
      }
    }

    virtual void gc()
    {
      foreach (i, in) {
        Sound *s = *i;
        if (s->done()) {
          i = in.erase(i);
          delete s;
        }
        else
          s->gc();
      }
    }
};


int Audio_render_thread(void *arg);
int Audio_gc_thread(void *arg);
void Audio_play_callback(void *userdata, Uint8 *stream, int len);

struct _Audio {
private:
  Mix playing;
public:
  volatile bool is_on;

  const int channels = 2;
  const int frame_size = channels * sizeof(Sint16);
  const int sample_size = sizeof(Sint16);
  int buf_len;
  int buf_samples;
  int buf_frames;
  int freq;
  double *buf_play;
  double *buf_render;
  double *buf_zeros;
  SDL_Thread *render_thread_token = 0;

  SDL_sem *please_render;

  SDL_AudioSpec audio_spec;

  list<Sound*> start_playing_queue;
  SDL_sem *change_playing_mutex;

  const char *write_to_path;
  FILE *write_to_f;
  
  _Audio() :
     is_on(false),
     buf_len(4096),
     freq(44100),
     write_to_path(NULL),
     write_to_f(NULL),
     buf_play(0),
     buf_render(0),
     buf_zeros(0),
     render_thread_token(0)
  {
    buf_samples = buf_len / sample_size;
    buf_frames = buf_len / frame_size;
    please_render = SDL_CreateSemaphore(0);
    change_playing_mutex = SDL_CreateSemaphore(1);
  }

  bool start()
  {
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

    // if a write_to_path is set, this opens the file.
    write_to(write_to_path);

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
    render_thread_token = SDL_CreateThread(Audio_render_thread, "audio", NULL);
    //gc_thread_token = SDL_CreateThread(Audio_gc_thread, NULL);
    SDL_PauseAudio(0);
    return true;
  }

  void play(Sound *sound)
  {
    start_playing_queue.push_back(sound);
  }

  void stop()
  {
    is_on = false;
    SDL_SemPost(please_render);
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
    if (write_to_f) {
      fclose(write_to_f);
      write_to_f = NULL;
    }
  }

  void next();

  void write_to(const char *path)
  {
    write_to_path = path;
    if (is_on && (write_to_f == NULL))
      write_to_f = fopen(write_to_path, "w");
  }

  void render_thread()
  {
    double *tmp;
#define MEASURE 0
#if MEASURE
    unsigned int t0, t1, rendered, waited;
#endif

    while(is_on) {

#if MEASURE
      t0 = SDL_GetTicks();
#endif

      tmp = buf_play;
      buf_play = buf_render;
      buf_render = tmp;

      next();

      if (write_to_f) {
        fwrite(buf_render, sizeof(double), buf_samples, write_to_f);
      }

      gc();

#if MEASURE
      t1 = SDL_GetTicks();
#endif

      SDL_SemWait(please_render);

#if MEASURE
      waited = SDL_GetTicks() - t1;
      rendered = t1 - t0;
      printf("audio rendered %d ms, waited %d ms\n", rendered, waited);
#endif
    }
  }

  void gc();

};

_Audio Audio;

void Audio_play_callback(void *userdata, Uint8 *stream, int len)
{
  if (len != Audio.buf_len)
    printf("AUDIO BUF LEN MISMATCH! %d %d\n", Audio.buf_len, len);
  double *p = (double*)Audio.buf_play;
  int n = len / Audio.sample_size;
  Sint16 *s = (Sint16*)stream;
  for (int i = 0; i < n; i++) {
    s[i] = min((Sint16)0x7fff, max((Sint16)0x8000, (Sint16)(p[i] * ((double)0x7fff))));
  }
  
  int semval = SDL_SemValue(Audio.please_render);
  if (semval)
    printf("underrun %d\n", semval);
  SDL_SemPost(Audio.please_render);
}

int Audio_render_thread(void *arg)
{
  printf("Audio render thread launched.\n");
  Audio.render_thread();
  printf("Audio render thread concluded.\n");
  return 0;
}

void Sound::drop()
{
  is_dropped = true;
}

void _Audio::next()
{
  double *dest = buf_render;
  for (int i = 0; i < buf_samples; i++)
    dest[i] = 0.;

  SDL_SemWait(change_playing_mutex);
  playing.add_next(dest);
  SDL_SemPost(change_playing_mutex);

  if (start_playing_queue.size()) {
    foreach (i, start_playing_queue) {
      Sound *s = *i;
      SDL_SemWait(change_playing_mutex);
      playing.add(s);
      SDL_SemPost(change_playing_mutex);
      i = start_playing_queue.erase(i);
    }
  }
}

void _Audio::gc()
{
  SDL_SemWait(change_playing_mutex);
  playing.gc();
  SDL_SemPost(change_playing_mutex);
}

double audible(double freq, double maxf=880, double minf=50)
{
  if (freq < 1e-6)
    return minf;
  if (freq < minf)
    return freq * ceil(minf/freq);
  if (freq > maxf)
    return freq / ceil(freq/maxf);
  return freq;
}

double octave(double freq, int octaves=1, double base=220)
{
  return audible(freq, base * (octaves + 1), base);
}


class Sine : public Sound {
  public:

    double amp;
    double freq_rad;
    double phase;

    Sine(double freq=440, double amp=.1, double phase = 0)
    {
      freq_rad = 2. * M_PI * freq;
      this->amp = amp;
      this->phase = phase;
    }

    virtual ~Sine() {};

    virtual void add_next(double *dest)
    {
      double f = (double)Audio.freq * Audio.channels;
      double t;
      double v;
      for (int i = 0; i < Audio.buf_samples; i+=2) {
        t = ((double)i) / f;
        v = amp * sin(phase + t * freq_rad);
        dest[i] += v;
        dest[i+1] += v;
      }
      t = ((double)Audio.buf_samples)/f;
      phase += t * freq_rad;
      phase -= trunc(phase / (2.*M_PI)) * (2.*M_PI);
    }
};

class Envelope : public Sound {
  public:

    Sound *in;

    int t;
    int a, h, d, ah, ahd;

    Envelope(double a_, double h_, double d_, Sound *in_)
    {
      in = in_;
      a = a_ * Audio.freq;
      h = h_ * Audio.freq;
      d = d_ * Audio.freq;
      ah = a + h;
      ahd = ah + d;
      t = 0;
    }

    virtual ~Envelope()
    {
      delete in;
      in = 0;
    }

    virtual void add_next(double *dest)
    {
      double buf[Audio.buf_samples] = { 0. };
      in->add_next(buf);

      int i = 0;
      double v;
      for (;(t < a) && (i < Audio.buf_samples); i+=2, t++) {
        v = ((double)t)/a;
        v *= v;
        dest[i] += buf[i] * v;
        dest[i+1] += buf[i+1] * v;
      }
      for (;(t < ah) && (i < Audio.buf_samples); i+=2, t++) {
        dest[i] += buf[i];
        dest[i+1] += buf[i+1];
      }
      for (;(t < ahd) && (i < Audio.buf_samples); i+=2, t++) {
        v = 1. - ((double)(t-ah))/(double)d;
        v *= v;
        dest[i] += buf[i] * v;
        dest[i+1] += buf[i+1] * v;
      }
      if (in->done() || (t >= ahd)) {
        drop();
      }
    }

};
