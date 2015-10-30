#include <list>

#include "foreach.h"

const int Audio_buf_len = 1024;
const int Audio_channels = 2;
const int Audio_frame_size = Audio_channels * sizeof(Sint16);
int Audio_freq = 44100;
double Audio_buf[Audio_buf_len];

class Sound {
  public:

    virtual bool add_next() = 0;
};

void audio_play_callback(void *userdata, Uint8 *stream, int len);

class Audio {
  public:
    SDL_AudioSpec audio_spec;

    list<Sound> playing;

    bool start()
    {
      SDL_AudioSpec audio_want_spec;
      audio_want_spec.freq = Audio_freq;
      audio_want_spec.format = AUDIO_S16;
      audio_want_spec.channels = Audio_channels;
      audio_want_spec.samples = Audio_buf_len;
      audio_want_spec.callback = audio_play_callback;
      audio_want_spec.userdata = static_cast<void*>(this);

      if (SDL_OpenAudio(&audio_want_spec, &audio_spec) < 0) {
        fprintf(stderr, "Cannot open audio: %s\n", SDL_GetError());
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

      if ((audio_spec.samples != Audio_buf_len)
          || (audio_spec.channels != Audio_channels)
          || (audio_spec.format != audio_want_spec.format))
      {
        printf("Audio: can't work with this format, sorry.\n");
        stop();
        return false;
      }
      Audio_freq = audio_spec.freq;

      SDL_PauseAudio(0);
      return true;
    }

    void stop()
    {
      SDL_CloseAudio();
    }

    void next(Uint8 *stream, int len)
    {
      static const double zeros[Audio_buf_len] = {0.};
      memcpy(Audio_buf, zeros, sizeof(Audio_buf));

      foreach(p, playing) {
      }
    }
};

void audio_play_callback(void *userdata, Uint8 *stream, int len) {
  Audio &audio = *static_cast<Audio*>(userdata);
  audio.next(stream, len);
}


class Sine {
  public:

    double amp;
    double freq_rad;
    double phase;

    Sine(double freq=440, double amp=.5)
    {
      freq_rad = 2. * M_PI * freq;
      this->amp = amp;
    }

    virtual bool add_next()
    {
      double f = Audio_freq << 1;
      double t;
      double v;
      const int samples = Audio_buf_len >> 1;
      for (int i = 0; i < samples; i+=2) {
        t = ((double)i) / f;
        v = amp * sin(phase + Audio_freq * t);
        Audio_buf[i] = v;
        Audio_buf[i+1] = v;
      }
    }
};

