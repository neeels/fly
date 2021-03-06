#include <SDL2/SDL.h>
#include <vector>
using namespace std;

#include "foreach.h"

class ControllerState {
  public:
    void *data;

    int selected_layer;
    bool clamp_button_held;

    ControllerState() {
      selected_layer = 0;
      clamp_button_held = 0;
    }
};

/* You need to define these functions: */
void on_joy_axis(ControllerState &ctrl, int axis, double axis_val);
void on_joy_button(ControllerState &ctrl, int button, bool down);


vector<ControllerState> controllers;

void controller_state_init(int n_joysticks, void *data=NULL) {
  controllers.resize(n_joysticks);
  foreach(c, controllers) {
    c->data = data;
  }
}

void handle_joystick_events(SDL_Event &event) {

  switch (event.type) 
  {
    case SDL_JOYHATMOTION:  /* Handle Hat Motion */
      if (event.jhat.which < controllers.size())
      {
        ControllerState &ctrl = controllers[event.jaxis.which];

        int l = ctrl.selected_layer;
        switch(event.jhat.value)
        {
          case 1:
            l = 0;
            break;
          case 2:
            l = 1;
            break;
          case 4:
            l = 2;
            break;
          case 8:
            l = 3;
            break;

          default:
            break;
        }

        if (l != ctrl.selected_layer) {
          // layer changed
          ctrl.selected_layer = l;
          printf("LAYER %d\n", ctrl.selected_layer);
        }
      }
      break;

    case SDL_JOYAXISMOTION:
      if (event.jaxis.which < controllers.size())
      {
        double axis_val = event.jaxis.value;
        axis_val /= 32768;

        ControllerState &ctrl = controllers[event.jaxis.which];

        int l = ctrl.selected_layer;
        switch(event.jaxis.axis)
        {
          case 6:
            if (axis_val < -.2) {
              l = 1;
            }
            else
            if (axis_val > .2) {
              l = 3;
            }
            break;

          case 7:
            if (axis_val < -.2) {
              l = 2;
            }
            else
            if (axis_val > .2) {
              l = 0;
            }
            break;

          default:
            break;
        }

        if (l != ctrl.selected_layer) {
          // layer changed
          ctrl.selected_layer = l;
          printf("LAYER %d %f\n", ctrl.selected_layer, axis_val);
        }

        on_joy_axis(ctrl,
                    event.jaxis.axis,
                    axis_val);
      }
      break;


    case SDL_JOYBUTTONDOWN:
    case SDL_JOYBUTTONUP:
      if (event.jbutton.which < controllers.size())
      {
        ControllerState &ctrl = controllers[event.jbutton.which];
        int l = ctrl.selected_layer;

        on_joy_button(ctrl, event.jbutton.button,
                      event.type == SDL_JOYBUTTONDOWN);
      }
      break;

    #if 0
    case SDL_JOYHATMOTION:  /* Handle Hat Motion */
      printf("%2d: hat %d = %d\n",
             event.jhat.which, event.jhat.hat, event.jhat.value);
      break;
    #endif

    case SDL_JOYBALLMOTION:  /* Handle Joyball Motion */
      printf("%2d: ball %d += %d, %d\n",
             event.jball.which, event.jball.ball,
             event.jball.xrel, event.jball.yrel);
      break;
  }

}

