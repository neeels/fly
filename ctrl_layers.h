
class ControllerState {
  public:
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

void controller_state_init(int n_joysticks) {
  ControllerState ctrlr_state;

  while (controllers.size() < n_joysticks)
    controllers.push_back(ctrlr_state);
}

void handle_joystick_events(SDL_Event &event) {

  switch (event.type) 
  {
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

    case SDL_JOYHATMOTION:  /* Handle Hat Motion */
      printf("%2d: hat %d = %d\n",
             event.jhat.which, event.jhat.hat, event.jhat.value);
      break;

    case SDL_JOYBALLMOTION:  /* Handle Joyball Motion */
      printf("%2d: ball %d += %d, %d\n",
             event.jball.which, event.jball.ball,
             event.jball.xrel, event.jball.yrel);
      break;
  }

}

