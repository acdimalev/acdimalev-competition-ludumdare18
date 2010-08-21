#include <cairo.h>
#include <math.h>

#include "SDL.h"

#define SDL_VIDEO_FLAGS 0
#define FRAMERATE 60

#define WHOMP_SPEED    3/2.0
#define SQUIBBLE_SPEED 2/1.0

#define FIELD_RES (1 << 4)
#define STUB_PATTERN { \
  0,0,0,0,0,0,0,0, \
  0,0,0,0,0,0,0,0, \
  0,0,0,0,0,0,0,0, \
  0,0,0,0,0,0,0,0, \
  0,0,0,0,1,1,0,0, \
  0,0,0,0,1,0,0,0, \
  0,0,0,0,0,0,0,0, \
  0,0,0,0,0,0,0,0  \
}

#define FIELD_RES_IS_INVALID (FIELD_RES % 2)

#define true -1
#define false 0

struct whomp {
  double x, y;
  struct squibble *squibble_held;
};

enum squibble_action {
  SQUIBBLE_IS_CHARGING_INTO_THE_FRAY,
  SQUIBBLE_IS_BEING_HELD
};

struct squibble {
  double x, y, vx, vy;
  enum   squibble_action action;
};

struct player {
  double x, y;
  int    action;
};

void generate_field(int *field) {
  int x, y;

  int s = FIELD_RES / 2;
  int pattern[] = STUB_PATTERN;

  /* load stub pattern */

  for (y = 0; y < s; y = y + 1) {
    for (x = 0; x < s; x = x + 1) {
      field[2 * y * s + x] = pattern[y * s + x];
    }
  }

  /* copy pattern to other quadrants */

  for (y = 0; y < s; y = y + 1) {
    for (x = 0; x < s; x = x + 1) {
      int yi = 2 * s - y - 1;
      int xi = 2 * s - x - 1;

      field[2 * y  * s + xi] = field[2 * y * s + x];
      field[2 * yi * s + x ] = field[2 * y * s + x];
      field[2 * yi * s + xi] = field[2 * y * s + x];
    }
  }
}

int collides_with_field(int *field, double x, double y, double r) {
  int xi, yi;

  int y0 = floor(y - r);
  int y1 = ceil(y + r);
  int x0 = floor(x - r);
  int x1 = ceil(x + r);

  int s  = FIELD_RES / 2;

  for (yi = y0; yi < y1; yi = yi + 1) {
    for (xi = x0; xi < x1; xi = xi + 1) {
      if ( field[(2 * (yi + s) * s + (xi + s))] ) {
        return true;
      }
    }
  }

  return false;
}

int main(int argc, char **argv) {
  int    hres, vres, width, height, pixels;
  double aspect, scale;

  SDL_Surface    *sdl_surface;
  Uint32         next_frame;
  cairo_t        *cr;
  cairo_matrix_t cm_display, cm_field;

  int    field[FIELD_RES * FIELD_RES];
  struct player   player;
  struct whomp    whomp;
  struct squibble squibble;

  int running = true;

  if (FIELD_RES_IS_INVALID) {
    fprintf(stderr, "Field resolution must be a multiple of two.\n");
    exit(0);
  }

  /* SETUP */

  SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER);

  { /* acquire video information */
    /* TODO: correct pixel format */
    SDL_Rect **modes = SDL_ListModes(NULL, SDL_VIDEO_FLAGS);

    if (NULL == modes) {
      fprintf(stderr, "no video modes available\n");
      exit(0);
    }

    if ( (SDL_Rect **) -1 == modes ) {
      hres = 640;
      vres = 480;
    } else {
      hres = modes[0]->w;
      vres = modes[0]->h;
    }

    width  = hres;
    height = vres;
    aspect = (double) hres / vres;
    pixels = hres * vres;
    scale  = sqrt(pixels);

    /* TODO: cap resolution if possible */
  }

  sdl_surface = SDL_SetVideoMode(hres, vres, 32, SDL_VIDEO_FLAGS);

  { /* Cairo */
    cairo_surface_t *cr_surface;
    cr_surface = cairo_image_surface_create_for_data(
        sdl_surface->pixels,
        CAIRO_FORMAT_RGB24,
        sdl_surface->w,
        sdl_surface->h,
        sdl_surface->pitch );
    cr = cairo_create(cr_surface);
    cairo_surface_destroy(cr_surface);
  }

  { /* screen-space transformation */
    cairo_matrix_t *m = &cm_display;

    cairo_matrix_init_identity(m);

    /* Cartesian */
    cairo_matrix_translate(m, hres/2.0, vres/2.0);
    cairo_matrix_scale(m, 1, -1);

    /* fixed scale */
    cairo_matrix_scale(m,
        scale * hres / width,
        scale * vres / height );
  }

  { /* field transformation */
    cairo_matrix_t *m = &cm_field;

    double scale = 1.0 / FIELD_RES / aspect;

    cairo_matrix_init_identity(m);
    cairo_matrix_multiply(m, m, &cm_display);

    cairo_matrix_scale(m, scale, scale);
  }

  { /* delay */
    Uint32 now = SDL_GetTicks();

    next_frame = now + 1000.0 / FRAMERATE;
  }

  generate_field(field);

  player.action = 0;

  whomp.x = 0;
  whomp.y = 0;

  whomp.squibble_held = NULL;

  {
    double x = -FIELD_RES / 2.0;
    double y = 4;

    double d = sqrt(x * x + y * y);

    squibble.x = x;
    squibble.y = y;
    squibble.action = SQUIBBLE_IS_CHARGING_INTO_THE_FRAY;
    squibble.vx = -x / d;
    squibble.vy = -y / d;
  }

  fprintf(stderr, "-- MARK -- setup complete\n");

  /* RUNNING */

  SDL_LockSurface(sdl_surface);

  while (running) {

    { /* Render Frame */
      int x, y;

      int s = FIELD_RES / 2;

      /* clear screen */
      cairo_set_operator(cr, CAIRO_OPERATOR_CLEAR);
      cairo_paint(cr);
      cairo_set_operator(cr, CAIRO_OPERATOR_OVER);

      /* draw field */
      for (y = 0; y < 2 * s; y = y + 1) {
        for (x = 0; x < 2 * s; x = x + 1) {
          if (! field[2 * y * s + x]) { continue; }

          cairo_set_matrix(cr, &cm_field);
          cairo_translate(cr, x - s + 1/2.0, y - s + 1/2.0);

          cairo_move_to(cr, -1/2.0, -1/2.0);
          cairo_line_to(cr,  1/2.0, -1/2.0);
          cairo_line_to(cr,  1/2.0,  1/2.0);
          cairo_line_to(cr, -1/2.0,  1/2.0);
          cairo_close_path(cr);
        }
      }
      cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
      cairo_fill(cr);

      /* draw whomp */
      cairo_set_matrix(cr, &cm_field);
      cairo_translate(cr, whomp.x, whomp.y);
      cairo_move_to(cr, -1/2.0, -1/2.0);
      cairo_line_to(cr,  1/2.0, -1/2.0);
      cairo_line_to(cr,  1/2.0,  1/2.0);
      cairo_line_to(cr, -1/2.0,  1/2.0);
      cairo_close_path(cr);
      cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
      cairo_set_line_width(cr, 1/8.0);
      cairo_stroke(cr);

      /* draw squibble */
      cairo_set_matrix(cr, &cm_field);
      cairo_translate(cr, squibble.x, squibble.y);
      cairo_arc(cr, 0, 0, 1/2.0, 0, 2 * M_PI);
      cairo_close_path(cr);
      cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
      cairo_set_line_width(cr, 1/8.0);
      cairo_stroke(cr);
    }

    { /* Update Display */
      SDL_UnlockSurface(sdl_surface);
      SDL_Flip(sdl_surface);
      SDL_LockSurface(sdl_surface);
    }

    { /* Delay */
      Uint32 now = SDL_GetTicks();

      if (now < next_frame) {
        SDL_Delay(next_frame - now);
        next_frame = next_frame + 1000.0 / FRAMERATE;
      } else {
        fprintf(stderr, "delay counter reset\n");
        next_frame = now + 1000.0 / FRAMERATE;
      }
    }

    { /* Handle Events */
      SDL_Event event;

      while ( SDL_PollEvent(&event) ) {
        switch (event.type) {
          case SDL_QUIT:
            running = 0;
            break;
          case SDL_KEYDOWN:
            if (SDLK_q == event.key.keysym.sym) {
              running = 0;
            }
            break;
        }
      }
    }

    { /* Gather Input */
      double len;

      Uint8 *keystate = SDL_GetKeyState(NULL);

      double x = 0;
      double y = 0;

      int action = 0;

      x = x - keystate[SDLK_LEFT];
      x = x + keystate[SDLK_RIGHT];
      y = y - keystate[SDLK_DOWN];
      y = y + keystate[SDLK_UP];

      action = action + keystate[SDLK_SPACE];
      action = action + keystate[SDLK_LCTRL];

      len = sqrt(x * x + y * y);

      if (1 < len) {
        x = x / len;
        y = y / len;
      }

      if (action) { action = 1; }

      if (! action) {
        player.action = 0;
      }
      if (2 != player.action) {
        player.action = player.action + action;
      }

      player.x = x;
      player.y = y;

      player.action = action;
    }

    { /* Animate */

      { /* Whomp */
        double x0 = whomp.x;
        double y0 = whomp.y;
        double r  = 1/2.0;

        double x = x0 + player.x * WHOMP_SPEED / FRAMERATE;
        double y = y0 + player.y * WHOMP_SPEED / FRAMERATE;

        if ( collides_with_field(field, x, y0, r) ) {
          x = x0;
        }
        if ( collides_with_field(field, x0, y, r) ) {
          y = y0;
        }
        if ( collides_with_field(field, x, y, r) ) {
          x = x0;
          y = y0;
        }

        whomp.x = x;
        whomp.y = y;

        if (whomp.squibble_held) {
          struct squibble *squibble = whomp.squibble_held;

          squibble->x = whomp.x;
          squibble->y = whomp.y;
        } else {
          if (1 == player.action) {
            double x = squibble.x - whomp.x;
            double y = squibble.y - whomp.y;
            double d = sqrt(x * x + y * y);

            if (1 > d) {
              squibble.x  = whomp.x;
              squibble.y  = whomp.y;
              squibble.vx = 0;
              squibble.vy = 0;

              squibble.action = SQUIBBLE_IS_BEING_HELD;

              whomp.squibble_held = &squibble;
            }
          }
        }
      }

      { /* Squibble */
        double x0 = squibble.x;
        double y0 = squibble.y;

        double x = x0 + squibble.vx * SQUIBBLE_SPEED / FRAMERATE;
        double y = y0 + squibble.vy * SQUIBBLE_SPEED / FRAMERATE;
        double r = 1/2.0;

        switch (squibble.action) {
          case SQUIBBLE_IS_CHARGING_INTO_THE_FRAY:
            if ( abs(x) < FIELD_RES / 4.0 && abs(y) < FIELD_RES / 4.0 ) {
              squibble.vx = 0;
              squibble.vy = 0;
            }
            break;
        }

        if ( collides_with_field(field, x, y0, r) ) {
          x = x0;
        }
        if ( collides_with_field(field, x0, y, r) ) {
          y = y0;
        }
        if ( collides_with_field(field, x, y, r) ) {
          x = x0;
          y = y0;
        }

        squibble.x = x;
        squibble.y = y;
      }

    }

  }

  SDL_UnlockSurface(sdl_surface);

  /* CLEANUP */

  cairo_destroy(cr);
  SDL_Quit();

  return 0;
}
