#include <cairo.h>
#include <math.h>

#include "SDL.h"

#define SDL_VIDEO_FLAGS 0
#define FRAMERATE 60

#define WHOMP_SPEED 3/2.0

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
};

struct player {
  double x, y;
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
  struct player player;
  struct whomp  whomp;

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
  whomp.x = 0;
  whomp.y = 0;

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
      Uint8 *keystate = SDL_GetKeyState(NULL);

      player.x = 0;
      player.y = 0;

      player.x = player.x - keystate[SDLK_LEFT];
      player.x = player.x + keystate[SDLK_RIGHT];
      player.y = player.y - keystate[SDLK_DOWN];
      player.y = player.y + keystate[SDLK_UP];

      if (-1 > player.x) { player.x = -1; }
      if ( 1 < player.x) { player.x =  1; }
      if (-1 > player.y) { player.y = -1; }
      if ( 1 < player.y) { player.y =  1; }
    }

    { /* Animate */
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
    }

  }

  SDL_UnlockSurface(sdl_surface);

  /* CLEANUP */

  cairo_destroy(cr);
  SDL_Quit();

  return 0;
}
