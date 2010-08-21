#include <cairo.h>
#include <math.h>

#include "SDL.h"

#include "path.h"

#define SDL_VIDEO_FLAGS 0
#define FRAMERATE 15
#define TOLERATE_N_MISSED_FRAMES 1

#define true -1
#define false 0

struct part {
  struct path paths;
};

void render_1x1_box(cairo_t *cr) {
  cairo_move_to(cr, -3/4.0, -1/2.0);
  cairo_line_to(cr,  3/4.0, -1/2.0);
  cairo_move_to(cr, -3/4.0,  1/2.0);
  cairo_line_to(cr,  3/4.0,  1/2.0);
  cairo_move_to(cr, -1/2.0, -3/4.0);
  cairo_line_to(cr, -1/2.0,  3/4.0);
  cairo_move_to(cr,  1/2.0, -3/4.0);
  cairo_line_to(cr,  1/2.0,  3/4.0);
  cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
  cairo_set_line_width(cr, 1/128.0);
  cairo_stroke(cr);
}

void path_render(
    cairo_t *cr,
    struct  path *path
  ) {
  struct path_segment *segment = (struct path_segment *) &path->segments;

  while (segment) {
    path_segment_render(cr, segment);
    segment = segment->next;
  }

  if (path->color) {
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_stroke_preserve(cr);
    cairo_set_source_rgb(cr, 1, 1, 1);
    if (path->fill) { cairo_fill(cr); }
  } else {
    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_stroke_preserve(cr);
    cairo_set_source_rgb(cr, 0, 0, 0);
    if (path->fill) { cairo_fill(cr); }
  }
}

void render_part(
    cairo_t *cr,
    struct part *part
  ) {
  struct path *path = &part->paths;

  cairo_set_line_width(cr, 1/32.0);
  while (path) {
    path_render(cr, path);
    path = path->next;
  }
}

void render_dots(cairo_t *cr, int subdiv) {
  int x, y;

  int a = 1 << subdiv;

  cairo_new_path(cr);

  for (y = -a; y <= a; y = y + 1) {
    for (x = -a; x <= a; x = x + 1) {
      cairo_new_sub_path(cr);
      cairo_arc(cr, (double) x / a, (double) y / a, 1/256.0, 0, 2 * M_PI);
    }
  }

  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_set_line_width(cr, 1/64.0);
  cairo_stroke_preserve(cr);
  cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
  cairo_set_line_width(cr, 1/128.0);
  cairo_stroke(cr);
}

void render_pen(cairo_t *cr, struct pos *pen) {
  cairo_arc(cr, pen->x, pen->y, 1/32.0, 0, 2 * M_PI);
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_set_line_width(cr, 1/64.0);
  cairo_stroke_preserve(cr);
  cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
  cairo_set_line_width(cr, 1/128.0);
  cairo_stroke(cr);
}

void render_outline(cairo_t *cr) {
  cairo_move_to(cr, 0, 0);
  cairo_save(cr);
  cairo_scale(cr, 1, -1);
  cairo_show_text(cr, "Testing");
  cairo_restore(cr);
}

void path_last_xy(
    struct path_segment *segment,
    struct pos          *pen
  ) {
  struct path_segment_move  *segment_move;
  struct path_segment_line  *segment_line;
  struct path_segment_curve *segment_curve;

  while (true) {
    if (PATH_MOVE == segment->type) {
      segment_move = (struct path_segment_move *) segment;

      pen->x = segment_move->p.x;
      pen->y = segment_move->p.y;
      break;
    }
    if (PATH_LINE == segment->type) {
      segment_line = (struct path_segment_line *) segment;

      pen->x = segment_line->p.x;
      pen->y = segment_line->p.y;
      break;
    }
    if (PATH_CURVE == segment->type) {
      segment_curve = (struct path_segment_curve *) segment;

      pen->x = segment_curve->p3.x;
      pen->y = segment_curve->p3.y;
      break;
    }
    segment = segment->prev;
  }
}

// segment_line = path_segment_line_append(segment_line, *pen);
struct path_segment_line *path_segment_line_append(
    struct path_segment *segment,
    struct pos          **pen
  ) {
  struct path_segment_line *next = (struct path_segment_line *)
    malloc( sizeof(struct path_segment_line) );

  *pen = &next->p;

  path_last_xy(segment, *pen);

  next->base.prev = segment;
  next->base.next = NULL;
  next->base.type = PATH_LINE;
  segment->next   = (struct path_segment *) next;

  return next;
}

/*
case SDLK_T:
  /* toggle segment type * /
  break;
*/

int main(int argc, char **argv) {
  int    hres, vres, width, height, pixels;
  double aspect, scale;

  SDL_Surface    *sdl_surface;
  Uint32         next_frame;
  cairo_t        *cr;
  cairo_matrix_t cm_display, cm_field, cm_panel;

  int missed_frames = 0;

  int running = true;
  int subdiv  = 1;

  struct part part = {
    NULL, NULL, NULL, NULL, PATH_MOVE, 0, 0, 0, 0
  };

  struct path_segment *segment = (struct path_segment *)
    &part.paths.segments.base; /* first segment of the first path */
  struct path         *path    = &part.paths;

  struct pos *pen = &part.paths.segments.p;

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

    double scale = 1.0 / (2.0) / aspect;

    cairo_matrix_init_identity(m);
    cairo_matrix_multiply(m, m, &cm_display);

    cairo_matrix_scale(m, scale, scale);
  }

  { /* delay */
    Uint32 now = SDL_GetTicks();

    next_frame = now + 1000.0 / FRAMERATE;
  }

  fprintf(stderr, "-- MARK -- setup complete\n");

  SDL_LockSurface(sdl_surface);

  while (running) {

    { /* Render Frame */

      /* clear screen */
      cairo_set_operator(cr, CAIRO_OPERATOR_CLEAR);
      cairo_paint(cr);
      cairo_set_operator(cr, CAIRO_OPERATOR_OVER);

      /* Render Field */

      cairo_set_matrix(cr, &cm_field);

      render_1x1_box(cr);
      render_part(cr, &part);
      render_dots(cr, subdiv);
      render_pen(cr, pen);
    }

    { /* Update Display */
      SDL_UnlockSurface(sdl_surface);
      SDL_Flip(sdl_surface);
      SDL_LockSurface(sdl_surface);
    }

    { /* Delay */
      Uint32 now = SDL_GetTicks();

      if (now < next_frame) {
        missed_frames = 0;
        SDL_Delay(next_frame - now);
        next_frame = next_frame + 1000.0 / FRAMERATE;
      } else {
        missed_frames = missed_frames + 1;
        if (1 + TOLERATE_N_MISSED_FRAMES == missed_frames) {
          fprintf(stderr, "delay counter reset\n");
          next_frame = now + 1000.0 / FRAMERATE;
        } else {
          next_frame = next_frame + 1000.0 / FRAMERATE;
        }
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
            switch (event.key.keysym.sym) {
              case SDLK_q:
                running = 0;
                break;
              case SDLK_1:
                subdiv = 1;
                break;
              case SDLK_2:
                subdiv = 2;
                break;
              case SDLK_3:
                subdiv = 3;
                break;
              case SDLK_4:
                subdiv = 4;
                break;
              case SDLK_LEFT:
                pen->x = pen->x - 1.0 / (1 << subdiv);
                break;
              case SDLK_RIGHT:
                pen->x = pen->x + 1.0 / (1 << subdiv);
                break;
              case SDLK_DOWN:
                pen->y = pen->y - 1.0 / (1 << subdiv);
                break;
              case SDLK_UP:
                pen->y = pen->y + 1.0 / (1 << subdiv);
                break;
              case SDLK_RETURN:
                segment = (struct path_segment *)
                  path_segment_line_append(segment, &pen);
                fprintf(stderr, "-- MARK -- new segment created\n");
                break;
              case SDLK_c:
                path->color = 1 - path->color;
                break;
              case SDLK_f:
                path->fill  = 1 - path->fill;
                break;
            }
            break;
        }
      }
    }

  }

  SDL_UnlockSurface(sdl_surface);

  /* CLEANUP */

  cairo_destroy(cr);
  SDL_Quit();

  return 0;
}
