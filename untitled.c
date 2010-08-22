#include <cairo.h>
#include <math.h>

#include "SDL.h"

#undef DEBUG_SQUIBBLE_CAN_SEE

#define SDL_VIDEO_FLAGS 0
/* #define SDL_VIDEO_FLAGS SDL_FULLSCREEN */
#define FRAMERATE 30

#define WHOMP_MAX 1
#define SQUIBBLE_MAX (1 << 4)

#define FIELD_RES (1 << 4)

#define WHOMP_SPEED (3/2.0)
#define SQUIBBLE_CHARGE_SPEED (2/1.0)
#define SQUIBBLE_WANDER_SPEED (1/1.0)
#define SQUIBBLE_BOREDOM_TIMER_DEFAULT (2/1.0)
#define SQUIBBLE_VIEW_DISTANCE (1/2.0 * FIELD_RES)
#define SQUIBBLE_VIEW_ANGLE (1/4.0 * M_PI)
#define SQUIBBLE_NEAR_PROXIMITY (3/2.0)
#define WHOMP_RADIUS (1/2.0)
#define SQUIBBLE_RADIUS (1/3.0)

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

int FIELD[FIELD_RES][FIELD_RES];
int OVERLAY[FIELD_RES][FIELD_RES];

int field(int x, int y) {
  int s = FIELD_RES / 2;

  if (x < -s) { return 0; }
  if (x >= s) { return 0; }
  if (y < -s) { return 0; }
  if (y >= s) { return 0; }

  return FIELD[y + s][x + s];
}

int field_mark(int x, int y) {
  int s = FIELD_RES / 2;

  if (x < -s) { return 0; }
  if (x >= s) { return 0; }
  if (y < -s) { return 0; }
  if (y >= s) { return 0; }

  OVERLAY[y + s][x + s] = 1;

  return FIELD[y + s][x + s];
}

void clear_overlay(void) {
  int x, y;

  for (y = 0; y < FIELD_RES; y = y + 1) {
    for (x = 0; x < FIELD_RES; x = x + 1) {
      OVERLAY[y][x] = 0;
    }
  }
}

struct pos {
  double x, y;
};

struct pos POS_CENTER = {0, 0};

double random_angle(void) {
  return 2 * M_PI * rand() / (RAND_MAX + 1.0);
}

double angle_between(struct pos *a, struct pos *b) {
  double x = b->x - a->x;
  double y = b->y - a->y;

  if (x) {
    if (x > 0) {
      return atan(y/x);
    } else {
      return atan(y/x) + M_PI;
    }
  } else {
    if (y > 0) {
      return 1/2.0 * M_PI;
    } else {
      return 3/2.0 * M_PI;
    }
  }
}

enum whomp_state {
  WHOMP_IS_DEAD,
  WHOMP_IS_ALIVE
};

struct whomp {
  struct pos p;
  double a;
  enum   whomp_state state;
  struct squibble    *squibble_held;
} WHOMPS[WHOMP_MAX];

enum squibble_state {
  SQUIBBLE_IS_DEAD,
  SQUIBBLE_IS_CHARGING,
  SQUIBBLE_IS_WANDERING,
  SQUIBBLE_IS_LOOKING_AROUND,
  SQUIBBLE_IS_SHOCKED,
  SQUIBBLE_IS_BEING_HELD,
  SQUIBBLE_IS_AIRBORNE
};

struct squibble {
  struct pos p;
  double a;
  enum   squibble_state state;
  double boredom_timer;
  double explosion_timer;
  struct whomp    *chasing_whomp;
  struct squibble *chasing_squibble;
} SQUIBBLES[SQUIBBLE_MAX];

void squibble_new(void) {
  struct squibble *squibble;
  int    i;

  for (i = 0; i < SQUIBBLE_MAX; i = i + 1) {
    squibble = &SQUIBBLES[i];
    if (SQUIBBLE_IS_DEAD == squibble->state) { break; }
  }

  if (SQUIBBLE_MAX == i) {
    fprintf(stderr, "out of squibbles\n");
    return;
  }

  {
    double d;

    double a = random_angle();

    double x = cos(a);
    double y = sin(a);

    if ( abs(x) > abs(y) ) {
      d = FIELD_RES / 2.0 / fabs(x);
    } else {
      d = FIELD_RES / 2.0 / fabs(y);
    }

    squibble->p.x = d * x;
    squibble->p.y = d * y;
    squibble->a = angle_between(&squibble->p, &POS_CENTER);
    squibble->state = SQUIBBLE_IS_CHARGING;
    squibble->boredom_timer = SQUIBBLE_BOREDOM_TIMER_DEFAULT;
  }
}

void squibble_move(struct squibble *squibble, double speed) {
  double x0 = squibble->p.x;
  double y0 = squibble->p.y;
  double vx = cos(squibble->a) * speed;
  double vy = sin(squibble->a) * speed;

  double x = x0 + vx / FRAMERATE;
  double y = y0 + vy / FRAMERATE;
  double r = SQUIBBLE_RADIUS;

  if ( collides_with_field(x, y0, r) ) {
    x = x0;
  }
  if ( collides_with_field(x0, y, r) ) {
    y = y0;
  }
  if ( collides_with_field(x, y, r) ) {
    x = x0;
    y = y0;
  }

  squibble->p.x = x;
  squibble->p.y = y;
}

int squibble_can_see(struct squibble *squibble, struct pos *pos) {
  int    u, ud, uf;
  double v, vd;

  double x1 = squibble->p.x;
  double y1 = squibble->p.y;
  double x2 = pos->x;
  double y2 = pos->y;
  double a  = angle_between(&squibble->p, pos) - squibble->a;

  double xd = x2 - x1;
  double yd = y2 - y1;
  int    b  = floor(a / (2*M_PI));

  double d = sqrt(xd * xd + yd * yd);

  a = a - b * 2*M_PI;
  if (a > M_PI) { a = 2*M_PI - a; }

  if (SQUIBBLE_VIEW_DISTANCE < d) { return false; }
  if (SQUIBBLE_VIEW_ANGLE    < a) { return false; }

#ifdef DEBUG_SQUIBBLE_CAN_SEE
  field_mark(-8, -8);
  field_mark( 7, -8);
  field_mark( 7,  7);
  field_mark(-8,  7);
#endif

  if ( fabs(xd) > fabs(yd) ) {
    ud = xd < 0 ? -1 : 1;
    uf = x2 + ud;
    vd = yd / fabs(xd);
    for (u = floor(x1), v = y1; u != uf; u = u + ud, v = v + vd) {
#ifdef DEBUG_SQUIBBLE_CAN_SEE
      if ( field_mark(u, floor(v)) ) { return false; }
#else
      if ( field(u, floor(v)) ) { return false; }
#endif
    }
  } else {
    ud = yd < 0 ? -1 : 1;
    uf = y2 + ud;
    vd = xd / fabs(yd);
    for (u = floor(y1), v = x1; u != uf; u = u + ud, v = v + vd) {
#ifdef DEBUG_SQUIBBLE_CAN_SEE
      if ( field_mark(v, floor(u)) ) { return false; }
#else
      if ( field(u, floor(v)) ) { return false; }
#endif
    }
  }

  return true;
}

int squibble_is_near(struct squibble *squibble, struct pos *pos) {
  double x = pos->x - squibble->p.x;
  double y = pos->y - squibble->p.y;

  double d = sqrt(x * x + y * y);

  if (SQUIBBLE_NEAR_PROXIMITY > d) { return true; }

  return false;
}

void squibble_do_nothing(struct squibble *squibble) {
}

void squibble_charge(struct squibble *squibble) {
  if (0 > squibble->boredom_timer) {
    squibble->chasing_whomp    = NULL;
    squibble->chasing_squibble = NULL;
    fprintf(stderr, "-- MARK -- squibble is wandering\n");
    squibble->state = SQUIBBLE_IS_WANDERING;
    squibble->boredom_timer = SQUIBBLE_BOREDOM_TIMER_DEFAULT;

    return;
  }

  if (squibble->chasing_whomp) {
    squibble_move(squibble, SQUIBBLE_CHARGE_SPEED);

    if ( squibble_can_see(squibble, &squibble->chasing_whomp->p) ) {
      squibble->boredom_timer = SQUIBBLE_BOREDOM_TIMER_DEFAULT;
      squibble->a = angle_between(&squibble->p, &squibble->chasing_whomp->p);
    }

    return;
  }

  if (squibble->chasing_squibble) {
    if ( squibble_is_near(squibble, &squibble->chasing_squibble->p) ) {
      squibble->a = angle_between(&squibble->p, &squibble->chasing_squibble->p);
      return;
    }

    squibble_move(squibble, SQUIBBLE_CHARGE_SPEED);

    if ( squibble_can_see(squibble, &squibble->chasing_squibble->p) ) {
      squibble->boredom_timer = SQUIBBLE_BOREDOM_TIMER_DEFAULT;
      squibble->a = angle_between(&squibble->p, &squibble->chasing_squibble->p);
    }

    return;
  }

  squibble_move(squibble, SQUIBBLE_CHARGE_SPEED);
}

void squibble_look_for_stuff_to_chase(struct squibble *squibble) {
  int i;

  for (i = 0; i < WHOMP_MAX; i = i + 1) {
    struct whomp *whomp = &WHOMPS[i];

    if (WHOMP_IS_DEAD == whomp->state) { continue; }

    if ( squibble_can_see(squibble, &whomp->p) ) {
      squibble->chasing_whomp = whomp;
      squibble->a = angle_between(&squibble->p, &whomp->p);
    fprintf(stderr, "-- MARK -- squibble is shocked\n");
      squibble->state = SQUIBBLE_IS_SHOCKED;
      squibble->boredom_timer = SQUIBBLE_BOREDOM_TIMER_DEFAULT;

      return;
    }
  }

  for (i = 0; i < SQUIBBLE_MAX; i = i + 1) {
    struct squibble *other = &SQUIBBLES[i];

    if (SQUIBBLE_IS_CHARGING != other->state) { continue; }

    if (squibble == other) { continue; }
    if ( squibble_can_see(squibble, &other->p) ) {
      squibble->chasing_squibble = other;
      squibble->a = angle_between(&squibble->p, &other->p);
    fprintf(stderr, "-- MARK -- squibble is shocked\n");
      squibble->state = SQUIBBLE_IS_SHOCKED;
      squibble->boredom_timer = SQUIBBLE_BOREDOM_TIMER_DEFAULT;

      return;
    }
  }
}

void squibble_wander(struct squibble *squibble) {
  int i;

  squibble_move(squibble, SQUIBBLE_WANDER_SPEED);

  if (0 > squibble->boredom_timer) {
    fprintf(stderr, "-- MARK -- squibble is looking around\n");
    squibble->state = SQUIBBLE_IS_LOOKING_AROUND;
    squibble->boredom_timer = SQUIBBLE_BOREDOM_TIMER_DEFAULT;

    return;
  }

  squibble_look_for_stuff_to_chase(squibble);
}

void squibble_look_around(struct squibble *squibble) {
  if (0 > squibble->boredom_timer) {
    squibble->a = random_angle();
    fprintf(stderr, "-- MARK -- squibble is wandering\n");
    squibble->state = SQUIBBLE_IS_WANDERING;
    squibble->boredom_timer = SQUIBBLE_BOREDOM_TIMER_DEFAULT;

    return;
  }

  squibble_look_for_stuff_to_chase(squibble);
}

void squibble_shock(struct squibble *squibble) {
  if (0 > squibble->boredom_timer) {
    fprintf(stderr, "-- MARK -- squibble is charging\n");
    squibble->state = SQUIBBLE_IS_CHARGING;
    squibble->boredom_timer = SQUIBBLE_BOREDOM_TIMER_DEFAULT;

    return;
  }
}

void (*HANDLE_SQUIBBLE_STATE[])(struct squibble *) = {
  squibble_do_nothing,  /* SQUIBBLE_IS_DEAD */
  squibble_charge,      /* SQUIBBLE_IS_CHARGING */
  squibble_wander,      /* SQUIBBLE_IS_WANDERING */
  squibble_look_around, /* SQUIBBLE_IS_LOOKING_AROUND */
  squibble_shock,       /* SQUIBBLE_IS_SHOCKED */
  squibble_do_nothing,  /* SQUIBBLE_IS_BEING_HELD */
  squibble_do_nothing   /* SQUIBBLE_IS_AIRBORNE */
};

struct player {
  double x, y;
  int    action;
};

void generate_field() {
  int x, y;

  int s = FIELD_RES / 2;
  int pattern[] = STUB_PATTERN;

  /* load stub pattern */

  for (y = 0; y < s; y = y + 1) {
    for (x = 0; x < s; x = x + 1) {
      FIELD[y][x] = pattern[y * s + x];
    }
  }

  /* copy pattern to other quadrants */

  for (y = 0; y < s; y = y + 1) {
    for (x = 0; x < s; x = x + 1) {
      int yi = 2 * s - y - 1;
      int xi = 2 * s - x - 1;

      FIELD[y][xi]  = FIELD[y][x];
      FIELD[yi][x]  = FIELD[y][x];
      FIELD[yi][xi] = FIELD[y][x];
    }
  }
}

int collides_with_field(double x, double y, double r) {
  int xi, yi;

  int y0 = floor(y - r);
  int y1 = ceil(y + r);
  int x0 = floor(x - r);
  int x1 = ceil(x + r);

  int s  = FIELD_RES / 2;

  for (yi = y0; yi < y1; yi = yi + 1) {
    for (xi = x0; xi < x1; xi = xi + 1) {
      if ( field(xi, yi) ) {
        return true;
      }
    }
  }

  return false;
}

int main(int argc, char **argv) {
  int    hres, vres, width, height, pixels, i;
  double aspect, scale;

  SDL_Surface    *sdl_surface;
  Uint32         next_frame;
  cairo_t        *cr;
  cairo_matrix_t cm_display, cm_field;

  struct player   player;
  struct whomp    *whomp    = &WHOMPS[0];
  struct squibble *squibble = &SQUIBBLES[0];

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

  /* hres = 640; vres = 480; */

  sdl_surface = SDL_SetVideoMode(hres, vres, 32, SDL_VIDEO_FLAGS);
  /* SDL_ShowCursor(false); */

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

  generate_field();
  clear_overlay();

  player.action = 0;

  for (i = 0; i < WHOMP_MAX; i = i + 1) {
    WHOMPS[i].state = WHOMP_IS_DEAD;
  }

  whomp->state = WHOMP_IS_ALIVE;
  whomp->p.x = 0;
  whomp->p.y = 0;
  whomp->squibble_held = NULL;

  for (i = 0; i < SQUIBBLE_MAX; i = i + 1) {
    SQUIBBLES[i].state = SQUIBBLE_IS_DEAD;
  }

  squibble_new();

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
      cairo_set_line_width(cr, 1/8.0);
      for (y = 0; y < 2 * s; y = y + 1) {
        for (x = 0; x < 2 * s; x = x + 1) {
          if (! FIELD[y][x]) { continue; }

          cairo_set_matrix(cr, &cm_field);
          cairo_translate(cr, x - s + 1/2.0, y - s + 1/2.0);

          cairo_move_to(cr, -1/2.0, -1/2.0);
          cairo_line_to(cr,  1/2.0, -1/2.0);
          cairo_line_to(cr,  1/2.0,  1/2.0);
          cairo_line_to(cr, -1/2.0,  1/2.0);
          cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
          cairo_fill(cr);

          cairo_move_to(cr, -1/2.0, -1/2.0);
          cairo_line_to(cr,  1/2.0, -1/2.0);
          cairo_move_to(cr, -1/2.0,  0/2.0);
          cairo_line_to(cr,  1/2.0,  0/2.0);
          cairo_move_to(cr, -1/2.0,  1/2.0);
          cairo_line_to(cr,  1/2.0,  1/2.0);
          cairo_move_to(cr, -1/2.0, -1/2.0);
          cairo_line_to(cr, -1/2.0,  0/2.0);
          cairo_move_to(cr,  1/2.0, -1/2.0);
          cairo_line_to(cr,  1/2.0,  0/2.0);
          cairo_move_to(cr,  0/2.0,  0/2.0);
          cairo_line_to(cr,  0/2.0,  1/2.0);
          cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
          cairo_stroke(cr);
        }
      }

      for (y = 0; y < 2 * s; y = y + 1) {
        for (x = 0; x < 2 * s; x = x + 1) {
          if (! OVERLAY[y][x]) { continue; }

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
      clear_overlay();

      /* draw whomp */
      cairo_set_matrix(cr, &cm_field);
      cairo_translate(cr, whomp->p.x, whomp->p.y);
      cairo_move_to(cr, -1/2.0,  0/2.0);
      cairo_line_to(cr,  1/2.0,  0/2.0);
      cairo_line_to(cr,  1/2.0,  2/2.0);
      cairo_line_to(cr, -1/2.0,  2/2.0);
      cairo_close_path(cr);
      cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
      cairo_set_line_width(cr, 1/4.0);
      cairo_stroke_preserve(cr);
      cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
      cairo_fill(cr);

      /* draw squibble */

      { /* Squibbles */
        int i;

        for (i = 0; i < SQUIBBLE_MAX; i = i + 1) {
          struct squibble *squibble = &SQUIBBLES[i];
          double r = SQUIBBLE_RADIUS;

          if (SQUIBBLE_IS_DEAD == squibble->state) { continue; }

          cairo_set_matrix(cr, &cm_field);
          cairo_translate(cr, squibble->p.x, squibble->p.y);
          while (0) cairo_rotate(cr, squibble->a);
          cairo_arc(cr, 0, r, r, 0, 2 * M_PI);
          cairo_close_path(cr);
          cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
          cairo_set_line_width(cr, 1/4.0);
          cairo_stroke_preserve(cr);
          cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
          cairo_fill(cr);
          while (0) cairo_move_to(cr, 0, 0);
          while (0) cairo_line_to(cr, 1/2.0, 0);
          while (0) cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
          while (0) cairo_stroke(cr);

          if (SQUIBBLE_IS_SHOCKED == squibble->state) {
            double a = 0;

            for (a = 0; a < 2 * M_PI; a = a + M_PI / 6.0) {
              double x = r * cos(a);
              double y = r * sin(a);
              cairo_move_to(cr, 5/4.0 * x, 5/4.0 * y + r);
              cairo_line_to(cr, 2/1.0 * x, 2/1.0 * y + r);
            }
            cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
            cairo_set_line_width(cr, 1/4.0);
            cairo_stroke_preserve(cr);
            cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
            cairo_set_line_width(cr, 1/8.0);
            cairo_stroke(cr);
          }
        }
      }

      /* draw field */
      cairo_set_line_width(cr, 1/8.0);
      for (y = 0; y < 2 * s; y = y + 1) {
        for (x = 0; x < 2 * s; x = x + 1) {
          if (! FIELD[y][x]) { continue; }

          cairo_set_matrix(cr, &cm_field);
          cairo_translate(cr, x - s + 1/2.0, y - s + 1/2.0);

          cairo_move_to(cr, -1/2.0,  1/2.0);
          cairo_line_to(cr,  1/2.0,  1/2.0);
          cairo_line_to(cr,  1/2.0,  3/2.0);
          cairo_line_to(cr, -1/2.0,  3/2.0);
          cairo_close_path(cr);
          cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
          cairo_stroke_preserve(cr);
          cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
          cairo_fill(cr);
        }
      }
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
            if (SDLK_n == event.key.keysym.sym) {
              squibble_new();
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
        double x0 = whomp->p.x;
        double y0 = whomp->p.y;
        double r  = WHOMP_RADIUS;

        double x = x0 + player.x * WHOMP_SPEED / FRAMERATE;
        double y = y0 + player.y * WHOMP_SPEED / FRAMERATE;

        if ( collides_with_field(x, y0, r) ) {
          x = x0;
        }
        if ( collides_with_field(x0, y, r) ) {
          y = y0;
        }
        if ( collides_with_field(x, y, r) ) {
          x = x0;
          y = y0;
        }

        whomp->p.x = x;
        whomp->p.y = y;

        if (whomp->squibble_held) {
          struct squibble *squibble = whomp->squibble_held;

          squibble->p.x = whomp->p.x;
          squibble->p.y = whomp->p.y;
        } else {
          if (1 == player.action) {
            double x = squibble->p.x - whomp->p.x;
            double y = squibble->p.y - whomp->p.y;
            double d = sqrt(x * x + y * y);

            if (1 > d) {
              squibble->p.x  = whomp->p.x;
              squibble->p.y  = whomp->p.y;
              squibble->state = SQUIBBLE_IS_BEING_HELD;

              whomp->squibble_held = squibble;
            }
          }
        }
      }

      { /* Squibbles */
        int i;

        for (i = 0; i < SQUIBBLE_MAX; i = i + 1) {
          struct squibble *squibble = &SQUIBBLES[i];

          if (SQUIBBLE_IS_DEAD == squibble->state) { continue; }

          squibble->boredom_timer = squibble->boredom_timer - 1.0 / FRAMERATE;

          HANDLE_SQUIBBLE_STATE[squibble->state](squibble);
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
