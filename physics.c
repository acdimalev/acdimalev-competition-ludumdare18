#include <math.h>
#include <stdio.h>

#include "physics.h"

#define BOUNCE_FRICTION (1/2.0)
#define SLIDING_FRICTION (32.0)
#define GRAVITY -16.0
#define LANDING_VELOCITY (1/8.0)
#define SLIDING_VELOCITY 4.0
#define BODY_IS_ON_THE_GROUND_ERROR_ALLOWANCE (1/256.0)

#define true -1
#define false 0

int FIELD[FIELD_RES][FIELD_RES];

struct vec VEC0 = {0, 0, 0};

int field(int x, int y) {
  int s = FIELD_RES / 2;

  if (x < -s) { return 0; }
  if (x >= s) { return 0; }
  if (y < -s) { return 0; }
  if (y >= s) { return 0; }

  return FIELD[y + s][x + s];
}

double dist(
    struct vec *a,
    struct vec *b
    ) {
  double x = b->x - a->x;
  double y = b->y - a->y;
  double z = b->z - a->z;

  return sqrt(x * x + y * y + z * z);
}

double dot(
    struct vec *a,
    struct vec *b
    ) {
  double x = a->x * b->x;
  double y = a->y * b->y;
  double z = a->z * b->z;

  return x + y + z;
}

void mulvc(
    struct vec *a,
    struct vec *b,
    double c
    ) {
  a->x = c * b->x;
  a->y = c * b->y;
  a->z = c * b->z;
}

void ldv(
    struct vec *a,
    struct vec *b
    ) {
  a->x = b->x;
  a->y = b->y;
  a->z = b->z;
}

void addv(
    struct vec *a,
    struct vec *b,
    struct vec *c
    ) {
  a->x = b->x + c->x;
  a->y = b->y + c->y;
  a->z = b->z + c->z;
}

void subv(
    struct vec *a,
    struct vec *b,
    struct vec *c
    ) {
  a->x = b->x - c->x;
  a->y = b->y - c->y;
  a->z = b->z - c->z;
}

double nearest(double a, double b) {
  if (a > b) {
    return a;
  } else if (a + 1 < b) {
    return a + 1;
  } else {
    return b;
  }
}

int body_collides_with_field(
    double *offset,
    struct vec  *normal,
    struct body *body
    ) {
  double d = body->r;
  int x0 = floor(body->p.x - 1/2.0);
  int y0 = floor(body->p.y - 1/2.0);
  int n = ceil(2 * body->r) + 1;

  struct vec c;
  int xi, yi, bits;

  for (xi = 0; xi < n; xi = xi + 1) {
    int x = xi + x0;

    for (yi = 0; yi < n; yi = yi + 1) {
      int y = yi + y0;
      int z = field(x, y) - 1;

      struct vec p;
      double dn;

      p.x = nearest(x, body->p.x);
      p.y = nearest(y, body->p.y);
      p.z = nearest(z, body->p.z);

      dn = dist(&p, &body->p);

      if (d > dn) {
        subv(&c, &p, &body->p);
        d = dn;
      }
    }
  }

  if (d == body->r) { return false; }

  *offset = body->r - d;
  mulvc(normal, &c, -1/d);

  return true;
}

int body_is_on_the_ground(
    struct body *body
    ) {
  int x = floor(body->p.x);
  int y = floor(body->p.y);

  double z = body->p.z - body->r - field(x, y);

  return (z < BODY_IS_ON_THE_GROUND_ERROR_ALLOWANCE);
}

void align_body_to_the_ground(
    struct body *body
    ) {
  int x = floor(body->p.x);
  int y = floor(body->p.y);

  body->p.z = field(x, y) + body->r;
}

void body_animate(struct body *body) {
  if (BODY_IS_CONTROLLED == body->state) {
    double x0 = body->p.x;
    double y0 = body->p.y;
    double dx = body->v.x / FRAMERATE;
    double dy = body->v.y / FRAMERATE;

    struct vec n;
    double o;

    body->p.x = x0 + dx;
    body->p.y = y0;
    if ( body_collides_with_field(&o, &n, body) ) {
      dx = 0;
    }

    body->p.x = x0;
    body->p.y = y0 + dy;
    if ( body_collides_with_field(&o, &n, body) ) {
      dy = 0;
    }

    body->p.x = x0 + dx;
    body->p.y = y0 + dy;
    if ( body_collides_with_field(&o, &n, body) ) {
      body->p.x = x0;
      body->p.y = y0;
    }

  } else { /* ! BODY_IS_CONTROLLED */
    struct vec dp, n;
    double o;

    mulvc(&dp, &body->v, 1.0/FRAMERATE);
    addv(&body->p, &body->p, &dp);

    if ( body_collides_with_field(&o, &n, body) ) {
      struct vec dp, dv;
      double d, f;

      mulvc( &dp, &n, (2 - BOUNCE_FRICTION) * o );
      addv(&body->p, &body->p, &dp);

      d = dot(&body->v, &n);
      /* use this for explosions and contact damage */
      /* f = d / body->m; */
      mulvc( &dv, &n, -(2 - BOUNCE_FRICTION) * d );
      addv(&body->v, &body->v, &dv);
    }

    if (BODY_IS_SLIDING == body->state) {
      double d = dist(&VEC0, &body->v);

      if (d < SLIDING_FRICTION / FRAMERATE) {
        body->v.x = 0;
        body->v.y = 0;
      } else {
        struct vec dv;

        d = -SLIDING_FRICTION / FRAMERATE / d;
        mulvc(&dv, &body->v, d);
        addv(&body->v, &body->v, &dv);
      }

    } else { /* BODY_IS_AIRBORNE */
      body->v.z = body->v.z + GRAVITY / FRAMERATE;
    }
  }

  if ( ! body_is_on_the_ground(body) ||
      LANDING_VELOCITY < fabs(body->v.z) ) {
    if (BODY_IS_AIRBORNE != body->state) {
      fprintf(stderr, "-- MARK -- body is airborne\n");
    }
    body->state = BODY_IS_AIRBORNE;
  } else {
    align_body_to_the_ground(body);
    body->v.z = 0;

    if ( SLIDING_VELOCITY <= dist(&VEC0, &body->v) ) {
      if (BODY_IS_SLIDING != body->state) {
        fprintf(stderr, "-- MARK -- body is sliding\n");
      }
      body->state = BODY_IS_SLIDING;
    } else {
      if (BODY_IS_CONTROLLED != body->state) {
        fprintf(stderr, "-- MARK -- body is controlled\n");
      }
      body->state = BODY_IS_CONTROLLED;
    }
  }
}
