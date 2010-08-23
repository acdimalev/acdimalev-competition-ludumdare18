#define FRAMERATE 30

#define FIELD_RES (1 << 4)

enum body_state {
  BODY_IS_CONTROLLED,
  BODY_IS_SLIDING,
  BODY_IS_AIRBORNE
};

struct vec {
  double x, y, z;
};

struct body {
  enum   body_state state;
  struct vec        p, v;
  double            r;
};

extern int FIELD[FIELD_RES][FIELD_RES];
extern struct vec VEC0;

int  field(int x, int y);
double dist(
    struct vec *a,
    struct vec *b
    );
void body_animate(struct body *body);
