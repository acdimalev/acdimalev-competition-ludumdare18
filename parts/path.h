#include <cairo.h>

struct pos {
  double x, y;
};

enum path_segment_type {
  PATH_MOVE,
  PATH_LINE,
  PATH_CURVE,
  PATH_CLOSE
};

struct path_segment {
  struct path_segment      *prev, *next;
  enum   path_segment_type type;
};

struct path_segment_move {
  struct path_segment base;
  struct pos  p;
};

struct path_segment_line {
  struct path_segment base;
  struct pos  p;
};

struct path_segment_curve {
  struct path_segment base;
  struct pos  p1, p2, p3;
};

struct path {
  struct path              *prev, *next;
  struct path_segment_move segments;
  int    color, fill;
};

void path_segment_move_render(
    cairo_t *cr,
    struct path_segment_move *segment
  );

void path_segment_line_render(
    cairo_t *cr,
    struct path_segment_line *segment
  );

void path_segment_curve_render(
    cairo_t *cr,
    struct path_segment_curve *segment
  );

void path_segment_render(
    cairo_t *cr,
    struct path_segment *segment
  );
