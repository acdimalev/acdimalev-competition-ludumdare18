#include <cairo.h>

#include "path.h"

void path_segment_move_render(
    cairo_t *cr,
    struct path_segment_move *segment
  ) {
  cairo_move_to(cr, segment->p.x, segment->p.y);
}

void path_segment_line_render(
    cairo_t *cr,
    struct path_segment_line *segment
  ) {
  cairo_line_to(cr, segment->p.x, segment->p.y);
}

void path_segment_curve_render(
    cairo_t *cr,
    struct path_segment_curve *segment
  ) {
  cairo_curve_to(cr,
      segment->p1.x, segment->p1.y,
      segment->p2.x, segment->p2.y,
      segment->p3.x, segment->p3.y
    );
}

void path_segment_render(
    cairo_t *cr,
    struct path_segment *segment
  ) {
  switch (segment->type) {
    case PATH_MOVE:
      path_segment_move_render(cr,
          (struct path_segment_move *) segment);
      break;
    case PATH_LINE:
      path_segment_line_render(cr,
          (struct path_segment_line *) segment);
      break;
    case PATH_CURVE:
      path_segment_curve_render(cr,
          (struct path_segment_curve *) segment);
      break;
    case PATH_CLOSE:
      cairo_close_path(cr);
      break;
  }
}

