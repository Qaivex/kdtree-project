
#ifndef RS_VIZ_H
#define RS_VIZ_H

#include "point.h"
#include <stddef.h>

typedef struct
{
    unsigned char r;
    unsigned char g;
    unsigned char b;
} rs_rgb_t;

rs_rgb_t rs_viz_color_for_cluster(int cluster_id);

void rs_viz_print_ansi(const rs_point_t *points,
                       size_t n,
                       int width,
                       int height);

#endif
