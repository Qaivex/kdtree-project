
#ifndef RS_VIZ_H
#define RS_VIZ_H

#include "point.h"
#include <stddef.h>

// 24-битный RGB-цвет
typedef struct
{
    unsigned char r; // канал красного, 0..255
    unsigned char g; // канал зелёного, 0..255
    unsigned char b; // канал синего, 0..255
} rs_rgb_t;

// цвет из палитры по идентификатору кластера
rs_rgb_t rs_viz_color_for_cluster(int cluster_id);

// печать цветной ASCII-карты точек в stdout
void rs_viz_print_ansi(const rs_point_t *points,
                       size_t n,
                       int width,
                       int height);

#endif
