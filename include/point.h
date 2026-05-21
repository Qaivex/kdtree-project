
#ifndef RS_POINT_H
#define RS_POINT_H

#include <stddef.h>

#define RS_MAX_DIM 8

typedef struct
{
    double coord[RS_MAX_DIM];
    size_t dim;
    int    cluster_id;
    size_t index;
} rs_point_t;

double rs_point_sqdist(const rs_point_t *a, const rs_point_t *b);

double rs_point_dist(const rs_point_t *a, const rs_point_t *b);

void rs_point_copy(rs_point_t *dst, const rs_point_t *src);

void rs_point_print(const rs_point_t *p);

#endif
