
#ifndef RS_KDTREE_H
#define RS_KDTREE_H

#include "point.h"
#include <stdbool.h>

typedef struct rs_kdnode
{
    rs_point_t point;
    int axis;
    struct rs_kdnode *left;
    struct rs_kdnode *right;
} rs_kdnode_t;

typedef struct
{
    rs_kdnode_t *root;
    size_t dim;
    size_t size;
} rs_kdtree_t;

rs_kdtree_t *rs_kdtree_create(size_t dim);

void rs_kdtree_destroy(rs_kdtree_t *tree);

rs_kdtree_t *rs_kdtree_build(const rs_point_t *points, size_t n, size_t dim);

bool rs_kdtree_insert(rs_kdtree_t *tree, const rs_point_t *p);

bool rs_kdtree_remove(rs_kdtree_t *tree, const rs_point_t *p);

bool rs_kdtree_nearest(const rs_kdtree_t *tree,
                       const rs_point_t *query,
                       rs_point_t *out,
                       double *out_sqdist);

size_t rs_kdtree_range(const rs_kdtree_t *tree,
                       const rs_point_t *query,
                       double radius,
                       rs_point_t *out,
                       size_t cap);

size_t rs_kdtree_size(const rs_kdtree_t *tree);

#endif
