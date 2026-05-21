
#ifndef RS_CMEANS_H
#define RS_CMEANS_H

#include "point.h"
#include <stddef.h>

typedef struct
{
    size_t c;
    double m;
    size_t max_iter;
    double tol;
    unsigned seed;
} rs_cmeans_params_t;

typedef struct
{
    size_t n;
    size_t c;
    size_t dim;
    double *membership;
    rs_point_t *centroids;
    double objective;
    size_t iterations;
} rs_cmeans_result_t;

rs_cmeans_params_t rs_cmeans_default_params(size_t c);

rs_cmeans_result_t *rs_cmeans_run(const rs_point_t *points,
                                  size_t n,
                                  size_t dim,
                                  const rs_cmeans_params_t *params);

void rs_cmeans_result_free(rs_cmeans_result_t *res);

void rs_cmeans_hard_assign(const rs_cmeans_result_t *res,
                           rs_point_t *out_points);

size_t rs_cmeans_elbow(const rs_point_t *points,
                       size_t n,
                       size_t dim,
                       size_t c_max,
                       double *out_J);

#endif
