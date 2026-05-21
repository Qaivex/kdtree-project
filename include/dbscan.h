
#ifndef RS_DBSCAN_H
#define RS_DBSCAN_H

#include "point.h"
#include <stddef.h>

#define RS_DBSCAN_NOISE (-1)

typedef struct
{
    double eps;
    size_t min_pts;
} rs_dbscan_params_t;

typedef struct
{
    size_t n;
    int *cluster_ids;
    size_t n_clusters;
    size_t n_noise;
} rs_dbscan_result_t;

rs_dbscan_result_t *rs_dbscan_run(rs_point_t *points,
                                  size_t n,
                                  size_t dim,
                                  const rs_dbscan_params_t *params);

void rs_dbscan_result_free(rs_dbscan_result_t *res);

#endif
