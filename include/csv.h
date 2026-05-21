
#ifndef RS_CSV_H
#define RS_CSV_H

#include "point.h"
#include <stddef.h>

int rs_csv_load(const char *path,
                rs_point_t **out_pts,
                size_t *out_n,
                size_t *out_dim);

#endif
