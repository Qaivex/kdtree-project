
#include "point.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

double rs_point_sqdist(const rs_point_t *a, const rs_point_t *b)
{
    double sum = 0.0;
    
    for (size_t i = 0; i < a->dim; ++i)
    {
        double d = a->coord[i] - b->coord[i];
        sum += d * d;
    }
    return sum;
}

double rs_point_dist(const rs_point_t *a, const rs_point_t *b)
{
    return sqrt(rs_point_sqdist(a, b));
}

void rs_point_copy(rs_point_t *dst, const rs_point_t *src)
{
    memcpy(dst, src, sizeof(*src));
}

void rs_point_print(const rs_point_t *p)
{
    putchar('(');
    for (size_t i = 0; i < p->dim; ++i)
    {
        if (i > 0)
        {
            fputs(", ", stdout);
        }
        printf("%.4f", p->coord[i]);
    }
    printf(")  cluster=%d  index=%zu\n", p->cluster_id, p->index);
}
