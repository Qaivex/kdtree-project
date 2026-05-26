
#include "point.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

// квадрат евклидова расстояния между двумя точками
double rs_point_sqdist(const rs_point_t *a, const rs_point_t *b)
{
    double sum = 0.0; // накопитель суммы квадратов разностей
    for (size_t i = 0; i < a->dim; ++i)
    {
        double d = a->coord[i] - b->coord[i]; // разность по оси i
        sum += d * d;
    }
    return sum;
}

// евклидово расстояние через корень от квадрата
double rs_point_dist(const rs_point_t *a, const rs_point_t *b)
{
    return sqrt(rs_point_sqdist(a, b));
}

// побитовое копирование точки через memcpy
void rs_point_copy(rs_point_t *dst, const rs_point_t *src)
{
    memcpy(dst, src, sizeof(*src));
}

// печать точки в формате (x, y[, ...])  cluster=...  index=...
void rs_point_print(const rs_point_t *p)
{
    putchar('(');
    for (size_t i = 0; i < p->dim; ++i)
    {
        if (i > 0)
        {
            fputs(", ", stdout); // разделитель между координатами
        }
        printf("%.4f", p->coord[i]);
    }
    printf(")  cluster=%d  index=%zu\n", p->cluster_id, p->index);
}
