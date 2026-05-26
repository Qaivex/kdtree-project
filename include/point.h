
#ifndef RS_POINT_H
#define RS_POINT_H

#include <stddef.h>

// максимальная поддерживаемая размерность пространства
#define RS_MAX_DIM 8

// точка в k-мерном пространстве
typedef struct
{
    double coord[RS_MAX_DIM]; // координаты по осям
    size_t dim;               // фактическая размерность точки
    int    cluster_id;        // идентификатор кластера или -1 (шум)
    size_t index;             // индекс точки во входном массиве
} rs_point_t;

// квадрат евклидова расстояния между двумя точками
double rs_point_sqdist(const rs_point_t *a, const rs_point_t *b);

// евклидово расстояние между двумя точками
double rs_point_dist(const rs_point_t *a, const rs_point_t *b);

// побайтовое копирование точки
void rs_point_copy(rs_point_t *dst, const rs_point_t *src);

// печать точки в stdout в читаемом виде
void rs_point_print(const rs_point_t *p);

#endif
