
#ifndef RS_DBSCAN_H
#define RS_DBSCAN_H

#include "point.h"
#include <stddef.h>

// специальное значение cluster_id для шумовых точек
#define RS_DBSCAN_NOISE (-1)

// параметры алгоритма DBSCAN
typedef struct
{
    double eps;     // радиус ε-окрестности
    size_t min_pts; // минимальное число соседей для ядра
} rs_dbscan_params_t;

// результат кластеризации DBSCAN
typedef struct
{
    size_t n;          // число точек
    int *cluster_ids;  // массив cluster_id для каждой точки
    size_t n_clusters; // итоговое число кластеров
    size_t n_noise;    // число точек, попавших в шум
} rs_dbscan_result_t;

// запуск DBSCAN с ε-окрестностью через k-d дерево
rs_dbscan_result_t *rs_dbscan_run(rs_point_t *points,
                                  size_t n,
                                  size_t dim,
                                  const rs_dbscan_params_t *params);

// освобождение результата кластеризации DBSCAN
void rs_dbscan_result_free(rs_dbscan_result_t *res);

#endif
