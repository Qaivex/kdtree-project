
#ifndef RS_CMEANS_H
#define RS_CMEANS_H

#include "point.h"
#include <stddef.h>

// параметры алгоритма Fuzzy C-means
typedef struct
{
    size_t c;        // число кластеров
    double m;        // коэффициент нечёткости (обычно 2.0)
    size_t max_iter; // верхняя граница числа итераций
    double tol;      // допуск для критерия останова
    unsigned seed;   // зерно ГПСЧ для воспроизводимости
} rs_cmeans_params_t;

// результат кластеризации Fuzzy C-means
typedef struct
{
    size_t n;              // число точек
    size_t c;              // число кластеров
    size_t dim;            // размерность точек
    double *membership;    // матрица принадлежности U[n*c], row-major
    rs_point_t *centroids; // центры кластеров (массив длины c)
    double objective;      // итоговое значение целевой функции J_m
    size_t iterations;     // реально выполненное число итераций
} rs_cmeans_result_t;

// параметры FCM по умолчанию: m=2, max_iter=200, tol=1e-5, seed=42
rs_cmeans_params_t rs_cmeans_default_params(size_t c);

// запуск алгоритма Fuzzy C-means
rs_cmeans_result_t *rs_cmeans_run(const rs_point_t *points,
                                  size_t n,
                                  size_t dim,
                                  const rs_cmeans_params_t *params);

// освобождение результата кластеризации FCM
void rs_cmeans_result_free(rs_cmeans_result_t *res);

// присваивание точкам жёсткой метки по максимуму принадлежности
void rs_cmeans_hard_assign(const rs_cmeans_result_t *res,
                           rs_point_t *out_points);

// метод локтя: возвращает рекомендованное число кластеров
size_t rs_cmeans_elbow(const rs_point_t *points,
                       size_t n,
                       size_t dim,
                       size_t c_max,
                       double *out_J);

#endif
