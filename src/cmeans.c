
#include "cmeans.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* псевдослучайный генератор xorshift32 */
static unsigned rs_xs_next(unsigned *state)
{
    unsigned x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    return x;
}

/* псевдослучайное число double в полуинтервале [0, 1) */
static double rs_xs_unit(unsigned *state)
{
    return (double)(rs_xs_next(state) & 0x00FFFFFFu) / (double)0x01000000;
}

rs_cmeans_params_t rs_cmeans_default_params(size_t c)
{
    rs_cmeans_params_t p;
    p.c = c;
    p.m = 2.0;
    p.max_iter = 200;
    p.tol = 1e-5;
    p.seed = 42u;
    return p;
}

/* инициализация матрицы принадлежности FCM случайными нормированными значениями */
static void cm_init_membership(double *U, size_t n, size_t c, unsigned *rng)
{
    for (size_t i = 0; i < n; ++i)
    {
        double sum = 0.0;
        for (size_t j = 0; j < c; ++j)
        {
            
            double r = rs_xs_unit(rng) + 0.001;
            U[i * c + j] = r;
            sum += r;
        }
        
        for (size_t j = 0; j < c; ++j)
        {
            U[i * c + j] /= sum;
        }
    }
}

/* пересчёт центроидов Fuzzy C-means: v_j = sum(u_ij^m * x_i) / sum(u_ij^m) */
static void cm_update_centroids(const rs_point_t *X,
                                const double *U,
                                size_t n,
                                size_t c,
                                size_t dim,
                                double m,
                                rs_point_t *V)
{
    for (size_t j = 0; j < c; ++j)
    {
        double w_sum = 0.0;
        double acc[RS_MAX_DIM];
        for (size_t d = 0; d < dim; ++d)
        {
            acc[d] = 0.0;
        }
        for (size_t i = 0; i < n; ++i)
        {
            double w = pow(U[i * c + j], m);
            w_sum += w;
            for (size_t d = 0; d < dim; ++d)
            {
                acc[d] += w * X[i].coord[d];
            }
        }
        if (w_sum < 1e-18)
        {
            
            continue;
        }
        for (size_t d = 0; d < dim; ++d)
        {
            V[j].coord[d] = acc[d] / w_sum;
        }
        V[j].dim = dim;
        V[j].cluster_id = (int)j;
        V[j].index = j;
    }
}

/* пересчёт матрицы принадлежности Fuzzy C-means по формуле Бездека */
static void cm_update_membership(const rs_point_t *X,
                                 const rs_point_t *V,
                                 size_t n,
                                 size_t c,
                                 double m,
                                 double *U)
{
    
    double exponent = 2.0 / (m - 1.0);
    for (size_t i = 0; i < n; ++i)
    {
        
        double d_iv[64];
        int zero_at = -1;
        for (size_t j = 0; j < c; ++j)
        {
            double d = rs_point_dist(&X[i], &V[j]);
            d_iv[j] = d;
            if (d < 1e-12)
            {
                zero_at = (int)j;
            }
        }
        if (zero_at >= 0)
        {
            for (size_t j = 0; j < c; ++j)
            {
                U[i * c + j] = (j == (size_t)zero_at) ? 1.0 : 0.0;
            }
            continue;
        }
        
        for (size_t j = 0; j < c; ++j)
        {
            double sum = 0.0;
            for (size_t k = 0; k < c; ++k)
            {
                double ratio = d_iv[j] / d_iv[k];
                sum += pow(ratio, exponent);
            }
            U[i * c + j] = 1.0 / sum;
        }
    }
}

/* вычисление целевой функции J_m алгоритма Fuzzy C-means */
static double cm_objective(const rs_point_t *X,
                           const rs_point_t *V,
                           const double *U,
                           size_t n,
                           size_t c,
                           double m)
{
    double J = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < c; ++j)
        {
            double d2 = rs_point_sqdist(&X[i], &V[j]);
            J += pow(U[i * c + j], m) * d2;
        }
    }
    return J;
}

/* max-норма разности двух матриц — критерий останова FCM */
static double cm_inf_diff(const double *A, const double *B, size_t n_total)
{
    double m = 0.0;
    for (size_t i = 0; i < n_total; ++i)
    {
        double d = A[i] - B[i];
        if (d < 0.0)
        {
            d = -d;
        }
        if (d > m)
        {
            m = d;
        }
    }
    return m;
}

rs_cmeans_result_t *rs_cmeans_run(const rs_point_t *points,
                                  size_t n,
                                  size_t dim,
                                  const rs_cmeans_params_t *params)
{
    if (points == NULL || n == 0 || params == NULL || params->c < 2
        || params->c > 64 || params->m <= 1.0 || dim == 0 || dim > RS_MAX_DIM)
    {
        return NULL;
    }
    size_t c = params->c;
    rs_cmeans_result_t *res = (rs_cmeans_result_t *)calloc(1, sizeof(*res));
    if (res == NULL)
    {
        return NULL;
    }
    res->n = n;
    res->c = c;
    res->dim = dim;
    res->membership = (double *)malloc(sizeof(double) * n * c);
    res->centroids = (rs_point_t *)calloc(c, sizeof(rs_point_t));
    double *U_old = (double *)malloc(sizeof(double) * n * c);
    if (res->membership == NULL || res->centroids == NULL || U_old == NULL)
    {
        rs_cmeans_result_free(res);
        free(U_old);
        return NULL;
    }
    
    unsigned rng = params->seed ? params->seed : 1u;
    for (size_t j = 0; j < c; ++j)
    {
        size_t idx = rs_xs_next(&rng) % n;
        rs_point_copy(&res->centroids[j], &points[idx]);
        res->centroids[j].cluster_id = (int)j;
    }
    
    cm_init_membership(res->membership, n, c, &rng);
    
    size_t it = 0;
    for (; it < params->max_iter; ++it)
    {
        memcpy(U_old, res->membership, sizeof(double) * n * c);
        cm_update_centroids(points, res->membership, n, c, dim, params->m, res->centroids);
        cm_update_membership(points, res->centroids, n, c, params->m, res->membership);
        double diff = cm_inf_diff(res->membership, U_old, n * c);
        if (diff < params->tol)
        {
            ++it;
            break;
        }
    }
    res->iterations = it;
    res->objective = cm_objective(points, res->centroids, res->membership, n, c, params->m);
    free(U_old);
    return res;
}

void rs_cmeans_result_free(rs_cmeans_result_t *res)
{
    if (res == NULL)
    {
        return;
    }
    free(res->membership);
    free(res->centroids);
    free(res);
}

/* жёсткое назначение кластеров FCM: cluster_id = argmax_j u_ij */
void rs_cmeans_hard_assign(const rs_cmeans_result_t *res,
                           rs_point_t *out_points)
{
    if (res == NULL || out_points == NULL)
    {
        return;
    }
    for (size_t i = 0; i < res->n; ++i)
    {
        size_t best = 0;
        double best_u = res->membership[i * res->c + 0];
        for (size_t j = 1; j < res->c; ++j)
        {
            double u = res->membership[i * res->c + j];
            if (u > best_u)
            {
                best_u = u;
                best = j;
            }
        }
        out_points[i].cluster_id = (int)best;
    }
}

/* метод локтя: выбор c по максимуму второй разности J(c-1)-2J(c)+J(c+1) */
size_t rs_cmeans_elbow(const rs_point_t *points,
                       size_t n,
                       size_t dim,
                       size_t c_max,
                       double *out_J)
{
    if (points == NULL || n < 4 || c_max < 3 || out_J == NULL)
    {
        return 2;
    }
    
    for (size_t c = 2; c <= c_max; ++c)
    {
        rs_cmeans_params_t p = rs_cmeans_default_params(c);
        rs_cmeans_result_t *r = rs_cmeans_run(points, n, dim, &p);
        if (r == NULL)
        {
            out_J[c] = 0.0;
            continue;
        }
        out_J[c] = r->objective;
        rs_cmeans_result_free(r);
    }
    
    size_t best_c = 2;
    double best_score = -1.0;
    for (size_t c = 3; c + 1 <= c_max; ++c)
    {
        double s = out_J[c - 1] - 2.0 * out_J[c] + out_J[c + 1];
        if (s > best_score)
        {
            best_score = s;
            best_c = c;
        }
    }
    return best_c;
}
