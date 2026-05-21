
#include "csv.h"
#include "kdtree.h"
#include "cmeans.h"
#include "dbscan.h"
#include "viz.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

typedef struct
{
    const char *csv_path;
    const char *op;
    const char *op_arg1;
    const char *op_arg2;
    int no_ansi;
    int width;
    int height;
} rs_args_t;

static void print_usage(const char *prog)
{
    fprintf(stderr,
            "kdtree — k-d дерево, Fuzzy C-means и DBSCAN\n\n"
            "Использование:\n"
            "  %s <csv> <операция> [аргументы] [опции]\n\n"
            "Операции:\n"
            "  -kd_insert    X,Y[,Z]\n"
            "  -kd_nearest   X,Y[,Z]\n"
            "  -cmeans       c\n"
            "  -cmeans_elbow c_max\n"
            "  -dbscan       eps,minPts\n\n"
            "Опции вывода:\n"
            "  --no-ansi     не печатать ANSI-карту в stdout\n"
            "  --width N     ширина холста\n"
            "  --height N    высота холста\n",
            prog);
}

static size_t parse_coords(const char *s, double *out, size_t max_dim)
{
    size_t k = 0;
    const char *p = s;
    while (*p && k < max_dim)
    {
        char *endp = NULL;
        double v = strtod(p, &endp);
        if (endp == p)
        {
            break;
        }
        out[k++] = v;
        p = endp;
        while (*p == ',' || isspace((unsigned char)*p))
        {
            ++p;
        }
    }
    return k;
}

static int parse_args(int argc, char **argv, rs_args_t *out)
{
    memset(out, 0, sizeof(*out));
    out->width = 80;
    out->height = 30;
    if (argc < 3)
    {
        return -1;
    }
    out->csv_path = argv[1];
    out->op = argv[2];
    int idx = 3;
    
    if (idx < argc && strncmp(argv[idx], "--", 2) != 0)
    {
        out->op_arg1 = argv[idx++];
    }
    if (idx < argc && strncmp(argv[idx], "--", 2) != 0)
    {
        out->op_arg2 = argv[idx++];
    }
    
    while (idx < argc)
    {
        const char *a = argv[idx];
        if (strcmp(a, "--no-ansi") == 0)
        {
            out->no_ansi = 1;
        }
        else if (strcmp(a, "--width") == 0 && idx + 1 < argc)
        {
            out->width = atoi(argv[++idx]);
        }
        else if (strcmp(a, "--height") == 0 && idx + 1 < argc)
        {
            out->height = atoi(argv[++idx]);
        }
        else
        {
            fprintf(stderr, "Неизвестный аргумент: %s\n", a);
            return -1;
        }
        ++idx;
    }
    return 0;
}

static int load_or_die(const rs_args_t *a,
                       rs_point_t **pts, size_t *n, size_t *dim)
{
    if (rs_csv_load(a->csv_path, pts, n, dim) != 0)
    {
        return -1;
    }
    fprintf(stderr, "Загружено %zu точек размерности %zu\n", *n, *dim);
    return 0;
}

static int op_kd_insert(const rs_args_t *a)
{
    rs_point_t *pts = NULL;
    size_t n = 0, dim = 0;
    if (load_or_die(a, &pts, &n, &dim) != 0)
    {
        return 1;
    }
    if (a->op_arg1 == NULL)
    {
        fprintf(stderr, "kd_insert: нужны координаты\n");
        free(pts);
        return 1;
    }
    double xq[RS_MAX_DIM];
    size_t k = parse_coords(a->op_arg1, xq, RS_MAX_DIM);
    if (k != dim)
    {
        fprintf(stderr, "kd_insert: ожидаю %zu координат, получил %zu\n", dim, k);
        free(pts);
        return 1;
    }
    rs_kdtree_t *t = rs_kdtree_build(pts, n, dim);
    rs_point_t q;
    memset(&q, 0, sizeof(q));
    for (size_t d = 0; d < dim; ++d)
    {
        q.coord[d] = xq[d];
    }
    q.dim = dim;
    q.cluster_id = -1;
    q.index = n;
    bool ok = rs_kdtree_insert(t, &q);
    printf("Вставка: %s\n", ok ? "ok" : "ошибка");
    printf("Размер дерева: %zu\n", rs_kdtree_size(t));
    rs_kdtree_destroy(t);
    free(pts);
    return 0;
}

static int op_kd_nearest(const rs_args_t *a)
{
    rs_point_t *pts = NULL;
    size_t n = 0, dim = 0;
    if (load_or_die(a, &pts, &n, &dim) != 0)
    {
        return 1;
    }
    if (a->op_arg1 == NULL)
    {
        fprintf(stderr, "kd_nearest: нужны координаты\n");
        free(pts);
        return 1;
    }
    double xq[RS_MAX_DIM];
    size_t k = parse_coords(a->op_arg1, xq, RS_MAX_DIM);
    if (k != dim)
    {
        fprintf(stderr, "kd_nearest: ожидаю %zu координат, получил %zu\n", dim, k);
        free(pts);
        return 1;
    }
    rs_kdtree_t *t = rs_kdtree_build(pts, n, dim);
    rs_point_t q;
    memset(&q, 0, sizeof(q));
    for (size_t d = 0; d < dim; ++d)
    {
        q.coord[d] = xq[d];
    }
    q.dim = dim;
    rs_point_t best;
    double sq = 0.0;
    if (rs_kdtree_nearest(t, &q, &best, &sq))
    {
        printf("Ближайший сосед: ");
        rs_point_print(&best);
        printf("Расстояние: %.6f\n", (sq > 0 ? sqrt(sq) : 0.0));
    }
    else
    {
        printf("Дерево пусто\n");
    }
    rs_kdtree_destroy(t);
    free(pts);
    return 0;
}

/* экспорт результатов кластеризации: ANSI-визуализация в stdout */
static void emit_results(const rs_args_t *a,
                         rs_point_t *pts, size_t n)
{
    if (!a->no_ansi)
    {
        rs_viz_print_ansi(pts, n, a->width, a->height);
    }
}

static int op_cmeans(const rs_args_t *a)
{
    rs_point_t *pts = NULL;
    size_t n = 0, dim = 0;
    if (load_or_die(a, &pts, &n, &dim) != 0)
    {
        return 1;
    }
    if (a->op_arg1 == NULL)
    {
        fprintf(stderr, "cmeans: нужно число кластеров\n");
        free(pts);
        return 1;
    }
    size_t c = (size_t)atoi(a->op_arg1);
    if (c < 2)
    {
        fprintf(stderr, "cmeans: c должно быть >= 2\n");
        free(pts);
        return 1;
    }
    rs_cmeans_params_t p = rs_cmeans_default_params(c);
    rs_cmeans_result_t *res = rs_cmeans_run(pts, n, dim, &p);
    if (res == NULL)
    {
        fprintf(stderr, "cmeans: ошибка запуска\n");
        free(pts);
        return 1;
    }
    rs_cmeans_hard_assign(res, pts);
    printf("Fuzzy C-means: c=%zu, итераций=%zu, J=%.4f\n",
           c, res->iterations, res->objective);
    printf("Центроиды:\n");
    for (size_t j = 0; j < c; ++j)
    {
        printf("  c%zu = ", j);
        rs_point_print(&res->centroids[j]);
    }
    emit_results(a, pts, n);
    rs_cmeans_result_free(res);
    free(pts);
    return 0;
}

static int op_cmeans_elbow(const rs_args_t *a)
{
    rs_point_t *pts = NULL;
    size_t n = 0, dim = 0;
    if (load_or_die(a, &pts, &n, &dim) != 0)
    {
        return 1;
    }
    if (a->op_arg1 == NULL)
    {
        fprintf(stderr, "cmeans_elbow: нужно c_max\n");
        free(pts);
        return 1;
    }
    size_t cmax = (size_t)atoi(a->op_arg1);
    if (cmax < 3 || cmax > 32)
    {
        fprintf(stderr, "cmeans_elbow: c_max должно быть в [3,32]\n");
        free(pts);
        return 1;
    }
    double *J = (double *)calloc(cmax + 1, sizeof(double));
    size_t best_c = rs_cmeans_elbow(pts, n, dim, cmax, J);
    printf("Метод локтя:\n");
    printf("   c |       J(c)\n");
    printf("   --+------------\n");
    double Jmax = 0.0;
    for (size_t c = 2; c <= cmax; ++c)
    {
        if (J[c] > Jmax)
        {
            Jmax = J[c];
        }
    }
    for (size_t c = 2; c <= cmax; ++c)
    {
        
        int bar = (Jmax > 0) ? (int)(40.0 * J[c] / Jmax) : 0;
        printf("  %2zu | %10.4f ", c, J[c]);
        for (int i = 0; i < bar; ++i)
        {
            putchar('#');
        }
        if (c == best_c)
        {
            printf("   <-- локоть");
        }
        putchar('\n');
    }
    printf("Рекомендованное c = %zu\n", best_c);
    
    rs_cmeans_params_t p = rs_cmeans_default_params(best_c);
    rs_cmeans_result_t *res = rs_cmeans_run(pts, n, dim, &p);
    if (res != NULL)
    {
        rs_cmeans_hard_assign(res, pts);
        emit_results(a, pts, n);
        rs_cmeans_result_free(res);
    }
    free(J);
    free(pts);
    return 0;
}

static int op_dbscan(const rs_args_t *a)
{
    rs_point_t *pts = NULL;
    size_t n = 0, dim = 0;
    if (load_or_die(a, &pts, &n, &dim) != 0)
    {
        return 1;
    }
    if (a->op_arg1 == NULL)
    {
        fprintf(stderr, "dbscan: нужны параметры eps,minPts\n");
        free(pts);
        return 1;
    }
    double pp[2];
    size_t k = parse_coords(a->op_arg1, pp, 2);
    if (k < 2)
    {
        fprintf(stderr, "dbscan: ожидаю eps,minPts\n");
        free(pts);
        return 1;
    }
    rs_dbscan_params_t prm;
    prm.eps = pp[0];
    prm.min_pts = (size_t)pp[1];
    rs_dbscan_result_t *res = rs_dbscan_run(pts, n, dim, &prm);
    if (res == NULL)
    {
        fprintf(stderr, "dbscan: ошибка запуска\n");
        free(pts);
        return 1;
    }
    printf("DBSCAN eps=%.4f minPts=%zu: кластеров=%zu, шум=%zu\n",
           prm.eps, prm.min_pts, res->n_clusters, res->n_noise);
    emit_results(a, pts, n);
    rs_dbscan_result_free(res);
    free(pts);
    return 0;
}

int main(int argc, char **argv)
{
    rs_args_t a;
    if (parse_args(argc, argv, &a) != 0)
    {
        print_usage(argv[0]);
        return 1;
    }
    
    if (strcmp(a.op, "-kd_insert") == 0)
    {
        return op_kd_insert(&a);
    }
    if (strcmp(a.op, "-kd_nearest") == 0)
    {
        return op_kd_nearest(&a);
    }
    if (strcmp(a.op, "-cmeans") == 0)
    {
        return op_cmeans(&a);
    }
    if (strcmp(a.op, "-cmeans_elbow") == 0)
    {
        return op_cmeans_elbow(&a);
    }
    if (strcmp(a.op, "-dbscan") == 0)
    {
        return op_dbscan(&a);
    }
    fprintf(stderr, "Неизвестная операция: %s\n", a.op);
    print_usage(argv[0]);
    return 1;
}
