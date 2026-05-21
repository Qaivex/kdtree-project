
#include "dbscan.h"
#include "kdtree.h"

#include <stdlib.h>
#include <string.h>

typedef struct
{
    rs_point_t *points;
    size_t n;
    char *visited;
    int *cluster_ids;
    rs_kdtree_t *tree;
    double eps;
    size_t min_pts;
} db_state_t;

typedef struct
{
    rs_point_t *data;
    size_t cap;
    size_t head;
    size_t tail;
} db_queue_t;

static int db_queue_init(db_queue_t *q, size_t cap)
{
    if (cap < 8)
    {
        cap = 8;
    }
    q->data = (rs_point_t *)malloc(sizeof(rs_point_t) * cap);
    if (q->data == NULL)
    {
        return -1;
    }
    q->cap = cap;
    q->head = 0;
    q->tail = 0;
    return 0;
}

static void db_queue_free(db_queue_t *q)
{
    free(q->data);
    q->data = NULL;
}

static int db_queue_push(db_queue_t *q, const rs_point_t *p)
{
    if (q->tail == q->cap)
    {
        size_t newcap = q->cap * 2;
        rs_point_t *nd = (rs_point_t *)realloc(q->data, sizeof(rs_point_t) * newcap);
        if (nd == NULL)
        {
            return -1;
        }
        q->data = nd;
        q->cap = newcap;
    }
    q->data[q->tail++] = *p;
    return 0;
}

static int db_queue_pop(db_queue_t *q, rs_point_t *out)
{
    if (q->head == q->tail)
    {
        return -1;
    }
    *out = q->data[q->head++];
    return 0;
}

/* поиск индекса точки в исходном массиве для DBSCAN */
static size_t db_locate(const db_state_t *st, const rs_point_t *p)
{
    if (p->index < st->n)
    {
        const rs_point_t *cand = &st->points[p->index];
        int match = 1;
        for (size_t d = 0; d < cand->dim; ++d)
        {
            double diff = cand->coord[d] - p->coord[d];
            if (diff < 0.0)
            {
                diff = -diff;
            }
            if (diff > 1e-12)
            {
                match = 0;
                break;
            }
        }
        if (match)
        {
            return p->index;
        }
    }
    for (size_t i = 0; i < st->n; ++i)
    {
        const rs_point_t *cand = &st->points[i];
        int match = 1;
        for (size_t d = 0; d < cand->dim; ++d)
        {
            double diff = cand->coord[d] - p->coord[d];
            if (diff < 0.0)
            {
                diff = -diff;
            }
            if (diff > 1e-12)
            {
                match = 0;
                break;
            }
        }
        if (match)
        {
            return i;
        }
    }
    return st->n;
}

/* ε-запрос окрестности DBSCAN через k-d дерево */
static rs_point_t *db_region_query(const db_state_t *st,
                                   const rs_point_t *p,
                                   size_t *out_count)
{
    size_t total = rs_kdtree_range(st->tree, p, st->eps, NULL, 0);
    if (total == 0)
    {
        *out_count = 0;
        return NULL;
    }
    rs_point_t *buf = (rs_point_t *)malloc(sizeof(rs_point_t) * total);
    if (buf == NULL)
    {
        *out_count = 0;
        return NULL;
    }
    size_t got = rs_kdtree_range(st->tree, p, st->eps, buf, total);
    *out_count = got;
    return buf;
}

/* расширение кластера DBSCAN через BFS-обход соседей */
static void db_expand(db_state_t *st,
                      size_t seed_idx,
                      rs_point_t *initial_neighbors,
                      size_t initial_count,
                      int cluster_id)
{
    st->cluster_ids[seed_idx] = cluster_id;
    st->points[seed_idx].cluster_id = cluster_id;
    db_queue_t q;
    if (db_queue_init(&q, initial_count + 8) != 0)
    {
        return;
    }
    for (size_t i = 0; i < initial_count; ++i)
    {
        db_queue_push(&q, &initial_neighbors[i]);
    }
    rs_point_t cur;
    while (db_queue_pop(&q, &cur) == 0)
    {
        size_t idx = db_locate(st, &cur);
        if (idx >= st->n)
        {
            continue;
        }
        
        if (st->cluster_ids[idx] == RS_DBSCAN_NOISE)
        {
            st->cluster_ids[idx] = cluster_id;
            st->points[idx].cluster_id = cluster_id;
        }
        if (st->visited[idx])
        {
            continue;
        }
        st->visited[idx] = 1;
        st->cluster_ids[idx] = cluster_id;
        st->points[idx].cluster_id = cluster_id;
        size_t cnt = 0;
        rs_point_t *nbrs = db_region_query(st, &cur, &cnt);
        if (cnt >= st->min_pts)
        {
            for (size_t i = 0; i < cnt; ++i)
            {
                db_queue_push(&q, &nbrs[i]);
            }
        }
        free(nbrs);
    }
    db_queue_free(&q);
}

rs_dbscan_result_t *rs_dbscan_run(rs_point_t *points,
                                  size_t n,
                                  size_t dim,
                                  const rs_dbscan_params_t *params)
{
    if (points == NULL || n == 0 || params == NULL
        || params->eps <= 0.0 || params->min_pts == 0)
    {
        return NULL;
    }
    
    for (size_t i = 0; i < n; ++i)
    {
        points[i].index = i;
        points[i].cluster_id = RS_DBSCAN_NOISE;
    }
    rs_dbscan_result_t *res = (rs_dbscan_result_t *)calloc(1, sizeof(*res));
    if (res == NULL)
    {
        return NULL;
    }
    res->n = n;
    res->cluster_ids = (int *)malloc(sizeof(int) * n);
    if (res->cluster_ids == NULL)
    {
        free(res);
        return NULL;
    }
    for (size_t i = 0; i < n; ++i)
    {
        res->cluster_ids[i] = RS_DBSCAN_NOISE;
    }
    db_state_t st;
    st.points = points;
    st.n = n;
    st.visited = (char *)calloc(n, 1);
    st.cluster_ids = res->cluster_ids;
    st.tree = rs_kdtree_build(points, n, dim);
    st.eps = params->eps;
    st.min_pts = params->min_pts;
    if (st.visited == NULL || st.tree == NULL)
    {
        free(st.visited);
        rs_kdtree_destroy(st.tree);
        free(res->cluster_ids);
        free(res);
        return NULL;
    }
    int next_cluster = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (st.visited[i])
        {
            continue;
        }
        st.visited[i] = 1;
        size_t cnt = 0;
        rs_point_t *nbrs = db_region_query(&st, &points[i], &cnt);
        if (cnt < st.min_pts)
        {
            free(nbrs);
            continue;
        }
        db_expand(&st, i, nbrs, cnt, next_cluster);
        free(nbrs);
        ++next_cluster;
    }
    res->n_clusters = (size_t)next_cluster;
    res->n_noise = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (res->cluster_ids[i] == RS_DBSCAN_NOISE)
        {
            res->n_noise += 1;
        }
    }
    free(st.visited);
    rs_kdtree_destroy(st.tree);
    return res;
}

void rs_dbscan_result_free(rs_dbscan_result_t *res)
{
    if (res == NULL)
    {
        return;
    }
    free(res->cluster_ids);
    free(res);
}
