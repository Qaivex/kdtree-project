
#include "kdtree.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* сравнение двух точек по координатам с допуском */
static int kd_points_equal(const rs_point_t *a, const rs_point_t *b)
{
    for (size_t i = 0; i < a->dim; ++i)
    {
        double diff = a->coord[i] - b->coord[i];
        if (diff < 0.0)
        {
            diff = -diff;
        }
        if (diff > 1e-12)
        {
            return 0;
        }
    }
    return 1;
}

static rs_kdnode_t *kd_node_new(const rs_point_t *p, int axis)
{
    rs_kdnode_t *n = (rs_kdnode_t *)malloc(sizeof(*n));
    if (n == NULL)
    {
        return NULL;
    }
    rs_point_copy(&n->point, p);
    n->axis = axis;
    n->left = NULL;
    n->right = NULL;
    return n;
}

static void kd_node_free(rs_kdnode_t *node)
{
    if (node == NULL)
    {
        return;
    }
    kd_node_free(node->left);
    kd_node_free(node->right);
    free(node);
}

static int g_kd_sort_axis = 0;

/* компаратор qsort по оси g_kd_sort_axis */
static int kd_cmp_axis(const void *a, const void *b)
{
    const rs_point_t *pa = (const rs_point_t *)a;
    const rs_point_t *pb = (const rs_point_t *)b;
    double da = pa->coord[g_kd_sort_axis];
    double db = pb->coord[g_kd_sort_axis];
    if (da < db)
    {
        return -1;
    }
    if (da > db)
    {
        return 1;
    }
    return 0;
}

/* рекурсивное построение сбалансированного k-d дерева через медиану */
static rs_kdnode_t *kd_build_recursive(rs_point_t *pts,
                                       size_t n,
                                       size_t dim,
                                       int depth)
{
    if (n == 0)
    {
        return NULL;
    }
    int axis = depth % (int)dim;
    
    g_kd_sort_axis = axis;
    qsort(pts, n, sizeof(*pts), kd_cmp_axis);
    size_t mid = n / 2;
    rs_kdnode_t *node = kd_node_new(&pts[mid], axis);
    if (node == NULL)
    {
        return NULL;
    }
    
    node->left  = kd_build_recursive(pts, mid, dim, depth + 1);
    node->right = kd_build_recursive(pts + mid + 1, n - mid - 1, dim, depth + 1);
    return node;
}

rs_kdtree_t *rs_kdtree_create(size_t dim)
{
    if (dim == 0 || dim > RS_MAX_DIM)
    {
        return NULL;
    }
    rs_kdtree_t *t = (rs_kdtree_t *)malloc(sizeof(*t));
    if (t == NULL)
    {
        return NULL;
    }
    t->root = NULL;
    t->dim = dim;
    t->size = 0;
    return t;
}

void rs_kdtree_destroy(rs_kdtree_t *tree)
{
    if (tree == NULL)
    {
        return;
    }
    kd_node_free(tree->root);
    free(tree);
}

rs_kdtree_t *rs_kdtree_build(const rs_point_t *points, size_t n, size_t dim)
{
    rs_kdtree_t *t = rs_kdtree_create(dim);
    if (t == NULL)
    {
        return NULL;
    }
    if (n == 0)
    {
        return t;
    }
    
    rs_point_t *buf = (rs_point_t *)malloc(sizeof(*buf) * n);
    if (buf == NULL)
    {
        rs_kdtree_destroy(t);
        return NULL;
    }
    memcpy(buf, points, sizeof(*buf) * n);
    t->root = kd_build_recursive(buf, n, dim, 0);
    t->size = n;
    free(buf);
    return t;
}

/* рекурсивная вставка точки в k-d дерево */
static rs_kdnode_t *kd_insert_rec(rs_kdnode_t *node,
                                  const rs_point_t *p,
                                  size_t dim,
                                  int depth,
                                  int *ok)
{
    if (node == NULL)
    {
        int axis = depth % (int)dim;
        rs_kdnode_t *fresh = kd_node_new(p, axis);
        if (fresh == NULL)
        {
            *ok = 0;
        }
        return fresh;
    }
    
    if (kd_points_equal(&node->point, p))
    {
        return node;
    }
    int axis = node->axis;
    if (p->coord[axis] < node->point.coord[axis])
    {
        node->left = kd_insert_rec(node->left, p, dim, depth + 1, ok);
    }
    else
    {
        node->right = kd_insert_rec(node->right, p, dim, depth + 1, ok);
    }
    return node;
}

bool rs_kdtree_insert(rs_kdtree_t *tree, const rs_point_t *p)
{
    if (tree == NULL || p == NULL || p->dim != tree->dim)
    {
        return false;
    }
    int ok = 1;
    size_t before = (tree->root != NULL) ? 1 : 0;
    (void)before;
    tree->root = kd_insert_rec(tree->root, p, tree->dim, 0, &ok);
    if (ok)
    {
        tree->size += 1;
    }
    return ok ? true : false;
}

/* поиск узла с минимальной координатой по оси для удаления Bentley */
static rs_kdnode_t *kd_find_min(rs_kdnode_t *node, int axis, size_t dim)
{
    if (node == NULL)
    {
        return NULL;
    }
    
    if (node->axis == axis)
    {
        if (node->left == NULL)
        {
            return node;
        }
        return kd_find_min(node->left, axis, dim);
    }
    
    rs_kdnode_t *lmin = kd_find_min(node->left,  axis, dim);
    rs_kdnode_t *rmin = kd_find_min(node->right, axis, dim);
    rs_kdnode_t *best = node;
    if (lmin != NULL && lmin->point.coord[axis] < best->point.coord[axis])
    {
        best = lmin;
    }
    if (rmin != NULL && rmin->point.coord[axis] < best->point.coord[axis])
    {
        best = rmin;
    }
    return best;
}

/* рекурсивное удаление узла k-d дерева по схеме Bentley */
static rs_kdnode_t *kd_remove_rec(rs_kdnode_t *node,
                                  const rs_point_t *p,
                                  size_t dim,
                                  int *removed)
{
    if (node == NULL)
    {
        return NULL;
    }
    int axis = node->axis;
    if (kd_points_equal(&node->point, p))
    {
        
        if (node->right != NULL)
        {
            
            rs_kdnode_t *minNode = kd_find_min(node->right, axis, dim);
            rs_point_copy(&node->point, &minNode->point);
            node->right = kd_remove_rec(node->right, &minNode->point, dim, removed);
            
            *removed = 1;
            return node;
        }
        if (node->left != NULL)
        {
            
            rs_kdnode_t *minNode = kd_find_min(node->left, axis, dim);
            rs_point_copy(&node->point, &minNode->point);
            node->right = kd_remove_rec(node->left, &minNode->point, dim, removed);
            node->left = NULL;
            *removed = 1;
            return node;
        }
        
        free(node);
        *removed = 1;
        return NULL;
    }
    if (p->coord[axis] < node->point.coord[axis])
    {
        node->left = kd_remove_rec(node->left, p, dim, removed);
    }
    else
    {
        node->right = kd_remove_rec(node->right, p, dim, removed);
    }
    return node;
}

bool rs_kdtree_remove(rs_kdtree_t *tree, const rs_point_t *p)
{
    if (tree == NULL || p == NULL || p->dim != tree->dim)
    {
        return false;
    }
    int removed = 0;
    tree->root = kd_remove_rec(tree->root, p, tree->dim, &removed);
    if (removed)
    {
        tree->size -= 1;
        return true;
    }
    return false;
}

typedef struct
{
    const rs_point_t *query;
    const rs_kdnode_t *best;
    double best_sq;
} kd_nn_state_t;

/* рекурсивный поиск ближайшего соседа с отсечением по гиперплоскости */
static void kd_nn_rec(const rs_kdnode_t *node, kd_nn_state_t *st)
{
    if (node == NULL)
    {
        return;
    }
    double d = rs_point_sqdist(&node->point, st->query);
    if (d < st->best_sq)
    {
        st->best_sq = d;
        st->best = node;
    }
    int axis = node->axis;
    double diff = st->query->coord[axis] - node->point.coord[axis];
    const rs_kdnode_t *near_side = (diff < 0.0) ? node->left  : node->right;
    const rs_kdnode_t *far_side  = (diff < 0.0) ? node->right : node->left;
    kd_nn_rec(near_side, st);
    
    if (diff * diff < st->best_sq)
    {
        kd_nn_rec(far_side, st);
    }
}

bool rs_kdtree_nearest(const rs_kdtree_t *tree,
                       const rs_point_t *query,
                       rs_point_t *out,
                       double *out_sqdist)
{
    if (tree == NULL || tree->root == NULL || query == NULL || out == NULL)
    {
        return false;
    }
    kd_nn_state_t st;
    st.query = query;
    st.best = NULL;
    st.best_sq = DBL_MAX;
    kd_nn_rec(tree->root, &st);
    if (st.best == NULL)
    {
        return false;
    }
    rs_point_copy(out, &st.best->point);
    if (out_sqdist != NULL)
    {
        *out_sqdist = st.best_sq;
    }
    return true;
}

typedef struct
{
    const rs_point_t *query;
    double r_sq;
    rs_point_t *out;
    size_t cap;
    size_t count;
} kd_range_state_t;

/* рекурсивный range-запрос k-d дерева: все точки в радиусе */
static void kd_range_rec(const rs_kdnode_t *node, kd_range_state_t *st)
{
    if (node == NULL)
    {
        return;
    }
    double d = rs_point_sqdist(&node->point, st->query);
    if (d <= st->r_sq)
    {
        if (st->count < st->cap)
        {
            st->out[st->count] = node->point;
        }
        st->count += 1;
    }
    int axis = node->axis;
    double diff = st->query->coord[axis] - node->point.coord[axis];
    
    if (diff < 0.0)
    {
        kd_range_rec(node->left, st);
        if (diff * diff <= st->r_sq)
        {
            kd_range_rec(node->right, st);
        }
    }
    else
    {
        kd_range_rec(node->right, st);
        if (diff * diff <= st->r_sq)
        {
            kd_range_rec(node->left, st);
        }
    }
}

size_t rs_kdtree_range(const rs_kdtree_t *tree,
                       const rs_point_t *query,
                       double radius,
                       rs_point_t *out,
                       size_t cap)
{
    if (tree == NULL || query == NULL)
    {
        return 0;
    }
    kd_range_state_t st;
    st.query = query;
    st.r_sq = radius * radius;
    st.out = out;
    st.cap = cap;
    st.count = 0;
    kd_range_rec(tree->root, &st);
    return st.count;
}

size_t rs_kdtree_size(const rs_kdtree_t *tree)
{
    return (tree != NULL) ? tree->size : 0;
}
