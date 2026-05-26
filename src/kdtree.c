
#include "kdtree.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// сравнение двух точек по координатам с допуском 1e-12
static int kd_points_equal(const rs_point_t *a, const rs_point_t *b)
{
    for (size_t i = 0; i < a->dim; ++i)
    {
        double diff = a->coord[i] - b->coord[i]; // разность по оси
        if (diff < 0.0)
        {
            diff = -diff; // модуль
        }
        if (diff > 1e-12)
        {
            return 0; // расхождение больше допуска — точки разные
        }
    }
    return 1;
}

// выделение нового узла с копией точки и заданной осью разбиения
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

// рекурсивное освобождение поддерева в post-order
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

// глобальная ось для kd_cmp_axis (qsort не принимает контекст)
static int g_kd_sort_axis = 0;

// компаратор qsort по оси g_kd_sort_axis
static int kd_cmp_axis(const void *a, const void *b)
{
    const rs_point_t *pa = (const rs_point_t *)a;
    const rs_point_t *pb = (const rs_point_t *)b;
    double da = pa->coord[g_kd_sort_axis]; // координата первой точки
    double db = pb->coord[g_kd_sort_axis]; // координата второй точки
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

// рекурсивное построение сбалансированного k-d дерева через медиану
static rs_kdnode_t *kd_build_recursive(rs_point_t *pts,
                                       size_t n,
                                       size_t dim,
                                       int depth)
{
    if (n == 0)
    {
        return NULL;
    }
    int axis = depth % (int)dim; // ось разбиения на этом уровне
    g_kd_sort_axis = axis;
    qsort(pts, n, sizeof(*pts), kd_cmp_axis); // сортируем по оси
    size_t mid = n / 2; // индекс медианы
    rs_kdnode_t *node = kd_node_new(&pts[mid], axis);
    if (node == NULL)
    {
        return NULL;
    }
    // левая половина (до медианы) — в левое поддерево
    node->left  = kd_build_recursive(pts, mid, dim, depth + 1);
    // правая половина (после медианы) — в правое поддерево
    node->right = kd_build_recursive(pts + mid + 1, n - mid - 1, dim, depth + 1);
    return node;
}

// создание пустого k-d дерева заданной размерности
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

// освобождение дерева вместе со всеми узлами
void rs_kdtree_destroy(rs_kdtree_t *tree)
{
    if (tree == NULL)
    {
        return;
    }
    kd_node_free(tree->root);
    free(tree);
}

// построение сбалансированного дерева из массива точек
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
    // копия входных данных, чтобы свободно их сортировать
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

// рекурсивная вставка точки в k-d дерево
static rs_kdnode_t *kd_insert_rec(rs_kdnode_t *node,
                                  const rs_point_t *p,
                                  size_t dim,
                                  int depth,
                                  int *ok)
{
    if (node == NULL)
    {
        // дошли до пустого места — создаём новый узел
        int axis = depth % (int)dim;
        rs_kdnode_t *fresh = kd_node_new(p, axis);
        if (fresh == NULL)
        {
            *ok = 0; // сигнализируем о неудаче выделения памяти
        }
        return fresh;
    }
    // дубликат игнорируем, чтобы не плодить идентичные точки
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

// публичная вставка: проверяет аргументы и вызывает рекурсивную версию
bool rs_kdtree_insert(rs_kdtree_t *tree, const rs_point_t *p)
{
    if (tree == NULL || p == NULL || p->dim != tree->dim)
    {
        return false;
    }
    int ok = 1; // флаг успешной вставки
    size_t before = (tree->root != NULL) ? 1 : 0;
    (void)before;
    tree->root = kd_insert_rec(tree->root, p, tree->dim, 0, &ok);
    if (ok)
    {
        tree->size += 1;
    }
    return ok ? true : false;
}

// поиск узла с минимальной координатой по оси для удаления Bentley
static rs_kdnode_t *kd_find_min(rs_kdnode_t *node, int axis, size_t dim)
{
    if (node == NULL)
    {
        return NULL;
    }
    // если в этом узле режут по нужной оси — минимум только слева
    if (node->axis == axis)
    {
        if (node->left == NULL)
        {
            return node;
        }
        return kd_find_min(node->left, axis, dim);
    }
    // иначе ищем минимум во всём поддереве
    rs_kdnode_t *lmin = kd_find_min(node->left,  axis, dim);
    rs_kdnode_t *rmin = kd_find_min(node->right, axis, dim);
    rs_kdnode_t *best = node; // текущий лучший кандидат
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

// рекурсивное удаление узла k-d дерева по схеме Bentley
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
        // случай 1: есть правое поддерево — берём минимум оттуда
        if (node->right != NULL)
        {
            rs_kdnode_t *minNode = kd_find_min(node->right, axis, dim);
            rs_point_copy(&node->point, &minNode->point);
            node->right = kd_remove_rec(node->right, &minNode->point, dim, removed);
            *removed = 1;
            return node;
        }
        // случай 2: только левое — переносим его направо и берём минимум
        if (node->left != NULL)
        {
            rs_kdnode_t *minNode = kd_find_min(node->left, axis, dim);
            rs_point_copy(&node->point, &minNode->point);
            node->right = kd_remove_rec(node->left, &minNode->point, dim, removed);
            node->left = NULL;
            *removed = 1;
            return node;
        }
        // случай 3: лист — просто освобождаем
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

// публичное удаление: проверяет аргументы и вызывает рекурсивную версию
bool rs_kdtree_remove(rs_kdtree_t *tree, const rs_point_t *p)
{
    if (tree == NULL || p == NULL || p->dim != tree->dim)
    {
        return false;
    }
    int removed = 0; // флаг — было ли реально удаление
    tree->root = kd_remove_rec(tree->root, p, tree->dim, &removed);
    if (removed)
    {
        tree->size -= 1;
        return true;
    }
    return false;
}

// состояние поиска ближайшего соседа
typedef struct
{
    const rs_point_t *query;   // запрашиваемая точка
    const rs_kdnode_t *best;   // текущий лучший узел
    double best_sq;            // квадрат расстояния до лучшего узла
} kd_nn_state_t;

// рекурсивный поиск ближайшего соседа с отсечением по гиперплоскости
static void kd_nn_rec(const rs_kdnode_t *node, kd_nn_state_t *st)
{
    if (node == NULL)
    {
        return;
    }
    double d = rs_point_sqdist(&node->point, st->query); // расстояние до текущего узла
    if (d < st->best_sq)
    {
        st->best_sq = d;
        st->best = node;
    }
    int axis = node->axis;
    double diff = st->query->coord[axis] - node->point.coord[axis]; // знак выбирает сторону
    const rs_kdnode_t *near_side = (diff < 0.0) ? node->left  : node->right; // ближняя половина
    const rs_kdnode_t *far_side  = (diff < 0.0) ? node->right : node->left;  // дальняя половина
    kd_nn_rec(near_side, st);
    // в дальнюю сторону идём только если плоскость ближе текущего лучшего
    if (diff * diff < st->best_sq)
    {
        kd_nn_rec(far_side, st);
    }
}

// публичный поиск ближайшего соседа: подготавливает состояние и стартует обход
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
    st.best_sq = DBL_MAX; // стартовое «бесконечно далеко»
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

// состояние range-запроса
typedef struct
{
    const rs_point_t *query; // запрашиваемая точка
    double r_sq;             // квадрат радиуса (чтобы не делать sqrt)
    rs_point_t *out;         // буфер для результатов
    size_t cap;              // размер буфера
    size_t count;            // сколько точек найдено всего (может быть > cap)
} kd_range_state_t;

// рекурсивный range-запрос k-d дерева: все точки в радиусе
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
            st->out[st->count] = node->point; // кладём в буфер, если есть место
        }
        st->count += 1;
    }
    int axis = node->axis;
    double diff = st->query->coord[axis] - node->point.coord[axis];
    // спускаемся в ближнюю половину, дальнюю — только если шар пересекает плоскость
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

// публичный range-запрос: подготавливает состояние и стартует обход
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
    st.r_sq = radius * radius; // храним квадрат, сравниваем без sqrt
    st.out = out;
    st.cap = cap;
    st.count = 0;
    kd_range_rec(tree->root, &st);
    return st.count;
}

// текущее число точек в дереве
size_t rs_kdtree_size(const rs_kdtree_t *tree)
{
    return (tree != NULL) ? tree->size : 0;
}
