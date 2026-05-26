
#ifndef RS_KDTREE_H
#define RS_KDTREE_H

#include "point.h"
#include <stdbool.h>

// узел k-d дерева
typedef struct rs_kdnode
{
    rs_point_t point;        // точка, хранимая в узле
    int axis;                // ось разбиения на данной глубине
    struct rs_kdnode *left;  // левое поддерево (координата < медианы)
    struct rs_kdnode *right; // правое поддерево (координата >= медианы)
} rs_kdnode_t;

// корневая структура k-d дерева
typedef struct
{
    rs_kdnode_t *root; // корневой узел
    size_t dim;        // размерность пространства
    size_t size;       // число точек в дереве
} rs_kdtree_t;

// создание пустого k-d дерева заданной размерности
rs_kdtree_t *rs_kdtree_create(size_t dim);

// освобождение всех ресурсов дерева
void rs_kdtree_destroy(rs_kdtree_t *tree);

// сбалансированное построение k-d дерева из массива точек
rs_kdtree_t *rs_kdtree_build(const rs_point_t *points, size_t n, size_t dim);

// вставка одной точки в k-d дерево
bool rs_kdtree_insert(rs_kdtree_t *tree, const rs_point_t *p);

// удаление точки из k-d дерева по координатам
bool rs_kdtree_remove(rs_kdtree_t *tree, const rs_point_t *p);

// поиск ближайшего соседа к query в k-d дереве
bool rs_kdtree_nearest(const rs_kdtree_t *tree,
                       const rs_point_t *query,
                       rs_point_t *out,
                       double *out_sqdist);

// поиск всех точек в радиусе radius от query
size_t rs_kdtree_range(const rs_kdtree_t *tree,
                       const rs_point_t *query,
                       double radius,
                       rs_point_t *out,
                       size_t cap);

// количество точек, хранимых в дереве
size_t rs_kdtree_size(const rs_kdtree_t *tree);

#endif
