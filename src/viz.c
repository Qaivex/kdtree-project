
#include "viz.h"

#include <stdio.h>
#include <stdlib.h>

static const rs_rgb_t g_palette[] = {
    {  31, 119, 180 },
    { 255, 127,  14 },
    {  44, 160,  44 },
    { 214,  39,  40 },
    { 148, 103, 189 },
    { 140,  86,  75 },
    { 227, 119, 194 },
    { 127, 127, 127 },
    { 188, 189,  34 },
    {  23, 190, 207 },
    {   0, 200, 255 },
    { 255, 215,   0 },
    { 200,   0, 200 },
    {   0, 150,  90 },
    { 255,  90,  90 },
    { 100,  60, 200 }
};
static const size_t g_palette_n = sizeof(g_palette) / sizeof(g_palette[0]);

rs_rgb_t rs_viz_color_for_cluster(int cluster_id)
{
    if (cluster_id < 0)
    {
        rs_rgb_t g = { 128, 128, 128 };
        return g;
    }
    return g_palette[(size_t)cluster_id % g_palette_n];
}

/* вычисление ограничивающего прямоугольника облака точек */
static void viz_bbox(const rs_point_t *pts,
                     size_t n,
                     double *xmin, double *ymin,
                     double *xmax, double *ymax)
{
    *xmin = *ymin = 1e300;
    *xmax = *ymax = -1e300;
    for (size_t i = 0; i < n; ++i)
    {
        double x = pts[i].coord[0];
        double y = (pts[i].dim >= 2) ? pts[i].coord[1] : 0.0;
        if (x < *xmin)
        {
            *xmin = x;
        }
        if (x > *xmax)
        {
            *xmax = x;
        }
        if (y < *ymin)
        {
            *ymin = y;
        }
        if (y > *ymax)
        {
            *ymax = y;
        }
    }
    double padx = (*xmax - *xmin) * 0.05 + 1e-9;
    double pady = (*ymax - *ymin) * 0.05 + 1e-9;
    *xmin -= padx;
    *ymin -= pady;
    *xmax += padx;
    *ymax += pady;
}

void rs_viz_print_ansi(const rs_point_t *points,
                       size_t n,
                       int width,
                       int height)
{
    if (points == NULL || n == 0 || width <= 0 || height <= 0)
    {
        return;
    }
    double xmin, ymin, xmax, ymax;
    viz_bbox(points, n, &xmin, &ymin, &xmax, &ymax);
    int *grid = (int *)malloc(sizeof(int) * (size_t)width * (size_t)height);
    if (grid == NULL)
    {
        return;
    }
    for (int i = 0; i < width * height; ++i)
    {
        grid[i] = -2;
    }
    for (size_t i = 0; i < n; ++i)
    {
        double x = points[i].coord[0];
        double y = (points[i].dim >= 2) ? points[i].coord[1] : 0.0;
        int cx = (int)((x - xmin) / (xmax - xmin) * (width - 1));
        int cy = (int)((ymax - y) / (ymax - ymin) * (height - 1));
        if (cx < 0)
        {
            cx = 0;
        }
        if (cx >= width)
        {
            cx = width - 1;
        }
        if (cy < 0)
        {
            cy = 0;
        }
        if (cy >= height)
        {
            cy = height - 1;
        }
        int *cell = &grid[cy * width + cx];
        if (*cell == -2 || (*cell == -1 && points[i].cluster_id >= 0))
        {
            *cell = points[i].cluster_id;
        }
    }
    for (int row = 0; row < height; ++row)
    {
        for (int col = 0; col < width; ++col)
        {
            int cid = grid[row * width + col];
            if (cid == -2)
            {
                fputc(' ', stdout);
            }
            else
            {
                rs_rgb_t c = rs_viz_color_for_cluster(cid);
                printf("\x1b[38;2;%u;%u;%um*\x1b[0m",
                       (unsigned)c.r, (unsigned)c.g, (unsigned)c.b);
            }
        }
        fputc('\n', stdout);
    }
    int max_cid = -1;
    for (size_t i = 0; i < n; ++i)
    {
        if (points[i].cluster_id > max_cid)
        {
            max_cid = points[i].cluster_id;
        }
    }
    printf("\nЛегенда:\n");
    for (int c = -1; c <= max_cid; ++c)
    {
        rs_rgb_t col = rs_viz_color_for_cluster(c);
        printf("  \x1b[38;2;%u;%u;%um###\x1b[0m  %s\n",
               (unsigned)col.r, (unsigned)col.g, (unsigned)col.b,
               (c == -1) ? "шум" : "кластер");
    }
    free(grid);
}
