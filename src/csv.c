
#include "csv.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char *csv_skip_ws(char *s)
{
    while (*s != '\0' && isspace((unsigned char)*s))
    {
        ++s;
    }
    return s;
}

/* разбор одной строки CSV в массив double */
static size_t csv_parse_line(char *line,
                             double *out,
                             size_t max_dim)
{
    size_t k = 0;
    char *p = line;
    while (*p != '\0' && k < max_dim)
    {
        p = csv_skip_ws(p);
        if (*p == '\0' || *p == '\n' || *p == '\r')
        {
            break;
        }
        char *endp = NULL;
        double v = strtod(p, &endp);
        if (endp == p)
        {
            
            break;
        }
        out[k++] = v;
        p = endp;
        p = csv_skip_ws(p);
        if (*p == ',')
        {
            ++p;
        }
    }
    return k;
}

int rs_csv_load(const char *path,
                rs_point_t **out_pts,
                size_t *out_n,
                size_t *out_dim)
{
    if (path == NULL || out_pts == NULL || out_n == NULL || out_dim == NULL)
    {
        return -1;
    }
    FILE *f = fopen(path, "r");
    if (f == NULL)
    {
        fprintf(stderr, "csv: не удалось открыть '%s'\n", path);
        return -1;
    }
    
    size_t cap = 64;
    size_t n = 0;
    rs_point_t *arr = (rs_point_t *)malloc(sizeof(rs_point_t) * cap);
    if (arr == NULL)
    {
        fclose(f);
        return -1;
    }
    size_t dim = 0;
    char buf[1024];
    size_t line_no = 0;
    while (fgets(buf, sizeof(buf), f) != NULL)
    {
        ++line_no;
        char *s = csv_skip_ws(buf);
        if (*s == '\0' || *s == '#' || *s == '\n' || *s == '\r')
        {
            continue;
        }
        double tmp[RS_MAX_DIM];
        size_t got = csv_parse_line(s, tmp, RS_MAX_DIM);
        if (got < 2)
        {
            fprintf(stderr, "csv: строка %zu проигнорирована (полей: %zu)\n",
                    line_no, got);
            continue;
        }
        
        size_t effective = got;
        if (got >= 3)
        {
            double last = tmp[got - 1];
            if (last == (double)(long long)last && last >= -1.0 && last < 1024.0)
            {
                
                if (dim != 0 && got == dim + 1)
                {
                    effective = got - 1;
                }
            }
        }
        if (dim == 0)
        {
            dim = effective;
        }
        if (effective != dim)
        {
            fprintf(stderr, "csv: строка %zu: размерность %zu != %zu, пропуск\n",
                    line_no, effective, dim);
            continue;
        }
        
        if (n == cap)
        {
            size_t newcap = cap * 2;
            rs_point_t *nb = (rs_point_t *)realloc(arr, sizeof(rs_point_t) * newcap);
            if (nb == NULL)
            {
                free(arr);
                fclose(f);
                return -1;
            }
            arr = nb;
            cap = newcap;
        }
        rs_point_t *p = &arr[n];
        memset(p, 0, sizeof(*p));
        for (size_t d = 0; d < dim; ++d)
        {
            p->coord[d] = tmp[d];
        }
        p->dim = dim;
        p->cluster_id = -1;
        p->index = n;
        ++n;
    }
    fclose(f);
    if (n == 0)
    {
        free(arr);
        fprintf(stderr, "csv: файл '%s' пуст или невалиден\n", path);
        return -1;
    }
    *out_pts = arr;
    *out_n = n;
    *out_dim = dim;
    return 0;
}
