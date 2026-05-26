
#include "csv.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// пропуск ведущих пробельных символов
static char *csv_skip_ws(char *s)
{
    while (*s != '\0' && isspace((unsigned char)*s))
    {
        ++s;
    }
    return s;
}

// разбор одной строки CSV в массив double
static size_t csv_parse_line(char *line,
                             double *out,
                             size_t max_dim)
{
    size_t k = 0;       // число успешно распарсенных чисел
    char *p = line;     // текущая позиция в строке
    while (*p != '\0' && k < max_dim)
    {
        p = csv_skip_ws(p);
        if (*p == '\0' || *p == '\n' || *p == '\r')
        {
            break; // конец строки
        }
        char *endp = NULL;
        double v = strtod(p, &endp); // парсим число
        if (endp == p)
        {
            break; // не удалось распарсить — выходим
        }
        out[k++] = v;
        p = endp;
        p = csv_skip_ws(p);
        if (*p == ',')
        {
            ++p; // пропускаем разделитель
        }
    }
    return k;
}

// загрузка CSV-файла с координатами в массив точек
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
    size_t cap = 64;    // начальная ёмкость динамического массива точек
    size_t n = 0;       // текущее число прочитанных точек
    rs_point_t *arr = (rs_point_t *)malloc(sizeof(rs_point_t) * cap);
    if (arr == NULL)
    {
        fclose(f);
        return -1;
    }
    size_t dim = 0;         // размерность точек (определится по первой строке)
    char buf[1024];         // буфер под одну строку файла
    size_t line_no = 0;     // номер текущей строки для сообщений об ошибках
    while (fgets(buf, sizeof(buf), f) != NULL)
    {
        ++line_no;
        char *s = csv_skip_ws(buf);
        if (*s == '\0' || *s == '#' || *s == '\n' || *s == '\r')
        {
            continue; // пустая строка или комментарий
        }
        double tmp[RS_MAX_DIM]; // временный буфер под одну точку
        size_t got = csv_parse_line(s, tmp, RS_MAX_DIM);
        if (got < 2)
        {
            fprintf(stderr, "csv: строка %zu проигнорирована (полей: %zu)\n",
                    line_no, got);
            continue;
        }
        // эвристика: если последнее поле — целое в допустимом диапазоне,
        // считаем его cluster_id и убираем из размерности
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
            dim = effective; // фиксируем размерность по первой непустой строке
        }
        if (effective != dim)
        {
            fprintf(stderr, "csv: строка %zu: размерность %zu != %zu, пропуск\n",
                    line_no, effective, dim);
            continue;
        }
        // расширение массива при необходимости (удвоение ёмкости)
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
