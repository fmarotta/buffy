#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "data_structures.h"

/* Dynamic array of fixed-length strings */

void InitStringsDynArray(StringsDynArray *sda)
{
    sda->size = STRINGSDYNARRAY_INITIAL_SIZE;
    sda->a = malloc(sizeof(char *) * sda->size);
    sda->used = 0;
}

void AddStringToDynArray(StringsDynArray *sda, const char *s)
{
    if (sda->used == sda->size) {
        sda->size *= 2;
        sda->a = realloc(sda->a, sizeof(char *) * sda->size);
    }
    sda->a[sda->used++] = malloc(sizeof(char) * (STRINGSDYNARRAY_MAX_LENGTH + 1));
    strncpy(sda->a[sda->used - 1], s, STRINGSDYNARRAY_MAX_LENGTH);
    sda->a[sda->used - 1][STRINGSDYNARRAY_MAX_LENGTH] = '\0';
}

int GetSizeOfStringsDynArray(StringsDynArray *sda)
{
    return sda->used;
}

char *GetNthString(StringsDynArray *sda, int n)
{
    return sda->a[n];
}

void FreeStringsDynArray(StringsDynArray *sda)
{
    for (int i = 0; i < sda->used; i++)
        free(sda->a[i]);
    free(sda->a);
    sda->a = NULL;
    sda->size = sda->used = 0;
}


/* Dynamic array of characters */

void InitCharsDynArray(CharsDynArray *sda)
{
    sda->size = CHARSDYNARRAY_INITIAL_SIZE;
    sda->a = malloc(sizeof(char) * sda->size);
    sda->used = 0;
}

void AddCharToDynArray(CharsDynArray *sda, const char c)
{
    if (sda->used == sda->size) {
        sda->size *= 2;
        sda->a = realloc(sda->a, sizeof(char) * sda->size);
    }
    sda->a[sda->used++] = c;
}

char GetNthChar(CharsDynArray *sda, int n)
{
    return sda->a[n];
}

void FreeCharsDynArray(CharsDynArray *sda)
{
    free(sda->a);
    sda->a = NULL;
    sda->size = sda->used = 0;
}


/* Dynamic array of ints */

void InitIntsDynArray(IntsDynArray *sda)
{
    sda->size = INTSDYNARRAY_INITIAL_SIZE;
    sda->a = malloc(sizeof(int) * sda->size);
    sda->used = 0;
}

void AddIntToDynArray(IntsDynArray *sda, const int i)
{
    if (sda->used == sda->size) {
        sda->size *= 2;
        sda->a = realloc(sda->a, sizeof(int) * sda->size);
    }
    sda->a[sda->used++] = i;
}

int GetNthInt(IntsDynArray *sda, int n)
{
    return sda->a[n];
}

int GetSizeOfIntsDynArray(IntsDynArray *sda)
{
    return sda->used;
}

void FreeIntsDynArray(IntsDynArray *sda)
{
    free(sda->a);
    sda->a = NULL;
    sda->size = sda->used = 0;
}


/* Dynamic array of bools */

void InitBoolsDynArray(BoolsDynArray *sda)
{
    sda->size = BOOLSDYNARRAY_INITIAL_SIZE;
    sda->a = malloc(sizeof(bool) * sda->size);
    sda->used = 0;
}

void AddBoolToDynArray(BoolsDynArray *sda, const bool b)
{
    if (sda->used == sda->size) {
        sda->size *= 2;
        sda->a = realloc(sda->a, sizeof(int) * sda->size);
    }
    sda->a[sda->used++] = b;
}

bool GetNthBool(BoolsDynArray *sda, int n)
{
    return sda->a[n];
}

unsigned long long int GetBoolsDynArrayAsInt(BoolsDynArray *sda)
{
    unsigned long long int decimal = 0;
    for (int i = sda->used - 1; i >= 0; i--)
        decimal += sda->a[i] * pow(2, sda->used - 1 - i);
    return decimal;
}

void FreeBoolsDynArray(BoolsDynArray *sda)
{
    free(sda->a);
    sda->a = NULL;
    sda->size = sda->used = 0;
}


/* Hash */

static int hash(Hash *h, unsigned long long int index)
{
    return index % h->size;
}

void InitHash(Hash *h, int n)
{
    h->map = malloc(sizeof(struct HashNode *) * n);
    for (int i = 0; i < n; i++)
        h->map[i] = NULL;
    h->size = n;
    h->used = 0;
}

void AddToHash(Hash *h, unsigned long long int index, double e)
{
    int i = hash(h, index);
    while (h->map[i] != NULL)
    {
        i++;
        i %= h->size;
    }
    h->map[i] = malloc(sizeof(struct HashNode));
    h->map[i]->index = index;
    h->map[i]->payload = e;
    h->used++;
}

bool AlreadyInHash(Hash *h, unsigned long long int index)
{
    int i = hash(h, index);
    while (h->map[i] != NULL)
    {
        if (h->map[i]->index == index)
            return true;
        i++;
        i %= h->size;
    }
    return false;
}

double GetHash(Hash *h, unsigned long long int index)
{
    int i = hash(h, index);
    while (h->map[i]->index != index)
    {
        i++;
        i %= h->size;
    }
    return h->map[i]->payload;
}

void PrintHash(Hash *h)
{
    for (int i = 0; i < h->size; i++)
        if (h->map[i] != NULL)
            printf("%d => %lf\n", i, h->map[i]->payload);
}

void FreeHash(Hash *h)
{
    for (int i = 0; i < h->size; i++)
    {
        free(h->map[i]);
        h->map[i] = NULL;
    }
    free(h->map);
    h->map = NULL;
    h->used = h->size = 0;
}
