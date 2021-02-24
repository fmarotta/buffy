#ifndef _data_structures_h
#define _data_structures_h

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

/* Dynamic array of fixed-length strings */

#define STRINGSDYNARRAY_INITIAL_SIZE 1024
#define STRINGSDYNARRAY_MAX_LENGTH 50

typedef struct {
    char **a;
    int size;
    int used;
} StringsDynArray;

void InitStringsDynArray(StringsDynArray *sda);

void AddStringToDynArray(StringsDynArray *sda, const char *s);

char *GetNthString(StringsDynArray *sda, int n);

int GetSizeOfStringsDynArray(StringsDynArray *sda);

void FreeStringsDynArray(StringsDynArray *sda);


/* Dynamic array of characters */

#define CHARSDYNARRAY_INITIAL_SIZE 64

typedef struct {
    char *a;
    int size;
    int used;
} CharsDynArray;

void InitCharsDynArray(CharsDynArray *sda);

void AddCharToDynArray(CharsDynArray *sda, const char c);

char GetNthChar(CharsDynArray *sda, int n);

void FreeCharsDynArray(CharsDynArray *sda);

/* Dynamic array of ints */

#define INTSDYNARRAY_INITIAL_SIZE 64

typedef struct {
    int *a;
    int size;
    int used;
} IntsDynArray;

void InitIntsDynArray(IntsDynArray *sda);

void AddIntToDynArray(IntsDynArray *sda, const int i);

int GetNthInt(IntsDynArray *sda, int n);

int GetSizeOfIntsDynArray(IntsDynArray *sda);

void FreeIntsDynArray(IntsDynArray *sda);


/* Dynamic array of bools */

#define BOOLSDYNARRAY_INITIAL_SIZE 64

typedef struct {
    bool *a;
    int size;
    int used;
} BoolsDynArray;

void InitBoolsDynArray(BoolsDynArray *sda);

void AddBoolToDynArray(BoolsDynArray *sda, const bool b);

unsigned long long int GetBoolsDynArrayAsInt(BoolsDynArray *sda);

bool GetNthBool(BoolsDynArray *sda, int n);

void FreeBoolsDynArray(BoolsDynArray *sda);


/* Hash */

struct HashNode {
    unsigned long long int index;
    double payload;
};

typedef struct {
    struct HashNode** map;
    int size;
    int used;
} Hash;

static int hash(Hash *h, unsigned long long int index);

void InitHash(Hash *h, int n);

void AddToHash(Hash *h, unsigned long long int index, double e);

bool AlreadyInHash(Hash *h, unsigned long long int index);

double GetHash(Hash *h, int unsigned long long index);

void PrintHash(Hash *h);

void FreeHash(Hash *h);

#endif
