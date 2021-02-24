#ifndef _bioinfo_data_structures_h
#define _bioinfo_data_structures_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data_structures.h"

#define R 1.9858775229213840e-3

/* BED files */

typedef struct {
    char chrom[10];
    unsigned long int start;
    unsigned long int end;
    char name[80];
    int size;
} Region;

Region* read_next_region(FILE *bed);


/* VCF files */

typedef struct {
    char chrom[10];
    unsigned long int pos;
    char id[150];
    char ref[150];
    char alt[150];
    char format[2];
} Variant;

typedef struct {
    char *prefseq; // pointer to the refseq
    int n_variants;
    int n_diffs;
    CharsDynArray differences; // dynarray with one character for each alt allele
    IntsDynArray positions; // dynarray with one position for each alt allele
    BoolsDynArray haplotype; // dynarray with one haplotype for each snp
} Seq;

void InitSeq(Seq *seq, char *refseq);

void AddVariant(Seq *seq, bool allele, int pos, char alt);

char GetNthCharOfSeq(Seq *seq, int n);

int GetPosOfNthDiff(Seq *seq, int n);

char* GetFullSeq(Seq *seq, char *s);

void FreeSeq(Seq *seq);


/* PWM files */

#define PWMID_MAX_LENGTH 80

struct PWM {
    char id[PWMID_MAX_LENGTH];
    double *matrix;
    int length;
    double concentration;
    struct PWM *next;
};

typedef struct {
    struct PWM *head;
    struct PWM *tail;
    int size;
    double T;
} PWMList;

void InitPWMList(PWMList *list, double T);

void AddPWMToList(PWMList *list, char *id, double *matrix, int l);

int ReadNextPWM(FILE *pwm, PWMList *list);

void FreePWMList(PWMList *list);

#endif
