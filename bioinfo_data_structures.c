#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data_structures.h"
#include "bioinfo_data_structures.h"

/* BED files */

Region* read_next_region(FILE *bed)
{
    Region *reg = malloc(sizeof(Region));
    if (fscanf(bed, "%s%ld%ld%s",
                reg->chrom, &reg->start, &reg->end, reg->name) == 4)
    {
        while (fgetc(bed) != '\n')
            ;
        reg->size = reg->end - reg->start;
        return reg;
    }
    return NULL;
}


/* VCF files */

void InitSeq(Seq *seq, char *refseq)
{
    seq->prefseq = refseq;
    seq->n_diffs = 0;
    seq->n_variants = 0;
    InitBoolsDynArray(&seq->haplotype);
    InitIntsDynArray(&seq->positions);
    InitCharsDynArray(&seq->differences);
}

void AddVariant(Seq *seq, bool allele, int pos, char alt)
{
    AddBoolToDynArray(&seq->haplotype, allele);
    seq->n_variants++;
    if (allele)
    {
        AddIntToDynArray(&seq->positions, pos);
        AddCharToDynArray(&seq->differences, alt);
        seq->n_diffs++;
    }
}

char GetNthCharOfSeq(Seq *seq, int n)
{
    // TODO maybe binary search
    for (int i = 0; i < seq->n_diffs; i++)
        if (n == GetNthInt(&seq->positions, i))
            return GetNthChar(&seq->differences, i);
    return seq->prefseq[n];
}

int GetPosOfNthDiff(Seq *seq, int n)
{
    if (n > seq->n_diffs)
        return -1;
    return GetNthInt(&seq->positions, n);
}

char* GetFullSeq(Seq *seq, char *s)
{
    int j = 0;

    if (!seq->n_diffs)
        strcpy(s, seq->prefseq);
    else
        for (int i = 0; i < strlen(seq->prefseq); i++)
        {
            if (j < seq->n_diffs && i == GetNthInt(&seq->positions, j))
            {
                s[i] = GetNthChar(&seq->differences, j++);
                continue;
            }
            s[i] = seq->prefseq[i];
        }
    s[strlen(seq->prefseq)] = '\0';

    return s;
}

void FreeSeq(Seq *seq)
{
    FreeBoolsDynArray(&seq->haplotype);
    FreeIntsDynArray(&seq->positions);
    FreeCharsDynArray(&seq->differences);
    seq->prefseq = NULL;
    seq->n_diffs = 0;
    seq->n_variants = 0;
}


/* PWM files */

void InitPWMList(PWMList *list, double T)
{
    list->head = list->tail = NULL;
    list->size = 0;
    list->T = T;
}

void AddPWMToList(PWMList *list, char *id, double *matrix, int l)
{
    struct PWM *pwm = malloc(sizeof(struct PWM));

    strncpy(pwm->id, id, PWMID_MAX_LENGTH);
    pwm->matrix = matrix;
    pwm->length = l;
    pwm->next = NULL;

    // Here we replace each entry in the matrix with its deltaG
    // deltaG for base b at position i is given by
    // log2(P_bi/background_b)) * 300 * R * log(2)
    for (int i = 0; i < 4 * l; i++)
        pwm->matrix[i] *= 300 * R * log(2);

    // Here we keep only the max deltaG for each position
    double k = 0;
    double max_deltaG;
    for (int m = 0; m < l; m++)
    {
        max_deltaG = pwm->matrix[4 * m];
        for (int i = 1; i < 4; i++)
            if (pwm->matrix[4 * m + i] > max_deltaG)
                max_deltaG = pwm->matrix[4 * m + i];
        k += max_deltaG;
    }
    pwm->concentration = exp(-k / (R * list->T));

    if (list->size == 0)
    {
        list->head = pwm;
        list->tail = pwm;
    }
    else
    {
        list->tail->next = pwm;
        list->tail = pwm;
    }
    list->size++;
}

int ReadNextPWM(FILE *pwm, PWMList *list)
{
    char id[PWMID_MAX_LENGTH];
    double *matrix = malloc(sizeof(double) * 4);
    int i = 0;
    char ch;

    if ((ch = fgetc(pwm)) != '>')
        return 0;
    fscanf(pwm, "%s", id);
    while(fgetc(pwm) != '\n')
        ;
    fscanf(pwm, "%lf%lf%lf%lf", &matrix[4*i], &matrix[4*i+1], &matrix[4*i+2], &matrix[4*i+3]);
    while ((ch = fgetc(pwm)) != '>' && ch != EOF)
    {
        matrix = realloc(matrix, sizeof(double) * 4 * (++i + 1));
        fscanf(pwm, "%lf%lf%lf%lf", &matrix[4*i], &matrix[4*i+1], &matrix[4*i+2], &matrix[4*i+3]);
    }
    ungetc(ch, pwm);

    AddPWMToList(list, id, matrix, i);
    return 1;
}

void FreePWMList(PWMList *list)
{
    struct PWM *p;
    struct PWM *pnext;
    p = list->head;
    while (p != NULL)
    {
        pnext = p->next;
        free(p->matrix);
        free(p);
        p = pnext;
    }
    list->head = list->tail = NULL;
    list->size = 0;
}
