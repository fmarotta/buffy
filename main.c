/* Take a VCF, a BED, and a reference FASTA, then output smaller FASTAS
 * (one for each individual in the VCF) with the respective sequences at
 * the regions in the BED.
 */

// TODO output stats like avg number of variants per regreg

// TODO maybe the pwm is log 2 not ln

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#include "data_structures.h"
#include "bioinfo_data_structures.h"

#define DEBUG false
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) < (Y)) ? (Y) : (X))

#define PLOIDY 2

char* rev_compl_string(char *forward)
{
    int n = strlen(forward);
    char *r = malloc(sizeof(char) * (n + 1));
    if (r == NULL)
        return NULL;
    for (int i = 0; i < n; i++)
    {
        switch (forward[n - (i + 1)])
        {
        case 'A':
            r[i] = 'T';
            break;
        case 'C':
            r[i] = 'G';
            break;
        case 'G':
            r[i] = 'C';
            break;
        case 'T':
            r[i] = 'A';
            break;
        default:
            return NULL;
        }
    }
    r[n] = '\0';
    return r;
}

char rev_compl_char(char base)
{
    switch (base)
    {
    case 'A':
        return 'T';
    case 'C':
        return 'G';
    case 'G':
        return 'C';
    case 'T':
        return 'A';
    default:
        return 'X';
    }
}

int DNA2int(char base)
{
    int i;
    switch(base)
    {
        case 'A':
            i = 0;
            break;
        case 'C':
            i = 1;
            break;
        case 'G':
            i = 2;
            break;
        case 'T':
            i = 3;
            break;
        default:
            return -1;
    }
    return i;
}

char* read_refseq(FILE *fasta, int n)
{
    unsigned long int start;
    char ch;
    int i = 0;
    char *refseq = malloc(sizeof(char) * (n + 1));
    if(refseq == NULL)
        return NULL;

    start = ftell(fasta);
    while (i < n)
    {
        ch = toupper(fgetc(fasta));
        if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'N')
            refseq[i++] = ch;
        else if (ch != '\n')
            refseq[i++] = 'X';
    }
    refseq[i] = '\0';

    fseek(fasta, start, SEEK_SET);
    return refseq;
}

double ComputeDB(double *k, int from, int to, double T, double conc)
{
    double DB = 0;
    for (int j = from; j < to; j++)
        DB += 1 / (1 + exp(-k[j] / (R * T)) / conc);
    return DB;
}

double ComputeDeltaDB(double *kd_xi, int reg_size, struct PWM *p, Seq *seq, double T)
{
    double kd_xi_tmp[reg_size - p->length + 1];
    for (int j = 0; j < reg_size - p->length + 1; j++)
        kd_xi_tmp[j] = kd_xi[j];
    for (int j = 0; j < seq->n_diffs; j++)
    {
        int k = GetPosOfNthDiff(seq, j);
        for (int l = MAX(k - p->length + 1, 0); l <= MIN(reg_size - p->length, k); l++)
        {
            kd_xi_tmp[l] = 0;
            for (int m = 0; m < p->length; m++)
                kd_xi_tmp[l] += p->matrix[4 * m + DNA2int(GetNthCharOfSeq(seq, l + m))];
        }
    }
    return ComputeDB(kd_xi_tmp, 0, reg_size - p->length + 1, T, p->concentration);
}

int main(int argc, char *argv[])
{
    FILE *bed, *fasta, *vcf, *pwm, *out;
    Region *reg;
    PWMList pwmlist;
    struct PWM *p;
    char *refseq_forward;
    char *refseq_reverse;
    StringsDynArray samples;
    Variant var;
    int n_variants;
    int n_samples;
    unsigned long int refpos = 0;
    unsigned long int vcfpos = 0;
    double *kd_xi_forward;
    double *kd_xi_reverse;
    Hash DB_forward;
    Hash DB_reverse;
    double T = 310;

    if ((bed = fopen(argv[1], "r")) == NULL)
    {
        fprintf(stderr, "Error: cannot open %s.\n", argv[1]);
        exit(EXIT_FAILURE);
    }
    if ((fasta = fopen(argv[2], "r")) == NULL)
    {
        fprintf(stderr, "Error: cannot open %s.\n", argv[2]);
        exit(EXIT_FAILURE);
    }
    if ((vcf = fopen(argv[3], "r")) == NULL)
    {
        fprintf(stderr, "Error: cannot open %s.\n", argv[3]);
        exit(EXIT_FAILURE);
    }
    if ((pwm = fopen(argv[4], "r")) == NULL)
    {
        fprintf(stderr, "Error: cannot open %s.\n", argv[4]);
        exit(EXIT_FAILURE);
    }
    if ((out = fopen(argv[5], "w")) == NULL)
    {
        fprintf(stderr, "Error: cannot open %s.\n", argv[5]);
        exit(EXIT_FAILURE);
    }

    // Get number of individuals
    char tmp_iid[STRINGSDYNARRAY_MAX_LENGTH];
    InitStringsDynArray(&samples);
    while (fgetc(vcf) == '#')
    {
        if (fgetc(vcf) == '#')
        {
            while (fgetc(vcf) != '\n')
                ;
        }
        else
        {
            fscanf(vcf, "%*s%*s%*s%*s%*s%*s%*s%*s%*s");
            char ch;
            while ((ch = fgetc(vcf)) != '\n')
            {
                fscanf(vcf, "%s", tmp_iid);
                AddStringToDynArray(&samples, tmp_iid);
            }
            break;
        }
    }
    vcfpos = ftell(vcf);
    // Allocate sequence for each individual. For each individual we use
    // PLOIDY consecutive seq elements, one for each allele
    n_samples = GetSizeOfStringsDynArray(&samples);
    Seq seq_forward[PLOIDY * n_samples];
    Seq seq_reverse[PLOIDY * n_samples];
#if DEBUG
    printf("SAMPLES: %d\n", n_samples);
    for (int i = 0; i < n_samples; i++)
        printf("%s\n", GetNthString(&samples, i));
#endif

    // Advance in the fasta
    refpos = 0;
    if (fgetc(fasta) != '>')
    {
        fprintf(stderr, "Error: malformed FASTA file %s.\n", argv[2]);
        exit(EXIT_FAILURE);
    }
    while (fgetc(fasta) != '\n')
        ;

    // Read all the PWM
    InitPWMList(&pwmlist, T);
    while (ReadNextPWM(pwm, &pwmlist)) // TODO good chance to compute some stats about e.g. length of pwms
        ;

    // Process each regreg
    while ((reg = read_next_region(bed)) != NULL)
    {
#if DEBUG
        printf("BED: %s\t%ld\t%ld\t%s\n",
                reg->chrom, reg->start, reg->end, reg->name);
#endif

        // read the ref fasta at the right position
        while (refpos < reg->start)
            if (fgetc(fasta) != '\n')
                refpos++;
        refseq_forward = read_refseq(fasta, reg->size);
        refseq_reverse = rev_compl_string(refseq_forward);
        if (refseq_forward == NULL || refseq_reverse == NULL)
        {
            fprintf(stderr, "Warning: problems reading the sequence; skipping.\n");
            continue;
        }
#if DEBUG
        //printf("REFPOS: %ld\n", refpos);
        printf("REFSEQ: %s\n", refseq_forward);
        printf("REFSEQ (RC): %s\n", refseq_reverse);
#endif
        for (int i = 0; i < PLOIDY * n_samples; i++)
        {
            InitSeq(&seq_forward[i], refseq_forward);
            InitSeq(&seq_reverse[i], refseq_reverse);
        }

        // read the vcf and add variants to individuals
        fseek(vcf, vcfpos, SEEK_SET);
        while (fscanf(vcf, "%s%ld%s%s%s%*s%*s%*s%s",
                       var.chrom, &var.pos, var.id,
                       var.ref, var.alt, var.format) == 6)
        {
            if (var.pos <= reg->start)
            {
                while (fgetc(vcf) != '\n')
                    ;
                vcfpos = ftell(vcf);
                continue;
            }
            if (var.pos > reg->end)
            {
                while (fgetc(vcf) != '\n')
                    ;
                break;
            }
            // Ignore indels
            if (strlen(var.ref) > 1 || strlen(var.alt) > 1)
            {
                while (fgetc(vcf) != '\n')
                    ;
                continue;
            }
            // TODO ignore variants that are closer than p->length to
            // the boundary of the regreg
#if DEBUG
            printf("VCF: %s %ld %s %s\n", var.id, var.pos, var.ref, var.alt);
            //printf("In the refseq, at position %d, there is %c.\n",
            //       var.pos, toupper(refseq[var.pos-1 - reg->start]));
            //if (toupper(refseq[var.pos-1 - reg->start]) != var.ref[0])
            //    printf("something is fishy here.");
#endif

            // Store the differences in seq among individuals
            bool allele[PLOIDY];
            for (int i = 0; i < n_samples; i++)
            {
            for (int j = 0; j < PLOIDY; j++)
            {
                fgetc(vcf);
                if (fgetc(vcf) == '0')
                    allele[j] = 0;
                else
                    allele[j] = 1;

                if (allele[j] == 0)
                {
                    if (refseq_forward[var.pos - 1 - reg->start] != *var.ref)
                        fprintf(stderr, "Warning: refseq != VCF\n");
                    if (refseq_reverse[reg->end - var.pos] != rev_compl_char(*var.ref))
                        fprintf(stderr, "Warning: refseq != VCF\n");
                }
                AddVariant(&seq_forward[PLOIDY * i + j], allele[j], var.pos - 1 - reg->start, *var.alt);
                AddVariant(&seq_reverse[PLOIDY * i + j], allele[j], reg->end - var.pos, rev_compl_char(*var.alt));
            }
            }
        }

        // Compute kd for each pos in the refseq
        p = pwmlist.head;
        n_variants = seq_forward[0].n_variants;
        if (n_variants >= 8 * (int) sizeof(unsigned long long int))
        {
            // We use unsigned long long ints to store the index from
            // GetBoolsDynArrayAsInt(); if there are too many variants
            // it will overflow
            fprintf(stderr, "Warning: %s has more than %d variants; skipping it.\n", reg->name, 8 * (int) sizeof(unsigned long long int) - 1);
            continue;
        }
        while (p != NULL)
        {
#if DEBUG
    printf("PWM: %s\n", p->id);
#endif
            InitHash(&DB_forward, (int) MIN(2 * n_samples, pow(2, n_variants)));
            InitHash(&DB_reverse, (int) MIN(2 * n_samples, pow(2, n_variants)));
            // Compute the Kd's
            kd_xi_forward = calloc(reg->size - p->length + 1, sizeof(double) * (reg->size - p->length + 1));
            kd_xi_reverse = calloc(reg->size - p->length + 1, sizeof(double));
            if (kd_xi_reverse == NULL || kd_xi_forward == NULL)
            {
                fprintf(stderr, "kd_xi: memory allocation failed.\n");
                exit(EXIT_FAILURE);
            }
            for (int i = 0; i < reg->size - p->length + 1; i++)
                for (int j = 0; j < p->length; j++)
                {
                    kd_xi_forward[i] += p->matrix[4 * j + DNA2int(refseq_forward[i + j])];
                    kd_xi_reverse[i] += p->matrix[4 * j + DNA2int(refseq_reverse[i + j])];
                }

            for (int i = 0; i < PLOIDY * n_samples; i++)
            {
                if (GetBoolsDynArrayAsInt(&seq_forward[i].haplotype) != GetBoolsDynArrayAsInt(&seq_reverse[i].haplotype))
                    return 1;
                if (!AlreadyInHash(&DB_forward, GetBoolsDynArrayAsInt(&seq_forward[i].haplotype)))
                {
#if DEBUG
    printf("Kd: computing forward for individual %d\n", i / PLOIDY);
#endif
                    AddToHash(&DB_forward,
                            GetBoolsDynArrayAsInt(&seq_forward[i].haplotype),
                            ComputeDeltaDB(kd_xi_forward, reg->size, p, &seq_forward[i], T));
                }
                if (!AlreadyInHash(&DB_reverse, GetBoolsDynArrayAsInt(&seq_reverse[i].haplotype)))
                {
#if DEBUG
    printf("Kd: computing reverse for individual %d\n", i / PLOIDY);
#endif
                    AddToHash(&DB_reverse,
                            GetBoolsDynArrayAsInt(&seq_reverse[i].haplotype),
                            ComputeDeltaDB(kd_xi_reverse, reg->size, p, &seq_reverse[i], T));
                }
                fprintf(out, "%s\t%s\t%s@%d\t%.*e\n", reg->name,
                        GetNthString(&samples, i / PLOIDY),
                        p->id,
                        i % PLOIDY,
                        DECIMAL_DIG,
                        GetHash(&DB_forward, GetBoolsDynArrayAsInt(&seq_forward[i].haplotype)) + GetHash(&DB_reverse, GetBoolsDynArrayAsInt(&seq_reverse[i].haplotype)));
            }

            free(kd_xi_forward);
            free(kd_xi_reverse);
            FreeHash(&DB_forward);
            FreeHash(&DB_reverse);
            p = p->next;
        }

        for (int i = 0; i < n_samples; i++)
        {
        for (int j = 0; j < PLOIDY; j++)
        {
            // Print the FASTA
            // TODO add an option for that
            // fprintf(out, ">%s@%s@%d\n", reg->name, GetNthString(&samples, i), j);
            // char tmp[strlen(seq[PLOIDY * i + j].prefseq) + 1];
            // GetFullSeq(&seq[PLOIDY * i + j], tmp);
            // fprintf(out, "%s\n", tmp);

            // Free the pointers
            FreeSeq(&seq_forward[PLOIDY * i + j]);
            FreeSeq(&seq_reverse[PLOIDY * i + j]);
        }
        }
        free(refseq_forward);
        free(refseq_reverse);
        free(reg);
    }

    FreeStringsDynArray(&samples);
    FreePWMList(&pwmlist);
    fclose(bed);
    fclose(fasta);
    fclose(vcf);
    fclose(pwm);
    fclose(out);
    return 0;
}
