#include <stdio.h>
#include <string.h>
#include "sam.h"
#include <assert.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define MAX_CHROMS 100
#define BUF_SIZE 1024

typedef struct {
  char name[BUF_SIZE];
  float *bins;
  size_t len;
} chrom_meta;

typedef struct {
  char name[BUF_SIZE];
  size_t len;
} chrom_size;

int read_chrom_sizes(char *fname, chrom_size *cs, const size_t cs_len) {
  char buf[BUF_SIZE];
  FILE *fp;
  int i = 0;

  fp = fopen(fname, "r");
  if (fp == NULL || fgets(buf, BUF_SIZE, fp) == NULL) {
    fprintf(stderr, "error opening chromsize file: %s\n", fname); 
    return -1;
  }
  
  while(fgets(buf, BUF_SIZE, fp)) {
    if (i >= cs_len) {
      fprintf(stderr, "Error: too many chromosomes (increase buffer size)\n");
      return -1;
    }
    int rval = sscanf(buf, "%s %lu\n", cs[i].name, &cs[i].len);
    if (rval != 2) {
      fprintf(stderr, "Error parsing line: %s\n", buf);
      fprintf(stderr, "Aborting\n");
      return -1;
    }
    i++;
  }
  fclose(fp);
  return i;
}

const chrom_size * find_chrom(const chrom_size *cs, const size_t cs_len, const char *name) {
  for (int i = 0; i < cs_len; ++i) {
    if (strcmp(cs[i].name, name) == 0) {
      return &cs[i];
    }
  }
  return NULL;
}

chrom_meta * find_chrom_meta(chrom_meta ** cs, size_t cs_len, const char * name) {
  for (int i = 0; i < cs_len; ++i) {
    if ((cs[i] != NULL) && (strcmp(cs[i]->name, name) == 0)) {
      return cs[i];
    }
  }
  return NULL;
}


int init_bins(const bam_header_t *h, chrom_meta **cm, const size_t cm_len, const size_t bin_size,
              chrom_size *chromnames, int nchroms) {
  
  int i;
  for (i=0; i < h->n_targets; i++) {
    assert(i<cm_len);
    const chrom_size *cur = find_chrom(chromnames, nchroms, h->target_name[i]);
    if (cur == NULL) {
      cm[i] = NULL;
    } else {
      cm[i] = malloc(sizeof(chrom_meta));
      cm[i]->len = cur->len / bin_size;
      assert(cm[i]->len < 3000);
      if (cur->len % 2 != 0) {
        cm[i]->len++;
      }
      
      cm[i]->bins = malloc(sizeof(float) * cm[i]->len);
      strncpy(cm[i]->name, cur->name, BUF_SIZE);
    }
  }
  
  return i;
}

int read_bam(samfile_t *fp, const size_t bin_size, chrom_meta **cm, const size_t cm_len) {
  bam1_t *b = bam_init1();
  size_t start, cov, sidx;
  
  while (samread(fp, b) >= 0) {
    const bam1_core_t *c = &b->core;
    if (!(c->flag & 0x4) && cm[c->tid] != NULL) {
      start = c->pos - 1; /* bam 1-based, bed is 0-based*/
      sidx = start / bin_size;
      cov = MIN(c->l_qseq, bin_size - start % bin_size);
      if (cov == c->l_qseq) {
        /* printf("%s [%lu] -> %d\n", cm[c->tid].name, sidx, 1); */
        cm[c->tid]->bins[sidx] += 1;
      } else {
        float r = cov / (float)c->l_qseq;
        /* printf("HOLA %d %f\n", bin_size, cov, r); */
        /* printf("%s [%lu] -> %.2f\n", cm[c->tid].name, sidx, r); */
        cm[c->tid]->bins[sidx] += r;
        cm[c->tid]->bins[sidx] += 1 - r;
      }
      //printf("%s\t%lu\t%lu\n", cm[c->tid].name, start, end);
    }
  }
  bam_destroy1(b);

  return 0;
}

void free_chrom_meta(chrom_meta **cm, size_t cm_len) {
  for (int i = 0; i < cm_len; ++i)
    {
      if (cm[i] != NULL) {
        free(cm[i]->bins);
        free(cm[i]);
      }
    }
  free(cm);
}

int write_bed(const char *fname, 
              chrom_meta **ip, const size_t ip_len,
              chrom_meta **input, const size_t input_len,
              const size_t bin_size) {
  FILE *fp;
  if ((fp=fopen(fname, "w")) == NULL) {
    fprintf(stderr, "unable to open output file: %s\n", fname);
  }
  for (int i=0; i < ip_len; i++) {
    if (ip[i] == NULL) {
      continue;
    }
    chrom_meta *input_chrom = find_chrom_meta(input, input_len, ip[i]->name);

    for (int j=0; j < ip[i]->len; j++) {
      float input_cnt = input_chrom == NULL ? 0 : input_chrom->bins[j];
      fprintf(fp, "%s\t%lu\t%lu\t%f\t%f\n", ip[i]->name,
              j * bin_size, (j + 1) * bin_size,
              ip[i]->bins[j], input_cnt);
    }
  }

  /* print all chroms exclusive to input */
  for (int i=0; i < input_len; i++) {
    if (input[i] == NULL) {
      continue;
    }
    chrom_meta *ip_chrom = find_chrom_meta(ip, ip_len, input[i]->name);
    if (ip_chrom != NULL) {
      /* both present chroms already written */
      continue;
    }
    for (int j=0; j < input[i]->len; j++) {
      fprintf(fp, "%s\t%lu\t%lu\t%f\t%f\n", input[i]->name,
              j * bin_size, (j + 1) * bin_size,
              0.0, input[i]->bins[j]);
    }
  }
  
  fclose(fp);
  return 0;
}

int count_bam(char *samfilename, chrom_meta ***cmp, const size_t bin_size,
              chrom_size *chromnames, int nchroms) {
  samfile_t *fp;
  if ((fp = samopen(samfilename, "rb", 0)) == 0) {
    fprintf(stderr, "Failed to open BAM file %s\n", samfilename);
    return 1;
  }
  size_t cm_len = fp->header->n_targets;
  chrom_meta **cm = malloc(sizeof(chrom_meta *) * cm_len);

  int rval = init_bins(fp->header, cm, cm_len, bin_size, chromnames, nchroms);
  if (rval < 0) {
    fprintf(stderr, "error initializing bins ...\n");
    return -1;
  }

  if (read_bam(fp, bin_size, cm, cm_len) != 0) {
    fprintf(stderr, "Error reading bam file. Aborting ...\n");
    return -1;
  }
  samclose(fp);
  *cmp = cm;
  return cm_len;
}

int main(int argc, char *argv[])
{
  if (argc != 6) {
    fprintf(stderr, "Usage: %s <chromsize> <ip.bam> <input.bam> <bin-size> <output_file>\n", argv[0]);
    return 1;
  }
  
  chrom_size *chromnames = malloc(sizeof(chrom_size) * MAX_CHROMS);
  int nchroms = read_chrom_sizes(argv[1], chromnames, MAX_CHROMS);
  if (nchroms < 0) {
    fprintf(stderr, "Unable to read chromsizes from file: %s\n.", argv[1]);
    return 1;
  }
  
  size_t bin_size;
  if (sscanf(argv[4], "%lu", &bin_size) != 1) {
    fprintf(stderr, "Error parsing bin size: %s\n", argv[4]);
    return 1;
  }
  
  chrom_meta **ip;
  int ip_len = count_bam(argv[2], &ip, bin_size, chromnames, nchroms);
  if (ip_len < 0) {
    fprintf(stderr, "problems reading ip file=%s\n", argv[2]);
    return 1;
  }
  chrom_meta **input;
  int input_len = count_bam(argv[3], &input, bin_size, chromnames, nchroms);
  if (ip_len < 0) {
    fprintf(stderr, "problems reading input file=%s\n", argv[3]);
    return 1;
  }    

  write_bed(argv[5], ip, ip_len, input, input_len, bin_size);
  free_chrom_meta(ip, ip_len);
  free_chrom_meta(input, input_len);
  free(chromnames);
  return 0;
}
