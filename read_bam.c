#include <stdio.h>
#include <string.h>
#include "sam.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define MAX_CHROMS 100
#define BUF_SIZE 1024

typedef struct {
  char name[BUF_SIZE];
  float *bins;
  size_t len;
} chrom_meta;


int init_bins(const bam_header_t *h, const size_t bin_size, chrom_meta **cmp) {
  chrom_meta *cm = malloc(sizeof(chrom_meta) * h->n_targets);
  int i;
  for (i=0; i < h->n_targets; i++) {
    cm[i].len = h->target_len[i] / bin_size;
    if (h->target_len[i] % 2 != 0) {
      cm[i].len++;
    }
    cm[i].bins = malloc(sizeof(float) * cm[i].len);
    strncpy(cm[i].name, h->target_name[i], BUF_SIZE);
  }
  *cmp = cm;
  return i;
}

int read_bam(samfile_t *fp, const size_t bin_size, chrom_meta *cm, const size_t cm_len) {
  bam1_t *b = bam_init1();
  size_t start, cov, sidx;
  
  while (samread(fp, b) >= 0) {
    const bam1_core_t *c = &b->core;
    if (!(c->flag & 0x4)) {
      start = c->pos - 1; /* bam 1-based, bed is 0-based*/
      sidx = start / bin_size;
      cov = MIN(c->l_qseq, bin_size - start % bin_size);
      if (cov == c->l_qseq) {
        /* printf("%s [%lu] -> %d\n", cm[c->tid].name, sidx, 1); */
        cm[c->tid].bins[sidx] += 1;
      } else {
        float r = cov / (float)c->l_qseq;
        /* printf("HOLA %d %f\n", bin_size, cov, r); */
        /* printf("%s [%lu] -> %.2f\n", cm[c->tid].name, sidx, r); */
        cm[c->tid].bins[sidx] += r;
        cm[c->tid].bins[sidx] += 1 - r;
      }
      //printf("%s\t%lu\t%lu\n", cm[c->tid].name, start, end);
    }
  }
  bam_destroy1(b);

  return 0;
}  

int write_bed(const char *fname, const chrom_meta *cm,
              const size_t cm_len, const size_t bin_size) {
  FILE *fp;
  if ((fp=fopen(fname, "w")) == NULL) {
    fprintf(stderr, "unable to open output file: %s\n", fname);
  }
  for (int i=0; i < cm_len; i++) {
    for (int j=0; j < cm[i].len; j++) {
      fprintf(fp, "%s\t%lu\t%lu\t%f\n", cm[i].name,
              j * bin_size, (j + 1) * bin_size,
              cm[i].bins[j]);
    }
  }
  fclose(fp);
  return 0;
}

int main(int argc, char *argv[])
{
	samfile_t *fp;
        size_t bin_size;
	if (argc != 4) {
          fprintf(stderr, "Usage: %s <in.bam> <bin-size>\n", argv[0]);
		return 1;
	}
	if ((fp = samopen(argv[1], "rb", 0)) == 0) {
          fprintf(stderr, "bam2bed: Fail to open BAM file %s\n", argv[1]);
		return 1;
	}
        if (sscanf(argv[2], "%lu", &bin_size) != 1) {
          fprintf(stderr, "Error parsing bin size: %s\n", argv[2]);
          return 1;
        }

        chrom_meta *chrom_info;
        size_t chrom_info_len;

        chrom_info_len = init_bins(fp->header, MAX_CHROMS, &chrom_info);
        if (chrom_info_len < 0) {
          fprintf(stderr, "error reading chrom file. aborting ...\n");
          return 1;
        }

        if (read_bam(fp, bin_size, chrom_info, chrom_info_len) != 0) {
          fprintf(stderr, "Error reading bam file. Aborting ...\n");
          return 1;
        }
        samclose(fp);
        write_bed(argv[3], chrom_info, chrom_info_len, bin_size);
        free(chrom_info);
	return 0;
}
