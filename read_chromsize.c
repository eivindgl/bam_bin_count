#include<stdio.h>

int main(int argc, char *argv[])
{
  const size_t BUF_SIZE = 1024;
  char buf[BUF_SIZE];
  char chrom_name[BUF_SIZE];
  size_t chrom_size;
  FILE *fp;
  if (argc == 1) {
    fprintf(stderr, "supply filename as arg\n");
    return 1;
  }
  fp = fopen(argv[1], "r");
  fgets(buf, BUF_SIZE, fp);
  
  while(fgets(buf, BUF_SIZE, fp)) {
    int rval = sscanf(buf, "%s %lu\n", chrom_name, &chrom_size);
    if (rval != 2) {
      fprintf(stderr, "Error parsing line: %s\n", buf);
      fprintf(stderr, "Aborting\n");
    }        
    printf("%d %s\t%lu\n", rval, chrom_name, chrom_size);
  }
  
  return 0;
}
