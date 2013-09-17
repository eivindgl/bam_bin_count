TARGET = read_bam
LIBS = -lm -lbam -lz -lpthread
CC = gcc
CFLAGS = -g -Wall -std=c99 -O2
INCLUDES = -I../samtools

.PHONY: clean all default

default: $(TARGET)
all: default

# OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
# HEADERS = $(wildcard *.h)

# %.o: %.c $(HEADERS)
# 	$(CC) $(CFLAGS) -c $< -o $@

# $(TARGET): $(OBJECTS)
# 	$(CC) $(OBJECTS) $(INCLUDES) -Wall $(LIBS) -o $@

read_bam: ../samtools/libbam.a read_bam.c
	$(CC) $(CFLAGS) -L../samtools $(INCLUDES) read_bam.c $(LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f $(TARGET)
