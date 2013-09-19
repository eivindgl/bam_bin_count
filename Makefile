TARGET = cabam
LIBS = -lm -lbam -lz -lpthread
CC = gcc
CFLAGS = -Wall -std=c99 -O2 # -DDEBUG -g
INCLUDES = -I samtools

.PHONY: clean all default

default: $(TARGET)
all: default

# OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
# HEADERS = $(wildcard *.h)

# %.o: %.c $(HEADERS)
# 	$(CC) $(CFLAGS) -c $< -o $@

# $(TARGET): $(OBJECTS)
# 	$(CC) $(OBJECTS) $(INCLUDES) -Wall $(LIBS) -o $@

$(TARGET): samtools/libbam.a read_bam.c
	$(CC) $(CFLAGS) -Lsamtools $(INCLUDES) read_bam.c $(LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f $(TARGET)
