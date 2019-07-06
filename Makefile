.PHONY: default help object executable all clean
CC = gcc

LD = $(CC)

SOURCE_C = $(wildcard src/*.c)
OBJECTS_C = $(patsubst src/%.c, bin/%, $(SOURCE_C))


CFLAGS = -std=gnu99 -O3 -Wall
EXECUTABLE = moleculardynamics.e

default: help

help:
	@echo $(SOURCE_C)
	@echo $(OBJECTS_C)
	@echo "\
Options:\n\n\
  make objects:       compiler makes objects for every *.c\n\
  make executable:    compiler makes executable\n\
  make all:           build all previous\n\
  make clean:         delete output files\n\
  make help:          display this help"

objects: $(OBJECTS_C)

executable: $(EXECUTABLE)

all: objects executable

bin/%: src/%.c
	$(CC) $(CFLAGS) -c $^ -o $@ -lm -largp

%.e: $(OBJECTS_C)
	$(LD) $^ -o $@ -lm -largp

clean:
	rm -rfv $(OBJECTS_C) $(EXECUTABLE1)
