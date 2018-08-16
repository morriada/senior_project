#
# Makefile
# Written by Adam Morris
#

### PROJECT SETTINGS ###
# Name of binary executable
BIN_NAME := senior_project
# Compiler used
CC ?= gcc
# Extention of source files used in the project
SRC_EXT = c
# Path to the source directory, relative to makefile
SRC_PATH = ./
# General compiler flags
COMPILE_FLAGS = -Wall -Wextra -O3 -lm -ffast-math -lrtlsdr -lpthread -lfftw3f -g -o

all: main.c control.c
	$(CC) $(COMPILE_FLAGS) $(BIN_NAME) main.c control.c

.PHONY: clean
clean:
	rm -f senior_project