CC = gcc
CFLAGS = -std=c99 -O3 -Wall

# Path in which the header and library files are installed
PREFIX = .
TARGET = libprand.a    # libprand.so for the dynamic library

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
INC_DIR = $(ROOT_DIR)/header
SRC_DIR = $(ROOT_DIR)/src
SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(SRC_DIR)/%.o, $(SRCS))

ifeq ($(suffix $(TARGET)), .a)
  TARGET_MOD = 644
else
  ifeq ($(suffix $(TARGET)), .so)
    TARGET_MOD = 755
  else
    $(error unrecognized building target: $(TARGET))
  endif
endif

all: $(TARGET)

libprand.a: $(OBJS)
	ar rcs $(SRC_DIR)/$@ $^

libprand.so: $(OBJS)
	$(CC) $(CFLAGS) -shared -o $(SRC_DIR)/$@ $^

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $^ -o $@

clean:
	rm $(SRC_DIR)/*.o $(SRC_DIR)/*.a $(SRC_DIR)/*.so

install: $(TARGET)
	install -d $(PREFIX)/lib/
	install -m $(TARGET_MOD) $(SRC_DIR)/$(TARGET) $(PREFIX)/lib/
	install -d $(PREFIX)/include/
	install -m 644 $(INC_DIR)/prand.h $(PREFIX)/include/
