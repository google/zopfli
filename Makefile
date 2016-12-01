CC = gcc
CXX = g++
AR = ar

CFLAGS = -W -Wall -Wextra -ansi -pedantic -lm -O2 -Wno-unused-function
CXXFLAGS = -W -Wall -Wextra -ansi -pedantic -O2
ARFLAGS = rv
ZOPFLI_TARGET = libzopfli.a
ZOPFLIPNG_TARGET = libzopflipng.a

ZOPFLILIB_SRC = src/zopfli/blocksplitter.c src/zopfli/cache.c\
                src/zopfli/deflate.c src/zopfli/gzip_container.c\
                src/zopfli/hash.c src/zopfli/katajainen.c\
                src/zopfli/lz77.c src/zopfli/squeeze.c\
                src/zopfli/tree.c src/zopfli/util.c\
                src/zopfli/zlib_container.c src/zopfli/zopfli_lib.c
ZOPFLILIB_OBJ := $(patsubst src/zopfli/%.c,%.o,$(ZOPFLILIB_SRC))
ZOPFLIBIN_SRC := src/zopfli/zopfli_bin.c
LODEPNG_SRC := src/zopflipng/lodepng/lodepng.cpp src/zopflipng/lodepng/lodepng_util.cpp
LODEPNG_OBJ := $(patsubst src/zopflipng/lodepng/%.cpp,%.o,$(LODEPNG_SRC))
ZOPFLIPNGLIB_SRC := src/zopflipng/zopflipng_lib.cc
ZOPFLIPNGLIB_OBJ := $(patsubst src/zopflipng/%.cc,%.o,$(ZOPFLIPNGLIB_SRC))
ZOPFLIPNGBIN_SRC := src/zopflipng/zopflipng_bin.cc

.PHONY: zopfli zopflipng

# Zopfli binary
zopfli:
	$(CC) $(ZOPFLILIB_SRC) $(ZOPFLIBIN_SRC) $(CFLAGS) -o zopfli
	$(AR) $(ARFLAGS) $(ZOPFLI_TARGET) $(ZOPFLILIB_OBJ)

# Zopfli shared library
libzopfli:
	$(CC) $(ZOPFLILIB_SRC) $(CFLAGS) -fPIC -c
	$(CC) $(ZOPFLILIB_OBJ) $(CFLAGS) -shared -Wl,-soname,libzopfli.so.1 -o libzopfli.so.1.0.1

# ZopfliPNG binary
zopflipng:
	$(CC) $(ZOPFLILIB_SRC) $(CFLAGS) -c
	$(CXX) $(ZOPFLIPNGLIB_SRC) $(LODEPNG_SRC) $(CXXFLAGS) -c
	$(CXX) $(ZOPFLILIB_OBJ) $(LODEPNG_SRC) $(ZOPFLIPNGLIB_SRC) $(ZOPFLIPNGBIN_SRC) $(CFLAGS) -o zopflipng
	$(AR) $(ARFLAGS) $(ZOPFLIPNG_TARGET) $(ZOPFLILIB_OBJ) $(LODEPNG_OBJ) $(ZOPFLIPNGLIB_OBJ)

# ZopfliPNG shared library
libzopflipng:
	$(CC) $(ZOPFLILIB_SRC) $(CFLAGS) -fPIC -c
	$(CXX) $(ZOPFLILIB_OBJ) $(LODEPNG_SRC) $(ZOPFLIPNGLIB_SRC) $(CFLAGS) -fPIC --shared -Wl,-soname,libzopflipng.so.1 -o libzopflipng.so.1.0.0

# Remove all libraries and binaries
clean:
	rm -f zopflipng zopfli $(ZOPFLILIB_OBJ) libzopfli*
