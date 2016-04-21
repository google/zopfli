CC = gcc
CXX = g++

CFLAGS = -W -Wall -Wextra -ansi -pedantic -O2
CXXFLAGS = -W -Wall -Wextra -ansi -pedantic -O2
LDFLAGS = -lm

ZOPFLI_SRC := 	src/zopfli/blocksplitter.c src/zopfli/cache.c\
		src/zopfli/deflate.c src/zopfli/gzip_container.c\
                src/zopfli/hash.c src/zopfli/katajainen.c\
                src/zopfli/lz77.c src/zopfli/squeeze.c\
                src/zopfli/tree.c src/zopfli/util.c\
                src/zopfli/zlib_container.c
ZOPFLILIB_SRC := $(ZOPFLI_SRC) src/zopfli/zopfli_lib.c
ZOPFLILIB_OBJ := $(addsuffix .lo,$(basename $(ZOPFLILIB_SRC)))
ZOPFLIBIN_OBJ := $(addsuffix .o,$(basename $(ZOPFLILIB_SRC)))
LODEPNG_SRC := src/zopflipng/lodepng/lodepng.cpp src/zopflipng/lodepng/lodepng_util.cpp
ZOPFLIPNGLIB_SRC := src/zopflipng/zopflipng_lib.cc $(LODEPNG_SRC)
ZOPFLIPNGLIB_OBJ := $(addsuffix .lo,$(basename $(ZOPFLIPNGLIB_SRC)))
ZOPFLIPNGBIN_OBJ := $(addsuffix .o,$(basename $(ZOPFLIPNGLIB_SRC)))

.PHONY: clean all

%.lo: %.c
	$(CC) $(CFLAGS) -fPIC -c -o $@ $^

%.lo: %.c?*
	$(CXX) $(CXXFLAGS) -fPIC -c -o $@ $^

# Zopfli binary
zopfli: $(ZOPFLIBIN_OBJ) src/zopfli/zopfli_bin.o
	$(CC) $^ $(LDFLAGS) -o zopfli

# Zopfli shared library
libzopfli: $(ZOPFLILIB_OBJ)
	$(CC) $(ZOPFLILIB_OBJ) -shared -Wl,-soname,libzopfli.so.1 -o libzopfli.so.1.0.1

# ZopfliPNG binary
zopflipng: $(ZOPFLIPNGBIN_OBJ) $(ZOPFLIBIN_OBJ) src/zopflipng/zopflipng_bin.o
	$(CXX) $^ $(LDFLAGS) -o zopflipng

# ZopfliPNG shared library
libzopflipng: $(ZOPFLIPNGLIB_OBJ) $(ZOPFLILIB_OBJ)
	$(CXX) $(ZOPFLILIB_OBJ) $(ZOPFLIPNGLIB_OBJ) --shared -Wl,-soname,libzopflipng.so.1 -o libzopflipng.so.1.0.0

# Remove all libraries and binaries
clean:
	rm -f $(ZOPFLILIB_OBJ) libzopfli*
	rm -f $(ZOPFLIBIN_OBJ) zopfli
	rm -f $(ZOPFLIPNGBIN_OBJ) zopflipng
	rm -f $(ZOPFLIPNGLIB_OBJ)

all: zopfli libzopfli zopflipng libzopflipng
