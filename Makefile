CC ?= gcc
CXX ?= g++

override CFLAGS := -W -Wall -Wextra -ansi -pedantic -lm -O2 -Wno-unused-function -fPIC $(CFLAGS)
override CXXFLAGS := -W -Wall -Wextra -ansi -pedantic -O2 -fPIC $(CXXFLAGS)

ZOPFLILIB_SRC = src/zopfli/blocksplitter.c src/zopfli/cache.c\
                src/zopfli/deflate.c src/zopfli/gzip_container.c\
                src/zopfli/hash.c src/zopfli/katajainen.c\
                src/zopfli/lz77.c src/zopfli/squeeze.c\
                src/zopfli/tree.c src/zopfli/util.c\
                src/zopfli/zlib_container.c src/zopfli/zopfli_lib.c
ZOPFLILIB_OBJ := $(patsubst %.c,obj/%.o,$(ZOPFLILIB_SRC))
ZOPFLIBIN_SRC := src/zopfli/zopfli_bin.c
ZOPFLIBIN_OBJ := $(patsubst %.c,obj/%.o,$(ZOPFLIBIN_SRC))
LODEPNG_SRC := src/zopflipng/lodepng/lodepng.cpp src/zopflipng/lodepng/lodepng_util.cpp
LODEPNG_OBJ := $(patsubst %.cpp,obj/%.o,$(LODEPNG_SRC))
ZOPFLIPNGLIB_SRC := src/zopflipng/zopflipng_lib.cc
ZOPFLIPNGLIB_OBJ := $(patsubst %.cc,obj/%.o,$(ZOPFLIPNGLIB_SRC))
ZOPFLIPNGBIN_SRC := src/zopflipng/zopflipng_bin.cc
ZOPFLIPNGBIN_OBJ := $(patsubst %.cc,obj/%.o,$(ZOPFLIPNGBIN_SRC))

.PHONY: all libzopfli libzopflipng

all: zopfli libzopfli libzopfli.a zopflipng libzopflipng libzopflipng.a

obj/%.o: %.c
	@mkdir -p `dirname $@`
	$(CC) $(CFLAGS) -c $< -o $@

obj/%.o: %.cc
	@mkdir -p `dirname $@`
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/%.o: %.cpp
	@mkdir -p `dirname $@`
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Zopfli binary
zopfli: $(ZOPFLILIB_OBJ) $(ZOPFLIBIN_OBJ)
	$(CC) $^ $(CFLAGS) -o $@ $(LDFLAGS)

# Zopfli shared library
libzopfli: $(ZOPFLILIB_OBJ)
	$(CC) $^ $(CFLAGS) -shared -Wl,-soname,libzopfli.so.1 -o libzopfli.so.1.0.2 $(LDFLAGS)

# Zopfli static library
libzopfli.a: $(ZOPFLILIB_OBJ)
	ar rcs $@ $^

# ZopfliPNG binary
zopflipng: $(ZOPFLILIB_OBJ) $(LODEPNG_OBJ) $(ZOPFLIPNGLIB_OBJ) $(ZOPFLIPNGBIN_OBJ)
	$(CXX) $^ $(CFLAGS) -o $@ $(LDFLAGS)

# ZopfliPNG shared library
libzopflipng: $(ZOPFLILIB_OBJ) $(LODEPNG_OBJ) $(ZOPFLIPNGLIB_OBJ)
	$(CXX) $^ $(CFLAGS) --shared -Wl,-soname,libzopflipng.so.1 -o libzopflipng.so.1.0.2 $(LDFLAGS)

# ZopfliPNG static library
libzopflipng.a: $(LODEPNG_OBJ) $(ZOPFLIPNGLIB_OBJ)
	ar rcs $@ $^

# Remove all libraries and binaries
clean:
	rm -f zopflipng zopfli $(ZOPFLILIB_OBJ) $(ZOPFLIBIN_OBJ) $(LODEPNG_OBJ) $(ZOPFLIPNGLIB_OBJ) $(ZOPFLIPNGBIN_OBJ) libzopfli*
