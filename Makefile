FLAGS := -Wall -Wextra -ansi -pedantic -Isrc/libzopfli
CFLAGS += $(FLAGS)
CXXFLAGS += $(FLAGS)

prefix ?= /usr/local
exec_prefix ?= $(prefix)
bindir ?= $(exec_prefix)/bin
libdir ?= $(exec_prefix)/lib
includedir ?= $(prefix)/include

LIBZOPFLI_SRC := src/libzopfli/blocksplitter.c src/libzopfli/cache.c\
                 src/libzopfli/deflate.c src/libzopfli/gzip_container.c\
                 src/libzopfli/hash.c src/libzopfli/katajainen.c\
                 src/libzopfli/lz77.c src/libzopfli/squeeze.c\
                 src/libzopfli/tree.c src/libzopfli/util.c\
                 src/libzopfli/zlib_container.c src/libzopfli/zopfli_lib.c
LIBZOPFLI_OBJ := $(LIBZOPFLI_SRC:.c=.o)
ZOPFLIBIN_SRC := src/zopfli/zopfli_bin.c
ZOPFLIBIN_OBJ := $(ZOPFLIBIN_SRC:.c=.o)
LODEPNG_SRC := src/zopflipng/lodepng/lodepng.cpp src/zopflipng/lodepng/lodepng_util.cpp
LODEPNG_OBJ := $(LODEPNG_SRC:.cpp=.o)
ZOPFLIPNG_SRC := src/zopflipng/zopflipng_bin.cc
ZOPFLIPNG_OBJ := $(ZOPFLIPNG_SRC:.cc=.o)
LIBZOPFLIPNG_SRC := src/zopflipng/zopflipng_lib.cc
LIBZOPFLIPNG_OBJ := $(LIBZOPFLIPNG_SRC:.cc=.o)

LIBZOPFLI := libzopfli.so.1.0.1
LIBZOPFLIPNG := libzopflipng.so.1.0.0
LIBZOPFLI_SONAME := libzopfli.so.1
LIBZOPFLIPNG_SONAME := libzopflipng.so.1

TARGETS := zopflipng zopfli libzopfli.so $(LIBZOPFLI_SONAME) $(LIBZOPFLIPNG_SONAME) $(LIBZOPFLI)

all: $(TARGETS)

# Zopfli shared libraries
libzopfli.so $(LIBZOPFLI_SONAME): $(LIBZOPFLI)
	ln -fs $< $@

libzopflipng.so $(LIBZOPFLIPNG_SONAME): $(LIBZOPFLIPNG)
	ln -fs $< $@

$(LIBZOPFLI): $(LIBZOPFLI_OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -shared -Wl,-soname,$(LIBZOPFLI_SONAME) -o $@ $^ -lm

$(LIBZOPFLIPNG): $(LIBZOPFLIPNG_OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -shared -Wl,-soname,$(LIBZOPFLIPNG_SONAME) -o $@ $^ -lm -L. -lzopfli

$(LIBZOPFLI_OBJ): %.o: %.c
	$(CC) $(CFLAGS) -fPIC -c -o $@ $<

$(LIBZOPFLIPNG_OBJ): %.o: %.cc
	$(CC) $(CFLAGS) -fPIC -c -o $@ $<

$(LODEPNG_SRC):
	git submodule update --init

# Zopfli binary
zopfli: $(ZOPFLIBIN_OBJ) libzopfli.so
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(ZOPFLIBIN_OBJ) -L. -lzopfli

# ZopfliPNG binary
zopflipng: $(ZOPFLIPNG_OBJ) $(LODEPNG_OBJ) libzopfli.so libzopflipng.so
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(ZOPFLIPNG_OBJ) $(LODEPNG_OBJ) -L. -lzopfli -lzopflipng

TEST_FILES := Makefile README README.zopflipng src/zopfli/zopfli_bin.c zopfli zopflipng libzopfli.so
test: $(TEST_FILES) zopfli zopflipng
	for file in $(TEST_FILES); do \
		LD_LIBRARY_PATH=. ./zopfli -c -v -- "$$file" | gzip -d | cmp - "$$file" || exit 1; \
	done
	@echo "Test succeeded."

# Remove all libraries and binaries
clean:
	$(RM) $(LIBZOPFLI_OBJ) $(ZOPFLIBIN_OBJ) $(LODEPNG_OBJ) $(ZOPFLIPNG_OBJ) $(LIBZOPFLIPNG_OBJ) $(TARGETS)

# install
install: all
	mkdir -p $(DESTDIR)$(bindir)
	install -m755 zopfli zopflipng $(DESTDIR)$(bindir)
	mkdir -p $(DESTDIR)$(libdir)
	install -m755 $(LIBZOPFLI) $(DESTDIR)$(libdir)
	cp -d libzopfli.so $(LIBZOPFLI) $(LIBZOPFLI_SONAME) \
	      libzopflipng.so $(LIBZOPFLIPNG) $(LIBZOPFLIPNG_SONAME) \
	      $(DESTDIR)$(libdir)
	mkdir -p $(DESTDIR)$(includedir)/zopfli
	install -m644 src/libzopfli/deflate.h src/libzopfli/zlib_container.h \
	              src/libzopfli/zopfli.h src/libzopfli/katajainen.h \
		      src/libzopfli/tree.h src/libzopfli/gzip_container.h \
		      src/libzopfli/cache.h src/libzopfli/squeeze.h \
		      src/libzopfli/lz77.h src/libzopfli/util.h \
		      src/libzopfli/blocksplitter.h src/libzopfli/hash.h \
		      $(DESTDIR)$(includedir)/zopfli
	mkdir -p $(DESTDIR)$(includedir)/zopflipng
	install -m644 src/zopflipng/zopflipng_lib.h $(DESTDIR)$(includedir)/zopflipng

.PHONY: clean test install
