# Description:
#   Zopfli is a zlib- and deflate-compatible compression algorithm.
#   Zopfli prefers compression ratio over compression speed.

STRICT_C_OPTIONS = [
    "-O3",
    "-W",
    "-Wall",
    "-Wextra",
    "-Wno-unused-function",
    "-ansi",
    "-fPIC",
    "-pedantic",
]

cc_binary(
    name = "zopfli",
    srcs = glob([
        "src/zopfli/*.c",
        "src/zopfli/*.h",
    ]),
    copts = STRICT_C_OPTIONS,
    visibility = ["//visibility:public"],
)

cc_binary(
    name = "zopflipng",
    srcs = glob(
        [
            "src/zopflipng/*.cc",
            "src/zopflipng/*.h",
            "src/zopfli/*.c",
            "src/zopfli/*.h",
            "src/zopflipng/lodepng/*.cpp",
            "src/zopflipng/lodepng/*.h",
        ],
        exclude = ["src/zopfli/zopfli_bin.c"],
    ),
    copts = STRICT_C_OPTIONS,
    visibility = ["//visibility:public"],
)
