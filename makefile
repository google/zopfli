make:
	gcc src/zopfli/*.c -O2 -W -Wall -Wextra -ansi -pedantic -lm -o zopfli

debug:
	gcc src/zopfli/*.c -g3 -lm -o zopfli
