make:
	gcc *.c -O2 -W -Wall -Wextra -ansi -pedantic -lm -o zopfli

debug:
	gcc *.c -g3 -lm -o zopfli
