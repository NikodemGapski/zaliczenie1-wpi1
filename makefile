CFLAGS=	-std=c17 -pedantic -Wall -Wextra -Wformat-security -Wduplicated-cond -Wfloat-equal\
		-Wshadow -Wconversion -Wjump-misses-init -Wlogical-not-parentheses -Wnull-dereference\
		-Wvla -Werror -fstack-protector-strong -fsanitize=undefined -fno-sanitize-recover -g\
		-fno-omit-frame-pointer -O1

test.e: test.c ary.c ary.h
		gcc ${CFLAGS} test.c ary.c -o test.e -lm

clean:
		rm -f *.e