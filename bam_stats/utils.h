#ifndef _UTILS_H_
#define _UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <stdint.h>
#include <stdarg.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>
#include <pthread.h>
#include <semaphore.h>
#include <getopt.h>
#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <sys/time.h>

/* color terminal */
#define VERBOSE_RED		"\033[1m\033[31m"
#define VERBOSE_GREEN		"\033[1m\033[32m"
#define VERBOSE_YELLOW		"\033[1m\033[33m"
#define VERBOSE_BLUE		"\033[1m\033[34m"
#define VERBOSE_CYAN		"\033[1m\033[36m"
#define VERBOSE_WHITE		"\033[0m"

#define MAX_INT32		2147483647
#define MIN_INT32		-2147483648

#define MASK32			4294967295ULL
#define SHIFT32			32

#define MAX_SZ_LINE		10000
#define MAX_DIR_LEN		4096

#define THREAD_STACK_SIZE	16777216

#define IS_LEFT			1
#define IS_RIGHT		0

/* Built in functions */
#define __swap(x, y) do {						       \
	assert(sizeof(x) == sizeof(y));					       \
	int8_t __swap_temp[sizeof(x)];					       \
	memcpy(__swap_temp, &(y), sizeof(x));				       \
	memcpy(&(y), &(x), sizeof(x));					       \
	memcpy(&(x), __swap_temp, sizeof(x));				       \
} while (0);

/* Built in macros */
#define __abs(x) 		((x) < 0 ? -(x) : (x))
#define __min(a, b) 		((a) < (b) ? (a) : (b))
#define __max(a, b) 		((a) > (b) ? (a) : (b))
#define __min3(a, b, c)		__min(__min((a), (b)), (c))
#define __max3(a, b, c)		__max(__max((a), (b)), (c))

#define __round_up_32(x) 	(--(x), (x) |= (x) >> 1,		       \
				 (x) |= (x) >> 2, (x) |= (x) >> 4,	       \
				 (x) |= (x) >> 8, (x) |= (x) >> 16, ++(x))

#define __verbose(fmt, args...) fprintf(stderr, fmt, ##args)

#define __warning(fmt, args...) do {					       \
	fprintf(stderr, "Warning: " fmt, ##args);			       \
	exit(EXIT_FAILURE);						       \
} while(0);

#define __error(fmt, args...) do {					       \
	fprintf(stderr, "Error: " fmt, ##args);				       \
	exit(EXIT_FAILURE);						       \
} while(0);

#define __ferror(fmt, args...) do {					       \
	fprintf(stderr, "Fatal error: " fmt, ##args);			       \
	exit(EXIT_FAILURE);						       \
} while(0);

#define __perror(fmt) do {						       \
	perror(fmt);							       \
	exit(EXIT_FAILURE);						       \
} while(0);

// extern int8_t ascii_table[256];
// extern char dna_char[5];

size_t xfread(void *ptr, size_t size, size_t nmemb, FILE *stream);
size_t xfwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);
char *xfgets(char *str, size_t n, FILE *stream);

double realtime();
void make_outdir(char *path);
char *reverse_str(char *seq);
char *reverse_complement(char *seq);

#endif /* _UTILS_H_ */