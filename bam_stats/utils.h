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
#define RED			"\033[1m\033[31m"
#define GREEN			"\033[1m\033[32m"
#define YELLOW			"\033[1m\033[33m"
#define BLUE			"\033[1m\033[34m"
#define CYAN			"\033[1m\033[36m"
#define WHITE			"\033[0m"

#define MAX_INT32		2147483647
#define MIN_INT32		-2147483648

#define MASK32			4294967295ULL

#define BUFSZ			4096

#define THREAD_STACK_SIZE	16777216

#define FORWARD			0
#define REVERSE			1
#define LEFT			0
#define RIGHT			1

/*
 * Built in macros
 */

#define __abs(x) 		((x) < 0 ? -(x) : (x))

#define __min(a, b) 		((a) < (b) ? (a) : (b))

#define __max(a, b) 		((a) > (b) ? (a) : (b))

#define __min3(a, b, c)		__min(__min((a), (b)), (c))

#define __max3(a, b, c)		__max(__max((a), (b)), (c))

#define __round_up_32(x) 	(--(x), (x) |= (x) >> 1,		       \
				 (x) |= (x) >> 2, (x) |= (x) >> 4,	       \
				 (x) |= (x) >> 8, (x) |= (x) >> 16, ++(x))

/*
 * Built-in macros function
 */

#define __ALLOC(ptr, sz)	(ptr) = xmalloc(sizeof(*(ptr)) * (sz))

#define __REALLOC(ptr, sz)	(ptr) = xrealloc((ptr), sizeof(*(ptr)) * (sz))

/* push back val to ptr, ptr has sz element, realloc + 1 */
#define __PUSH_BACK(ptr, sz, val) do {					       \
	assert((sz) >= 0);						       \
	__REALLOC((ptr), (sz) + 1);					       \
	(ptr)[(sz)++] = (val);						       \
} while(0)

#define __FREE_AND_NULL(ptr) do {					       \
	free(p);							       \
	(p) = NULL;							       \
} while (0)

#define __VERBOSE(fmt, args...) fprintf(stderr, fmt, ##args)

#define __DEBUG(fmt, x)		fprintf(stderr, "%s = " fmt "\n", #x, x)

#define __DEBUG_ARR(fmt, x, sz) do {					       \
	int iter;							       \
	for (iter = 0; iter < sz; ++iter)				       \
		fprintf(stderr, "%s[%d] = " fmt "\n", #x, iter, x[iter]);      \
} while (0)

#define __WARNING(fmt, args...) do {					       \
	fprintf(stderr, "WARNING: " fmt, ##args);			       \
} while(0)

#define __ERROR(fmt, args...) do {					       \
	fprintf(stderr, "ERROR: " fmt, ##args);				       \
	exit(EXIT_FAILURE);						       \
} while(0)

#define __PERROR(fmt) do {						       \
	fprintf(stderr, "ERROR: ");					       \
	perror(fmt);							       \
	exit(EXIT_FAILURE);						       \
} while(0)

#define __SWAP(x, y) do {						       \
	assert(sizeof(x) == sizeof(y));					       \
	int8_t __swap_temp[sizeof(x)];					       \
	memcpy(__swap_temp, &(y), sizeof(x));				       \
	memcpy(&(y), &(x), sizeof(x));					       \
	memcpy(&(x), __swap_temp, sizeof(x));				       \
} while (0)

/*
 * Built-in function
 */

/* check fread function read enough nmemb */
size_t xfread(void *ptr, size_t size, size_t nmemb, FILE *stream);

/* check fwrite function write enough nmemb */
size_t xfwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);

/* auto remove /n character if found */
ssize_t xgetline(char **str, size_t *n, FILE *stream);

/* get time */
double realtime();

/* make directory if is not exist */
void make_dir(char *path);

/* reverse string */
char *reverse_str(char *str);

/* reverse compelemnt */
char *reverse_complement(char *str);

/*
 * Global variable
 */

extern int8_t nt4_table[256];
extern char nt4_char[5];

#endif /* _UTILS_H_ */