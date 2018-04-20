#include "utils.h"

int8_t ascii_table[256] = {
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 0, 4, 1,   4, 4, 4, 2,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   3, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 0, 4, 1,   4, 4, 4, 2,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   3, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4
};

char dna_char[5] = "ACGTN";

size_t xfread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fread(ptr, size, nmemb, stream);
	if (ret != nmemb)
		__error("fread error, file is wrong format or is corrupted");
	return ret;
}

size_t xfwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fwrite(ptr, size, nmemb, stream);
	if (ret != nmemb)
		__error("fwrite error, could not write data to file");
	return ret;
}

char *xfgets(char *str, size_t n, FILE *stream)
{
	char *ret = fgets(str, n, stream);
	if (!ret)
		return NULL;
	char *s = str;
	char found = 0;

	for (; *s; ++s) {
		if (*s == '\n') {
			*s = '\0';
			found = 1;
			break;
		}
	}

	if (!found)
		if (fgets(str, n, stream))
			__error("Number character of line exceed limit!");

	return ret;
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

void make_outdir(char *path)
{
	struct stat st = {0};
	if (stat(path, &st) == -1) {
		if (mkdir(path, 0700)) {
			perror("Could not make output directory");
			exit(EXIT_FAILURE);
		}
	}
}

char *reverse_str(char *seq)
{
	int i, j;
	int len = strlen(seq);
	char *ret = malloc(len + 1);
	for (i = 0, j = len - 1; i < len; ++i, --j)
		ret[i] = seq[j];
	ret[len] = '\0';
	return ret;
}

char *reverse_complement(char *seq)
{
	int i, j;
	int len = strlen(seq);
	char *ret = malloc(len + 1);
	for (i = 0, j = len - 1; i < len; ++i, --j) {
		int8_t c = ascii_table[(int)seq[j]];
		if (c == 4)
			ret[i] = dna_char[4];
		else
			ret[i] = dna_char[3 - c];
	}
	ret[len] = '\0';
	return ret;
}