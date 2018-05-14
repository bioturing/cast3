#include "utils.h"

int8_t nt4_table[256] = {
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

char nt4_char[5] = "ACGTN";

size_t xfread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fread(ptr, size, nmemb, stream);
	if (ret != nmemb)
		__ERROR("fread error, file is wrong format or is corrupted");
	return ret;
}

size_t xfwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fwrite(ptr, size, nmemb, stream);
	if (ret != nmemb)
		__ERROR("fwrite error, could not write data to file");
	return ret;
}

ssize_t xgetline(char **str, size_t *n, FILE *stream)
{
	ssize_t ret = getline(str, n, stream);
	if (ret == 0 || ret == -1)
		return ret;
	if ((*str)[ret - 1] == '\n')
		(*str)[--ret] = '\0';
	return ret;
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

void make_dir(char *path)
{
	struct stat st = {0};
	if (stat(path, &st) == -1)
		if (mkdir(path, 0700))
			__PERROR("Could not make output directory");
}

char *reverse_str(char *str)
{
	int i, j, len = strlen(str);
	char *ret = malloc(len + 1);
	for (i = 0, j = len - 1; i < len; ++i, --j)
		ret[i] = str[j];
	ret[len] = '\0';
	return ret;
}

char *reverse_complement(char *str)
{
	int i, j, len = strlen(str);
	char *ret = malloc(len + 1);
	for (i = 0, j = len - 1; i < len; ++i, --j) {
		int8_t c = nt4_table[(int)str[j]];
		if (c == 4)
			ret[i] = nt4_char[4];
		else
			ret[i] = nt4_char[3 - c];
	}
	ret[len] = '\0';
	return ret;
}