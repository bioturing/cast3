#include "readbam.h"

#define DIFF_ALLOW		1000

int64_t ncorrect, nread;
FILE *f_fail;

struct record {
	char *ref;
	int pos;
	int len;
	char *type;
};

struct record parse_record(char *s, int *pos, int *end)
{
	int u = *pos, v = *pos;
	struct record ret;
	while (s[v] != ':')
		++v;
	s[v] = '\0';
	ret.ref = s + u;
	u = ++v;
	while (s[v] != ':')
		++v;
	s[v] = '\0';
	ret.pos = atoi(s + u);
	u = ++v;
	while (s[v] != ':')
		++v;
	s[v] = '\0';
	ret.len = atoi(s + u);
	u = ++v;
	while (s[v] != '_' && s[v] != '|')
		++v;
	if (s[v] == '|')
		*end = 1;
	else
		*end = 0;
	s[v] = '\0';
	ret.type = s + u;
	*pos = ++v;

	return ret;
}

void check_align(bam1_t *b)
{
	struct record record1[100], record2[100];
	int flag, pos, tid, i, i1, i2, spos, end, found;
	char *seq;

	seq = bam1_qname(b);
	char *tmp = strdup(seq);
	flag = b->core.flag;
	pos = b->core.pos;
	tid = b->core.tid;

	i1 = 0, spos = 0, end = 0;
	while (1) {
		record1[i1++] = parse_record(seq, &spos, &end);
		if (end)
			break;
	}
	i2 = 0, end = 0;
	while (1) {
		record2[i2++] = parse_record(seq, &spos, &end);
		if (end)
			break;
	}

	/* first read */
	if (flag & (1 << 6)) {
		found = 0;
		for (i = 0; i < i1; ++i) {
			if (get_tid(record1[i].ref) == tid &&
			    abs(record1[i].pos - pos) < DIFF_ALLOW) {
				found = 1;
				break;
			}
		}

		if (found) {
			++ncorrect;
		} else {
			fprintf(f_fail, "%s\n", tmp);
			sam_write(f_fail, b);
			fprintf(f_fail, "--------------------------------\n");
		}
	/* second read */
	} else {
		found = 0;
		for (i = 0; i < i2; ++i) {
			if (get_tid(record2[i].ref) == tid &&
			    abs(record2[i].pos - pos) < DIFF_ALLOW) {
				found = 1;
				break;
			}
		}

		if (found) {
			++ncorrect;
		} else {
			fprintf(f_fail, "%s\n", tmp);
			sam_write(f_fail, b);
			fprintf(f_fail, "--------------------------------\n");
		}
	}

	free(tmp);
}

static void process(char *filepath)
{
	int rc;

	init_bam(filepath);
	bam1_t *b = bam_init1();

	do {
		if (nread % 123456 == 0)
			fprintf(stderr, "\r%ld", nread);
		rc = read_bam(b);
		if (!rc)		
			break;
		check_align(b);
		++nread;
	} while (1);
	fprintf(stderr, "\r%ld\n", nread);
}

int main(int argc, char *argv[])
{
	char path[4096];
	sprintf(path, "%s.fail.txt", argv[1]);
	f_fail = fopen(path, "w");

	process(argv[1]);

	__verbose("Total reads: %ld\n", nread);
	__verbose("Total correct: %ld\n", ncorrect);
	__verbose("Ratio: %.2f%%\n", 1.0 * ncorrect / nread);

	return 0;
}
