#include "generate.h"

#include <unistd.h>
#include <sys/stat.h>

#define BARCODE_PATH                "../4M-with-alts-february-2016.txt"
#define OUT_DIR                     "./"

struct argument_t {
    std::string outdir;
    std::string barcode_path; 
    int64_t total_read_pair;
    int mean_mlc_per_bx;
    int read_len;
    char *hap1_path;
    char *hap2_path;
    char *bed1_path;
    char *bed2_path;
    char *fai_path;
} args;

void print_usage()
{
    fprintf(stderr, "./gen_read <option> ref.fai hap1.bed hap1.fa hap2.bed hap2.fa\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "option:\n");
    fprintf(stderr, "  -n INT           # million reads pairs in total to simulated [400]\n");
    fprintf(stderr, "  -m INT           average # of molecule per barcode [10]\n");
    fprintf(stderr, "  -l INT           read length [151]\n");
    fprintf(stderr, "  -o STR           output directory [./]\n");
    fprintf(stderr, "  -b STR           file that contains barcode database [./]\n");
    fprintf(stderr, "  -h               print usage and exit\n");
}

void make_outdir(const char *path)
{
    struct stat st = {0};
    if (stat(path, &st) == -1) {
        if (mkdir(path, 0700)) {
            perror("Could not make output directory");
            exit(EXIT_FAILURE);
        }
    }
}

void parse_argument(int argc, char *argv[])
{
    args.outdir = OUT_DIR;
    args.barcode_path = BARCODE_PATH;
    args.mean_mlc_per_bx = 10;
    args.total_read_pair = 400000000;
    args.read_len = 151;
    int c;

    while ((c = getopt(argc, argv, "n:m:l:ho:")) != -1) {
        switch (c) {
        case 'n':
            args.total_read_pair = atoi(optarg) * 1000000;
            break;
        case 'm':
            args.mean_mlc_per_bx = atoi(optarg);
            break;
        case 'o':
            args.outdir = std::string(optarg);
            break;
        case 'l':
            args.read_len = atoi(optarg);
            break;
        case 'b':
        	args.barcode_path = std::string(optarg);
        	break;
        case 'h':
            print_usage();
            exit(0);
        default:
            fprintf(stderr, "Argument error: unknow option\n\n");
            print_usage();
            exit(1);
        }
    }

    if (optind + 5 == argc) {
        args.fai_path = argv[optind];
        args.bed1_path = argv[optind + 1];
        args.hap1_path = argv[optind + 2];
        args.bed2_path = argv[optind + 3];
        args.hap2_path = argv[optind + 4];
    } else {
        fprintf(stderr, "Argument error: wrong format\n\n");
        print_usage();
        exit(1);
    }

    make_outdir(args.outdir.c_str());
}

int main(int argc, char *argv[])
{
    parse_argument(argc, argv);
    generate_t generate(args.barcode_path.c_str(), args.hap1_path, args.bed1_path,
                        args.hap2_path, args.bed2_path, args.fai_path);
    generate.generate_read(args.read_len, args.total_read_pair, args.mean_mlc_per_bx);
    return 0;
}
