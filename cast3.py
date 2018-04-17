#! /usr/bin/env python

from argparse import RawTextHelpFormatter

import os, errno, sys, argparse, subprocess

VERSION = '0.9.0'

sys.stdout = sys.stderr

def perform_main(args):
    if 'func' in args:
        try:
            args.func(args)
        except Exception as ex:
            print(ex)

# define parse
prog_description = 'Cast3 is a program developed by Bioturing INC to simulate shotgun sequencing data in linked-read technology\n'
prog_description += 'It also attempts to validate the alignment of its simulated data\n'
prog_description += 'Please contact info@bioturing.com if you need further support.\n\n'
prog_description += 'Version: '
prog_description += VERSION + '\n\n'
prog_description += 'Note: To use Cast3, you have to simulate two haplotypes first, then generate the sequences after\'.'

parser = argparse.ArgumentParser(
                    formatter_class=RawTextHelpFormatter,
                    description=prog_description)
subparsers = parser.add_subparsers(title='commands')

# index subcommand
def cast3_sim(args):
    ### Generates outdir
    try:
        os.makedirs(args.outdir)
    except OSError as err:
        if err.errno != errno.EEXIST:
            raise
    ### Check the program is available
    if not os.path.exists("./bin/sim.exe"):
        print "Please install Cast3 first!"
        exit()
        
    #### Simulate the haplotypes
    try:
        print "Generating haplotype A"
        val = subprocess.Popen(['./bin/sim.exe', args.hapSV[0], args.ref, args.retro, \
                args.hapSNP[0], args.outdir.strip("/"), ".".join([args.prefix, "hapA"])])
        output = val.communicate()[0]
        print "Generating haplotype B"
        val = subprocess.Popen(['./bin/sim.exe', args.hapSV[1], args.ref, args.retro, \
                args.hapSNP[1], args.outdir.strip("/"), ".".join([args.prefix, "hapB"])])
        output = val.communicate()[0]
    except Exception:
        raise
    if val.returncode != 0:
        exit()

parser_sim = subparsers.add_parser(
                	'sim',
                	formatter_class=RawTextHelpFormatter,
                	description= 'Simulate one haplotype\n\n'
                		         'Example: python cast3.py sim --hapSV ./test_data/SV/NA12878.hap.hetA.SV.tsv ./test_data/SV/NA12878.hap.hetB.SV.tsv --ref\n'
                		         '         Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa --retro ./test_data/retro/retros.bed\n'
                		         '--hapSNP ./test_data/SNP/NA12878.hap.hetA.SNP.tsv ./test_data/SNP/NA12878.hap.hetB.SNP.tsv --outdir sim',
                	help='builds a mermaid index')

parser_sim.add_argument('--hapSV', help='Structural variants tsv file, please see the example in test_data folder', required=True, metavar='STR', nargs = '+')
parser_sim.add_argument('--ref',  help='Path to genome file, GRCh37 version', required=True, metavar='STR')
parser_sim.add_argument('--retro', help='Retrotransposon tsv file, please see the example in test_data folder', required=True, metavar='STR')
parser_sim.add_argument('--hapSNP', help='Single Nucleotide tsv file, please see the example in test_data folder', required=True, metavar='STR', nargs = '+')
parser_sim.add_argument('--outdir', help='Output folder', required=True, metavar='STR')
parser_sim.add_argument('--prefix', help='Prefix of output file', required=False, metavar='STR', default = "out")

parser_sim.set_defaults(func=cast3_sim)

# fusion subcommand
def cast3_gen(args):
    ### Check the program is available
    if not os.path.exists("./bin/gen_read"):
        print "Please install Cast3 first!"
        exit()

    cmd = ["./bin/gen_read", "-n %d" %args.mreads, "-m %d"%args.nmols, \
            "-l %d" %args.length, "-o %s" %args.outdir, "-b %s" %args.barcode]
    cmd.append(args.faidx)
    prefix = args.prefix
    for f in [".hapA.bed", ".hapA.fasta", ".hapB.bed", ".hapB.fasta"]:
        cmd.append("%s/%s%s" %(args.outdir, args.prefix, f))

    print " ".join(cmd)
    #### Generate the sequences
    try:
        val = subprocess.Popen(cmd)
        output = val.communicate()[0]
    except Exception:
        raise
    if val.returncode != 0:
        exit()

parser_gen = subparsers.add_parser(
                    'gen_read',
                	formatter_class=RawTextHelpFormatter,
                	description= 'Generate short read sequencing data in 10X Genomics technology.\n\n'
                		         'Example: python cast3.py gen_read --faidx ../longread/reference/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa.fai --prefix tan --outdir test --nmols 10 --length 151 --mreads 400 --barcode test_data/barcode/4M-with-alts-february-2016.txt',
                                         help='Simulate the pair-end sequencing data in 10X Genomics technology with mean molecule length is 50kb\n'
                                         'and median insert size is 300. Barcode is the 16 base pairs at the beginning at read 1',
                    add_help=False)

parser_gen.add_argument('--faidx', help='Fasta index file of human genome grch37 verison', required=True, metavar='STR')
parser_gen.add_argument('--prefix', help='Prefix of the output file, must be the same as the haplotype sequences', required=True, metavar='STR')
parser_gen.add_argument('--outdir', help='Output folder, which contains the haplotype sequences in the simulation step', required=True, metavar='STR')
parser_gen.add_argument('--barcode', help='File that contains barcode database', required=True, metavar='STR')
parser_gen.add_argument('--nmols', help='Average number of molecules per barcode', required=False, metavar='N', default = 10, type = int)
parser_gen.add_argument('--length', help='Read length', required=False, metavar='N', default = 151, type = int)
parser_gen.add_argument('--mreads', help='Number of million reads pairs in total to be simulated', required=False, metavar='N', default = 400, type = int)


parser_gen.set_defaults(func=cast3_gen)

# version subcommand
def cast3_version(args):
    print('')
    print('Cast3 version %s' % (VERSION))
    print('')
    print('Changes from previous version:')

parser_version = subparsers.add_parser(
                    'version',
                    help='show version information')

parser_version.set_defaults(func=cast3_version)

#Parse args
args = parser.parse_args()
perform_main(args)
