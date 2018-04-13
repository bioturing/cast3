<img src="static/cast3_logo.png" width="300" title="CAST3">

# Cloud Alignment Simulation Tool

CAST3 simulates whole genome sequencing data in read-cloud techonology and evaluates the alignment result of such data.
We attempt to simulate real-life genomes by using public SNPs and Structural variants data from 1000 genomes consortium.
Insertion sequences are randomly chosen from retrotransposon database for GRCh37 version.

In version 0.9.0, we only provide the genomic variants of the sample NA12878. You can configure the number of read-pairs, 
number of barcodes, and read-length.

## Installation

```
git clone https://github.com/bioturing/cast3
./build.sh
```

## Usage

### Haplotype simulation

```
Usage:

./bin/sim ./test_data/SV/NA12878.hap.hetA.SV.tsv grch37_reference.fa ./test_data/retro/retros.bed ./test_data/SNP/NA12878.hap.hetA.SNP.tsv hapA.bed > hapA.fasta
./bin/sim ./test_data/SV/NA12878.hap.hetB.SV.tsv grch37_reference.fa ./test_data/retro/retros.bed ./test_data/SNP/NA12878.hap.hetB.SNP.tsv hapB.bed > hapB.fasta

```
### Cloud reads simulation
```
Usage:

./gen_read <option> grch37_reference.fai hapA.bed hapA.fa hapB.bed hapB.fa

option:
  -n INT           # million reads pairs in total to simulated [400]
  -m INT           average # of molecule per barcode [10]
  -l INT           read length [151]
  -o STR           output directory [./]
  -h               print usage and exit

```

### Alignment validation

Updating...

## Output

CAST3 produces two fastq files that include the truth positions of each read in the read names
. User can also investigate the structural variants that spanned by each read by the encoded read name.

Example: 

@4:56907618:135:NOR|4:56907964:151:NOR/1: 

Read 1 from a read pair that span a normal region, which does not include any structural variant.
This read should be aligned at the position 56907618, chromosome 4

