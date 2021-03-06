#!/bin/bash

git submodule init && git submodule update

cd htslib && autoheader && autoconf && ./configure && make && cd ../

mkdir -p bin

cd gen_read && make && cd ../

cd sim_hap && make && cd ../

cd validate_align && make && cd ../

cd bam_stats && make && cd ../
