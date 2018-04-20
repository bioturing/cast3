#!/bin/bash

mkdir -p bin

cd lib/htslib && autoheader && autoconf && ./configure && make && cd ../../

cd gen_read && make && cd ../

cd sim_hap && make && cd ../
