#!/bin/bash

mkdir -p bin

cd gen_read && make && cd ../
cd sim_hap && make && cd ../
