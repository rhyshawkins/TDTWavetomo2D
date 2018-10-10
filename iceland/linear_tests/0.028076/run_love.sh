#!/bin/bash

mkdir -p results_love

../../../wavetomo2dfrequencysliceinvert -i \
	~/PhD/PhDS/tools/ANTEstimate/Transcale/SingleObservations/singleobs_love_0.028076.txt \
	-n -25.60 -N -12.40 \
	-a 61.90 -A 67.90 \
	-M prior.txt \
	-o results_love/ \
	-x 6 -y 5 \
	-w 4 \
	-t 1000000 \
	-v 10000 \
	-L 0.2 \
	-l 200.0 -H 0.1 \
	-k 200 \
	-E

../../../postprocess_slice_mean -i results_love/ch.dat \
	-Z results_love/zoffset.txt \
	-x 6 -y 5 -w 4 \
	-s 1500000 -t 30 \
	-o results_love/mean.txt



