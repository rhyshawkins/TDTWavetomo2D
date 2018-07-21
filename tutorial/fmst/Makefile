
TARGETS = fmst-sources.dat \
	gridi.vtx

all : $(TARGETS)

fmst-sources.dat : ../tutorial_data.txt
	python2 ../../scripts/converttofmst.py -i $< -o fmst

gridi.vtx :
	python2 ../../scripts/imagetovtx.py -m 3 -o $@ \
	-n -10.0 -N 10.0 \
        -a -10.0 -A 10.0


