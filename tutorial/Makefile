
TARGETS = tutorial_raw_0.1.txt \
	tutorial_stations.txt \
	tutorial_observations.txt \
	tutorial_data.txt

all : $(TARGETS)

tutorial_raw_0.1.txt :
	../mksynthetic -o $@ \
	-m CosineCheck0 \
	-I tutorial_image.img \
	-n -10.0 -N 10.0 \
	-a -10.0 -A 10.0 \
	-f 0.05 -F 0.5 \
	-x 8 -y 8 -z 1 \
	-C \
	-W 4 -H 4 \
	-s 0.1

tutorial_stations.txt : tutorial_raw_0.1.txt
	python2 extract_stations.py -i $< -o $@

tutorial_observations.txt : tutorial_raw_0.1.txt
	python2 extract_observations.py -i $< -o $@ -s 0

tutorial_data.txt : tutorial_stations.txt tutorial_observations.txt
	python2 ../scripts/convertsingleobservations.py \
	-s tutorial_stations.txt \
	-v tutorial_observations.txt \
	-f 1.0 \
	-o $@


clean :
	rm -f $(TARGETS)
