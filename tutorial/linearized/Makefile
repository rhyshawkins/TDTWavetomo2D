
TARGETS = results_linearized/ch.dat

LONMIN=-10.0
LONMAX=10.0
LATMIN=-10.0
LATMAX=10.0

all : $(TARGETS)

#
# Generic prior
#
prior.txt :
	echo "laplace" > $@
	echo "0.1" >> $@
	echo "priorbirth" >> $@
	echo "depthgaussianperturb" >> $@
	echo "7" >> $@
	echo "0.020" >> $@
	echo "0.020" >> $@
	echo "0.025" >> $@
	echo "0.035" >> $@
	echo "0.050" >> $@
	echo "0.050" >> $@
	echo "0.050" >> $@


#
# Uniform model
#
../synthetic_models/uniform_444_0.05.txt :
	make -C ../synthetic_models uniform_444_0.05.txt

results_linearized/ch.dat : ../tutorial_data.txt prior.txt
	mkdir -p results_linearized
	../../wavetomo2dfrequencysliceinvert -i ../tutorial_data.txt \
	-M prior.txt \
	-o results_linearized/ \
	-s 0 \
	-x 5 -y 5 \
	-u 2 \
	-n $(LONMIN) -N $(LONMAX) \
	-a $(LATMIN) -A $(LATMAX) \
	-w 4 \
	-t 100000 \
	-l 1.0 \
	-k 100 \
	-H 0.1 \
	-E

results_linearized/mean.txt : results_linearized/ch.dat
	../../postprocess_slice_mean -i $< -Z results_linearized/zoffset.txt \
	-o $@ -s 50000 -t 10 \
	-x 5 -y 5 -w 4

results_linearized_hp/ch.dat : ../tutorial_data.txt prior.txt
	mkdir -p results_linearized_hp
	../../wavetomo2dfrequencysliceinvert -i ../tutorial_data.txt \
	-M prior.txt \
	-o results_linearized_hp/ \
	-s 0 \
	-x 5 -y 5 \
	-u 2 \
	-n $(LONMIN) -N $(LONMAX) \
	-a $(LATMIN) -A $(LATMAX) \
	-w 4 \
	-t 100000 \
	-l 1.0 \
	-k 100 \
	-H 0.1 \
	-L 0.2 \
	-E

results_linearized_hp/mean.txt : results_linearized_hp/ch.dat
	../../postprocess_slice_mean -i $< -Z results_linearized_hp/zoffset.txt \
	-o $@ -s 50000 -t 10 \
	-x 5 -y 5 -w 4

results_linearized_hp/like.txt : results_linearized_hp/ch.dat
	../../postprocess_slice_likelihood -i $< -o $@ \
	-H results_linearized_hp/hierarchical.txt \
	-L results_linearized_hp/prior.txt \
	-K results_linearized_hp/khistory.txt \
	-t 100

#
# Parallel version
#
results_linearized_pt/ch.dat : ../tutorial_data.txt prior.txt
	mkdir -p results_linearized_pt
	mpirun -np 8 ../../wavetomo2dfrequencysliceinvert_pt -i ../tutorial_data.txt \
	-M prior.txt \
	-o results_linearized_pt/ \
	-s 0 \
	-x 5 -y 5 \
	-u 2 \
	-n $(LONMIN) -N $(LONMAX) \
	-a $(LATMIN) -A $(LATMAX) \
	-w 4 \
	-t 100000 \
	-l 1.0 \
	-k 100 \
	-H 0.1 \
	-E \
	-c 2
