
TARGETS = constant_0.1.txt

all : $(TARGETS)

constant_0.1.txt :
	../mksynthetic -o $@ \
	-m Uniform \
	-I constant_image.img \
	-T constant_true.txt \
	-n -10.0 -N 10.0 \
	-a -10.0 -A 10.0 \
	-f 0.05 -F 0.5 \
	-x 8 -y 8 -z 1 \
	-C \
	-W 2 -H 1 \
	-s 1.0

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

test : constant_0.1.txt prior.txt
	../sliceoptimizer -i constant_0.1.txt -M prior.txt \
	-w 4 -x 8 -y 8 \
	-n -10.0 -N 10.0 \
	-a -10.0 -A 10.0 \
	-s 0 \
	-E 

clean :
	rm -f $(TARGETS)
