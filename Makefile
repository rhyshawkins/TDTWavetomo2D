
TDTBASE=../TDTbase
TRAVELTIMEBASE=..

INCLUDES = \
	-I$(TDTBASE)/log \
	-I$(TDTBASE)/hnk \
	-I$(TDTBASE)/tracking \
	-I$(TDTBASE)/wavetree \
	-I$(TDTBASE)/oset \
	-I$(TDTBASE)/wavelet \
	-I$(TRAVELTIMEBASE)/traveltime2d

EXTRA_LIBS = \
	-L$(TDTBASE)/hnk -lhnk \
	-L$(TDTBASE)/tracking -ltracking \
	-L$(TDTBASE)/log -llog \
	-L$(TDTBASE)/wavetree -lwavetree \
	-L$(TDTBASE)/oset -loset \
	-L$(TDTBASE)/wavelet -lwavelet \
	-L$(TRAVELTIMEBASE)/traveltime2d -ltraveltime2d

CXX ?= mpicxx
CXXFLAGS = -c -g -Wall --std=c++11 

CXXFLAGS += -O3

CXXFLAGS += $(INCLUDES)

INSTALL = install
INSTALLFLAGS = -D

LIBS = $(EXTRA_LIBS) -lm $(shell gsl-config --libs) -lgmp 
#MPI_LIBS = $(shell mpicxx -showme:link)
MPI_LIBS = 

OBJS = wavetomo2dexception.o \
	wavetomo2dutil.o \
	wavetomo2dobservations.o \
	hierarchical.o \
	hierarchicalmodel.o \
	hierarchicalprior.o \
	rng.o \
	global.o \
	birth.o \
	death.o \
	value.o \
	ptexchange.o \
	ptexchangeslice.o \
	resample.o \
	volume.o \
	globalslice.o \
	birthslice.o \
	deathslice.o \
	valueslice.o \
	hierarchicalslice.o \
	hierarchicalpriorslice.o

MPI_OBJS = 

SRCS = Makefile \
	$(wildcard *.?pp)

EXTRADIST = LICENSE \
	README.md \
	CHANGELOG \
	documentation/manual.tex \
	documentation/manual.pdf \
	documentation/bibliography.bib \
	tutorial/tutorial_data.txt \
        tutorial/tutorial_observations.txt \
        tutorial/tutorial_prior.txt \
        tutorial/tutorial_stations.txt \
	tutorial/fmst/fm2dss.in \
	tutorial/fmst/Makefile \
	tutorial/fmst/residualss.in \
	tutorial/fmst/subinvss.in \
	tutorial/fmst/ttomoss.in \
	scripts/convertsingleobservations.py \
	scripts/converttofmst.py \
	scripts/imagetovtx.py \
	scripts/mkstarting.py \
	scripts/vtxtoimage.py \
	slurm/example_parallel.slurm \
	pbs/pbs_example_parallel.sh

TARGETS = wavetomo2dfrequencyinvert \
	wavetomo2dfrequencyinvert_pt \
	wavetomo2dfrequencysliceinvert \
	wavetomo2dfrequencysliceinvert_pt \
	postprocess_coeff_marginal \
	postprocess_likelihood \
	postprocess_mean \
	postprocess_mean_mpi \
	postprocess_validate_likelihood \
	postprocess_acceptance \
	postprocess_khistory \
	mksynthetic \
	extractslice \
	extractdispersion \
	postprocess_slice_mean \
	postprocess_slice_likelihood \
	postprocess_slice_mean_mpi \
	analyseslicemodel \
	slicemodellikelihood \
	generatemask \
	generatematrices \
	sliceoptimizer

all : $(TARGETS)

wavetomo2dfrequencyinvert : wavetomo2dfrequencyinvert.o $(OBJS)
	$(CXX) -o wavetomo2dfrequencyinvert wavetomo2dfrequencyinvert.o $(OBJS) $(LIBS) $(MPI_LIBS)

wavetomo2dfrequencyinvert_pt : wavetomo2dfrequencyinvert_pt.o $(OBJS)
	$(CXX) -o wavetomo2dfrequencyinvert_pt wavetomo2dfrequencyinvert_pt.o $(OBJS) $(LIBS) $(MPI_LIBS)

wavetomo2dfrequencysliceinvert : wavetomo2dfrequencysliceinvert.o $(OBJS)
	$(CXX) -o wavetomo2dfrequencysliceinvert wavetomo2dfrequencysliceinvert.o $(OBJS) $(LIBS) $(MPI_LIBS)

wavetomo2dfrequencysliceinvert_pt : wavetomo2dfrequencysliceinvert_pt.o $(OBJS)
	$(CXX) -o wavetomo2dfrequencysliceinvert_pt wavetomo2dfrequencysliceinvert_pt.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_coeff_marginal : postprocess_coeff_marginal.o $(OBJS)
	$(CXX) -o postprocess_coeff_marginal postprocess_coeff_marginal.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_likelihood : postprocess_likelihood.o $(OBJS)
	$(CXX) -o postprocess_likelihood postprocess_likelihood.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_mean : postprocess_mean.o $(OBJS)
	$(CXX) -o postprocess_mean postprocess_mean.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_mean_mpi : postprocess_mean_mpi.o $(OBJS)
	$(CXX) -o postprocess_mean_mpi postprocess_mean_mpi.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_validate_likelihood : postprocess_validate_likelihood.o $(OBJS)
	$(CXX) -o postprocess_validate_likelihood postprocess_validate_likelihood.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_acceptance : postprocess_acceptance.o $(OBJS)
	$(CXX) -o postprocess_acceptance postprocess_acceptance.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_coeff_history : postprocess_coeff_history.o $(OBJS)
	$(CXX) -o postprocess_coeff_history postprocess_coeff_history.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_khistory : postprocess_khistory.o $(OBJS)
	$(CXX) -o postprocess_khistory postprocess_khistory.o $(OBJS) $(LIBS) $(MPI_LIBS)

mksynthetic : mksynthetic.o $(OBJS)
	$(CXX) -o mksynthetic mksynthetic.o $(OBJS) $(LIBS) $(MPI_LIBS)

extractslice : extractslice.o $(OBJS)
	$(CXX) -o extractslice extractslice.o $(OBJS) $(LIBS) $(MPI_LIBS)

extractdispersion : extractdispersion.o $(OBJS)
	$(CXX) -o extractdispersion extractdispersion.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_slice_mean : postprocess_slice_mean.o $(OBJS)
	$(CXX) -o postprocess_slice_mean postprocess_slice_mean.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_slice_likelihood : postprocess_slice_likelihood.o $(OBJS)
	$(CXX) -o postprocess_slice_likelihood postprocess_slice_likelihood.o $(OBJS) $(LIBS) $(MPI_LIBS)

postprocess_slice_mean_mpi : postprocess_slice_mean_mpi.o $(OBJS)
	$(CXX) -o postprocess_slice_mean_mpi postprocess_slice_mean_mpi.o $(OBJS) $(LIBS) $(MPI_LIBS)

analyseslicemodel : analyseslicemodel.o $(OBJS)
	$(CXX) -o analyseslicemodel analyseslicemodel.o $(OBJS) $(LIBS) $(MPI_LIBS)

slicemodellikelihood : slicemodellikelihood.o $(OBJS)
	$(CXX) -o slicemodellikelihood slicemodellikelihood.o $(OBJS) $(LIBS) $(MPI_LIBS)

generatemask : generatemask.o $(OBJS)
	$(CXX) -o generatemask generatemask.o $(OBJS) $(LIBS) $(MPI_LIBS)

generatematrices : generatematrices.o $(OBJS)
	$(CXX) -o generatematrices generatematrices.o $(OBJS) $(LIBS) $(MPI_LIBS)

sliceoptimizer : sliceoptimizer.o $(OBJS)
	$(CXX) -o sliceoptimizer sliceoptimizer.o $(OBJS) $(LIBS) $(MPI_LIBS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -o $*.o $*.cpp

DATE = $(shell date +"%Y%m%d%H%M")
DIR = TDTWavetomo2D
TGZ = $(DIR).tar.gz

documentation/manual.pdf : documentation/manual.tex
	cd documentation && pdflatex manual.tex && bibtex manual && pdflatex manual.tex && pdflatex manual.tex 

dist : $(EXTRADIST)
	mkdir -p $(DIR)
	echo $(DATE) > $(DIR)/Version
	for f in Makefile $(SRCS) $(EXTRADIST); do \
	    $(INSTALL) $(INSTALLFLAGS) $$f $(DIR)/$$f ; \
	done
	tar -czf $(TGZ) $(DIR)/*
	rm -rf $(DIR)

clean :
	rm -f $(TARGETS) *.o




