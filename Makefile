
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
	analyseslicemodel.cpp \
	birth.cpp \
	birthslice.cpp \
	death.cpp \
	deathslice.cpp \
	extractslice.cpp \
	extractdispersion.cpp \
	global.cpp \
	globalslice.cpp \
	hierarchical.cpp \
	hierarchicalmodel.cpp \
	hierarchicalprior.cpp \
	hierarchicalpriorslice.cpp \
	hierarchicalslice.cpp \
	linearweights.hpp \
	mksynthetic.cpp \
	postprocess_acceptance.cpp \
	postprocess_coeff_history.cpp \
	postprocess_coeff_marginal.cpp \
	postprocess_khistory.cpp \
	postprocess_likelihood.cpp \
	postprocess_mean.cpp \
	postprocess_mean_mpi.cpp \
	postprocess_slice_mean.cpp \
	postprocess_slice_likelihood.cpp \
	postprocess_slice_mean_mpi.cpp \
	postprocess_validate_likelihood.cpp \
	ptexchange.cpp \
	ptexchangeslice.cpp \
	resample.cpp \
	rng.cpp \
	slicemodellikelihood.cpp \
	value.cpp \
	valueslice.cpp \
	volume.cpp \
	wavetomo2dexception.cpp \
	wavetomo2dfrequencyinvert.cpp \
	wavetomo2dfrequencyinvert_pt.cpp \
	wavetomo2dfrequencysliceinvert.cpp \
	wavetomo2dfrequencysliceinvert_pt.cpp \
	wavetomo2dobservations.cpp \
	wavetomo2dutil.cpp \
	birth.hpp \
	birthslice.hpp \
	death.hpp \
	deathslice.hpp \
	global.hpp \
	globalslice.hpp \
	hierarchical.hpp \
	hierarchicalmodel.hpp \
	hierarchicalprior.hpp \
	hierarchicalpriorslice.hpp \
	hierarchicalslice.hpp \
	ptexchange.hpp \
	ptexchangeslice.hpp \
	resample.hpp \
	rng.hpp \
	value.hpp \
	valueslice.hpp \
	volume.hpp \
	wavetomo2dexception.hpp \
	wavetomo2dobservations.hpp \
	wavetomo2dutil.hpp \
	generatemask.cpp \
	generatematrices.cpp


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
	scripts/convertsingleobservations.py

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
	postprocess_coeff_history \
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
	generatematrices

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




