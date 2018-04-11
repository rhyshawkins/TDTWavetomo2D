//
//    TDTWavetomo2d : Software for the inversion of surface wave datasets using the
//    trans-dimensional tree approach using a wavelet parameterisation. See
//
//      R Hawkins and M Sambridge, "Geophysical imaging using trans-dimensional trees",
//      Geophysical Journal International, 2015, 203:2, 972 - 1000,
//      https://doi.org/10.1093/gji/ggv326
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <getopt.h>

extern "C" {
#include "chain_history.h"
  
};

#include "global.hpp"

#include "wavetomo2dutil.hpp"
#include "wavetomo2dexception.hpp"

#include <map>

struct cid {
  int depth;
  int index;

  cid(int d, int i) :
    depth(d),
    index(i)
  {
  }

  friend bool operator<(const cid &a, const cid &b)
  {
    if (a.depth == b.depth) {
      return a.index < b.index;
    } else {
      return a.depth < b.depth;
    }
  }
};

struct coefficient_counter {

  coefficient_counter(int _depth, int _index) :
    idx(_depth, _index),
    pv(0),
    av(0),
    pb(0),
    ab(0),
    pd(0),
    ad(0),
    meann(0),
    mean(0.0),
    var(0.0)
  {
  }

  void update_mean(double v)
  {
    meann ++;
    double delta = v - mean;
    mean += delta/(double)(meann);
    var += delta * (v - mean);
  }

  double safepercent(int p, int a)
  {
    if (p == 0) {
      return 0.0;
    }
    return 100.0 * (double)a/(double)p;
  }

  double safesigma()
  {
    if (meann < 2) {
      return 0.0;
    } else {
      return sqrt(var/(double)(meann - 1));
    }
  }
  
  void write(FILE *fp)
  {
    fprintf(fp, "%3d %5d B %5d %7.3f D %5d %7.3f V %5d %7.3f mu %10.6f sigma %10.6f N %5d\n",
	    idx.depth, idx.index,
	    pb, safepercent(pb, ab),
	    pd, safepercent(pd, ad),
	    pv, safepercent(pv, av),
	    mean,
	    safesigma(),
	    meann);
  }

  cid idx;
  
  int pv;
  int av;

  int pb;
  int ab;

  int pd;
  int ad;

  int meann;
  double mean;
  double var;
};

  

struct user_data {
  int thincounter;
  int thin;
  int skip;
  
  int counter;

  coefficient_counter *get_or_create(int depth, int idx)
  {
    std::map<cid, coefficient_counter*>::iterator i;
    cid k(depth, idx);
    
    i = coefficients.find(k);
    if (i == coefficients.end()) {
      i = coefficients.insert(std::pair<cid, coefficient_counter*>(k, new coefficient_counter(k.depth, k.index))).first;
    }

    return i->second;
  }

  bool save(const char *filename)
  {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
      return false;
    }
    
    for (auto &i : coefficients) {
      i.second->write(fp);
    }
  }
  
  std::map<cid, coefficient_counter*> coefficients;
};

static int process(int i,
		   void *user,
		   const chain_history_change_t *step,
		   const multiset_int_double_t *S_v);

static char short_options[] = "i:o:t:s:S:h";
static struct option long_options[] = {

  {"input", required_argument, 0, 'i'},
  {"output", required_argument, 0, 'o'},

  {"thin", required_argument, 0, 't'},
  {"skip", required_argument, 0, 's'},

  {"maxsteps", required_argument, 0, 'S'},

  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;
  
  chain_history_t *ch;
  
  std::vector<std::string> input_file;
  
  char *output_file;

  int thin;
  int skip;
  int maxsteps;

  FILE *fp_in;
  FILE *fp_out;

  struct user_data data;
  multiset_int_double_t *S_v;

  int i;
  int j;

  int processesperchain;
  
  /*
   * Default values
   */
  fp_in = NULL;
  fp_out = NULL;
  
  output_file = NULL;

  thin = 0;
  skip = 0;

  maxsteps = 1000000;

  processesperchain = 1;
  
  while (1) {
    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch(c) {
    case 'i':
      input_file.push_back(optarg);
      break;

    case 'o':
      output_file = optarg;
      break;

    case 't':
      thin = atoi(optarg);
      break;

    case 's':
      skip = atoi(optarg);
      break;
      
    case 'S':
      maxsteps = atoi(optarg);
      if (maxsteps < 1000) {
	fprintf(stderr, "error: maxsteps should be 1000 or greater\n");
	return -1;
      }
      break;

    case 'h':
    default:
      usage(argv[0]);
      return -1;
      
    }
  }

  if (input_file.size() == 0) {
    fprintf(stderr, "error: required parameter input file missing\n");
    return -1;
  }

  if (output_file == NULL) {
    fprintf(stderr, "error: required parameter output file missing\n");
    return -1;
  }

  ch = chain_history_create(maxsteps);
  if (ch == NULL) {
    fprintf(stderr, "error: failed to create chain history\n");
    return -1;
  }
  
  data.thincounter = 0;
  data.thin = thin;
  data.skip = skip;

  data.counter = 0;

  S_v = multiset_int_double_create();
  if (S_v == NULL) {
    fprintf(stderr, "error: failed to create multiset\n");
    return -1;
  }

  for (auto &infile: input_file) {

    fp_in = fopen(infile.c_str(), "r");
    if (fp_in == NULL) {
      fprintf(stderr, "error: failed to open input file\n");
      return -1;
    }
    printf("Loaded: %s\n", infile.c_str());
    
    /*
     * Process the chain history
     */
    while (!feof(fp_in)) {
      
      if (chain_history_read(ch,
			     (ch_read_t)fread,
			     fp_in) < 0) {
	if (feof(fp_in)) {
	  break;
	}
	
	fprintf(stderr, "error: failed to read chain history\n");
	return -1;
      }
      
      if (chain_history_replay(ch,
			       S_v,
			       (chain_history_replay_function_t)process,
			       &data) < 0) {
	fprintf(stderr, "error: failed to replay\n");
	return -1;
      }
    }
    printf("%d records\n", data.counter);
    fclose(fp_in);
  }

  if (!data.save(output_file)) {
    fprintf(stderr, "error: failed to save history\n");
    return -1;
  }
  
  chain_history_destroy(ch);
  multiset_int_double_destroy(S_v);

  return 0;
}

static int process(int stepi,
		   void *user,
		   const chain_history_change_t *step,
		   const multiset_int_double_t *S_v)
{
  struct user_data *d = (struct user_data *)user;

  coefficient_counter *cc;
  double delta;
  
  if ((d->thincounter >= d->skip) && (d->thin <= 1 || (d->thincounter % d->thin) == 0)) {

    switch (step->header.type) {
    case CH_NOCHANGE:
      printf("No change ignored\n");
      break;
      
    case CH_MOVE:
      printf("Move ignored\n");
      break;

    case CH_PTEXCHANGE:
      printf("PT Exchange ignored\n");
      break;

    case CH_HIERARCHICAL:
      break;

    case CH_INITIALISE:
      {
	int depth = 0;
	int count = multiset_int_double_depth_count(S_v, depth);
	while (count > 0) {
	  for (int i = 0; i < count; i ++) {

	    int idx;
	    double value;
	      
	    if (multiset_int_double_nth_element(S_v, depth, i, &idx, &value) < 0) {
	      fprintf(stderr, "failed to get nth element");
	      return -1;
	    }

	    cc = d->get_or_create(depth, idx);
	    cc->update_mean(value);
	    
	  }
	  
	  depth ++;
	  count = multiset_int_double_depth_count(S_v, depth);
	}
      }
      break;

    case CH_BIRTH:
      
      cc = d->get_or_create(step->perturbation.birth.node_depth,
			    step->perturbation.birth.node_id);

      cc->pb ++;
      if (step->header.accepted) {
	cc->ab ++;
	cc->update_mean(step->perturbation.birth.new_value);
      }
      break;
      
    case CH_DEATH:
      cc = d->get_or_create(step->perturbation.death.node_depth,
			    step->perturbation.death.node_id);

      cc->pd ++;
      if (step->header.accepted) {
	cc->ad ++;
      }
      break;

    case CH_VALUE:
      cc = d->get_or_create(step->perturbation.value.node_depth,
			    step->perturbation.value.node_id);
      
      cc->pv ++;
      if (step->header.accepted) {
	cc->av ++;

	cc->update_mean(step->perturbation.value.new_value);
      }
      break;
    }
  }
  d->thincounter ++;
  
  return 0;
}


static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -i|--input <file>                Input ch file\n"
	  " -o|--output <file>               Output mean model file\n"
	  "\n"
	  " -t|--thin <int>                  Only processing every ith sample\n"
	  " -s|--skip <int>                  Skip n samples from beginning\n"
	  "\n"
	  " -S|--maxsteps <int>              Chain history max steps\n"
	  "\n"
	  " -h|--help            Show usage\n"
	  "\n",
	  pname);
}

 
