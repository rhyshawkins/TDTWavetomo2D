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

#include <getopt.h>

extern "C" {
#include "chain_history.h"
};

static const int CHAIN_MAXSTEPS = 1000000;

struct user_data {
  int thin;
  int thincounter;
  
  int counter;

  FILE *fp_out;
  FILE *fp_hierarchical_out;
  FILE *fp_prior_out;

  double hyper;
};

static int process(int i,
		   void *user,
		   const chain_history_change_t *step,
		   const multiset_int_double_t *S_v);

static char short_options[] = "i:o:H:P:t:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"output", required_argument, 0, 'o'},
  {"hierarchical", required_argument, 0, 'H'},
  {"prior", required_argument, 0, 'P'},
  {"thin", required_argument, 0, 't'},

  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;
  
  chain_history_t *ch;
  
  char *input_file;
  char *output_file;
  char *hierarchical_file;
  char *prior_file;

  FILE *fp_in;
  FILE *fp_out;

  int thin;

  struct user_data data;
  multiset_int_double_t *S_v;

  input_file = NULL;
  output_file = NULL;
  hierarchical_file = NULL;
  prior_file = NULL;
  
  thin = 0;
  
  while (1) {
    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch(c) {
    case 'i':
      input_file = optarg;
      break;

    case 'o':
      output_file = optarg;
      break;

    case 'H':
      hierarchical_file = optarg;
      break;

    case 'P':
      prior_file = optarg;
      break;

    case 't':
      thin = atoi(optarg);
      break;
      
    case 'h':
    default:
      usage(argv[0]);
      return -1;
      
    }
  }

  if (input_file == NULL) {
    fprintf(stderr, "error: required parameter input file missing\n");
    return -1;
  }

  if (output_file == NULL) {
    fprintf(stderr, "error: required parameter output file missing\n");
    return -1;
  }

  ch = chain_history_create(CHAIN_MAXSTEPS);
  if (ch == NULL) {
    fprintf(stderr, "error: failed to create chain history\n");
    return -1;
  }
  
  fp_in = fopen(input_file, "r");
  if (fp_in == NULL) {
    fprintf(stderr, "error: failed to open input file\n");
    return -1;
  }

  fp_out = fopen(output_file, "w");
  if (fp_out == NULL) {
    fprintf(stderr, "error: failed to open output file\n");
    return -1;
  }


  data.thin = thin;
  data.thincounter = 0;
  data.counter = 0;
  data.fp_out = fp_out;

  if (hierarchical_file != NULL) {
    data.fp_hierarchical_out = fopen(hierarchical_file, "w");
    if (data.fp_hierarchical_out == NULL) {
      fprintf(stderr, "error: failed to create hierarchical output file\n");
      return -1;
    }
  } else {
    data.fp_hierarchical_out = NULL;
  }
  
  if (prior_file != NULL) {
    data.fp_prior_out = fopen(prior_file, "w");
    if (data.fp_prior_out == NULL) {
      fprintf(stderr, "error: failed to create prior output file\n");
      return -1;
    }
  } else {
    data.fp_prior_out = NULL;
  }
  data.hyper = 1.0;

  S_v = multiset_int_double_create();
  if (S_v == NULL) {
    fprintf(stderr, "error: failed to create multiset\n");
    return -1;
  }
  
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
  fclose(fp_out);

  if (hierarchical_file != NULL) {
    fclose(data.fp_hierarchical_out);
  }
  
  return 0;
}

static int process(int i,
		   void *user,
		   const chain_history_change_t *step,
		   const multiset_int_double_t *S_v)
{
  struct user_data *d = (struct user_data *)user;

  if (step->header.type == CH_HYPER) {
    if (step->header.accepted) {
      d->hyper = step->perturbation.hyper.new_value;
    }
  }

  if (d->thin <= 1 || (d->thincounter % d->thin == 0)) {
    fprintf(d->fp_out, "%.6f\n", step->header.likelihood);

    if (d->fp_hierarchical_out != NULL) {
      fprintf(d->fp_hierarchical_out, "%.6f\n", step->header.hierarchical);
    }
      
    if (d->fp_prior_out != NULL) {
      fprintf(d->fp_prior_out, "%10.6f\n", d->hyper);
    }
      
    d->counter ++;
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
	  " -i|--input <file>        Input chain history file\n"
	  " -o|--output <file>       Output likelihood file\n"
	  " -H|--hierarchical <file> Output hierarchical parameter to file\n"
	  " -P|--prior <file>        Output hierarchical prior to file\n"
	  " -t|--thin <int>          Thinning\n"
	  "\n"
	  " -h|--help              Usage information\n"
	  "\n",
	  pname);
}

