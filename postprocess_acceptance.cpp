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
  
#include "cdf97_lift.h"
#include "cdf97_lift_periodic.h"
#include "haar_lift.h"
#include "daub4_dwt.h"
#include "daub6_dwt.h"
#include "daub8_dwt.h"
#include "generic_lift.h"
  
};

constexpr int NPROPOSALS = 8;

struct user_data {
  int counter;

  int interval;
  int max;

  int accept[NPROPOSALS];
  int reject[NPROPOSALS];
  
  int total_accept[NPROPOSALS];
  int total_reject[NPROPOSALS];
  
};

static const double CREDIBLE_INTERVAL = 0.95;

static int process(int i,
		   void *user,
		   const chain_history_change_t *step,
		   const multiset_int_double_t *S_v);

static void reset(int *ar);
static double safear(int accept, int reject);

static char short_options[] = "i:I:m:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},

  {"interval", required_argument, 0, 'I'},

  {"max", required_argument, 0, 'm'},

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
  int interval;
  int max;
  
  FILE *fp_in;

  struct user_data data;
  multiset_int_double_t *S_v;

  int maxsteps;

  /*
   * Default values
   */
  fp_in = NULL;

  input_file = NULL;
  interval = 1000;
  
  maxsteps = 1000000;

  while (1) {
    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch(c) {

    case 'i':
      input_file = optarg;
      break;

    case 'I':
      interval = atoi(optarg);
      if (interval < 1) {
	fprintf(stderr, "error: interval must be a postive integer\n");
	return -1;
      }
      break;

    case 'm':
      max = atoi(optarg);
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

  ch = chain_history_create(maxsteps);
  if (ch == NULL) {
    fprintf(stderr, "error: failed to create chain history\n");
    return -1;
  }

  data.interval = interval;
  data.max = max;
  data.counter = 0;

  reset(data.total_accept);
  reset(data.total_reject);
  reset(data.accept);
  reset(data.reject);

  S_v = multiset_int_double_create();
  if (S_v == NULL) {
    fprintf(stderr, "error: failed to create multiset\n");
    return -1;
  }
  
  fp_in = fopen(input_file, "r");
  if (fp_in == NULL) {
    fprintf(stderr, "error: failed to open input file\n");
    return -1;
  }

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
  fclose(fp_in);
  printf("Total %d\n  ", data.counter);
  for (int i = 0; i < NPROPOSALS; i ++) {
    printf("%6d ", data.total_accept[i]);
  }
  printf("\n  ");
  for (int i = 0; i < NPROPOSALS; i ++) {
    printf("%6d ", data.total_accept[i] + data.total_reject[i]);
  }
  printf("\n  ");
  for (int i = 0; i < NPROPOSALS; i ++) {
    printf("%6.2f ", safear(data.total_accept[i], data.total_reject[i]));
  }
  printf("\n");

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
  
  if (d->counter < d->max) {
    if (step->header.accepted) {
      d->accept[(int)step->header.type] ++;
      d->total_accept[(int)step->header.type] ++;
    } else {
      d->reject[(int)step->header.type] ++;
      d->total_reject[(int)step->header.type] ++;
    }
    
    d->counter ++;
    
    if (d->counter % d->interval == 0) {
      printf("%d\n  ", d->counter);
      for (int i = 0; i < NPROPOSALS; i ++) {
	printf("%6d ", d->accept[i]);
      }
      printf("\n  ");
      for (int i = 0; i < NPROPOSALS; i ++) {
	printf("%6d ", d->accept[i] + d->reject[i]);
      }
      printf("\n  ");
      for (int i = 0; i < NPROPOSALS; i ++) {
	printf("%6.2f ", safear(d->accept[i], d->reject[i]));
      }
      printf("\n");
      
      reset(d->accept);
      reset(d->reject);
    }
  }
  
  return 0;
}

static void reset(int *ar)
{
  for (int i = 0; i < NPROPOSALS; i ++) {
    ar[i] = 0;
  }
}

static double safear(int accept, int reject)
{
  int total = accept + reject;
  if (total == 0) {
    return 0.0;
  } else {
    return 100.0 * (double)accept/(double)total;
  }
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -O|--observations <filename>     Input observations\n"
	  " -S|--stm <filename>              Input stm file\n"
	  " -H|--hierarchical <filename>     Input hierarchical file\n"
	  "\n"
	  " -d|--degree-depth <int>    Number of layers as power of 2\n"
	  " -l|--degree-lateral <int>  Number of horizontal samples as power of 2\n"
	  "\n"
	  " -i|--input <file>                Input ch file\n"
	  "\n"
	  " -t|--thin <int>                  Only processing every ith sample\n"
	  " -s|--skip <int>                  Skip n samples from beginning\n"
	  "\n"
	  " -w|--wavelet-vertical <int>      Wavelet for vertical direction\n"
	  " -W|--wavelet-horizontal <int>    Wavelet for horizontal direction\n"
	  "\n"
	  " -h|--help            Show usage\n"
	  "\n",
	  pname);
}

 
