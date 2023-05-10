#ifndef truth_assign_h
#define truth_assign_h

#include "vector_ops.h"
#include <stdio.h>

#define ORIGIN_BUFFER 22

typedef unsigned int uint;

typedef struct event_ {
  uint number;
  double energy;
  double deposited;
  vec3d location;
  double tof;
  int particle;
  char orgin[ORIGIN_BUFFER];
  int id;
} event;

typedef struct scatter_truth_ {
	int true_n;
	double true_eng;
	double true_time;
} scatter_truth;

typedef struct scatter_ {
	vec3d loc;
	vec3d dir;
	char has_dir;
	double deposit;
	double eng_uncert;
	double space_uncert;
	double time;
	double time_uncert;
	scatter_truth* truth;
} scatter;


double energy_cutoff;

double first_scat_hypot;
double second_scat_hypot;

int predicted_vs_real[4] = {0,0,0,0};



#endif