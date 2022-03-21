#ifndef truth_assign_h
#define truth_assign_h

#include "vector_ops.h"
#include <stdio.h>

typedef unsigned int uint;

typedef struct event_ {
  uint number;
  double energy;
  double depoisted;
  vec3* location;
  double tof;
  int particle;
  char orgin[22];
  int id;
} event;

typedef struct scatter_ {
	vec3* loc;
	vec3* dir;
	double deposit;
	double eng_uncert;
	double space_uncert;
	double time;
	double time_uncert;
} scatter;

double energy_cutoff;

double first_scat_hypot;
double second_scat_hypot;

int predicted_vs_real[4] = {0,0,0,0};



#endif