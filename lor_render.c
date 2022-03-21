#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "vector_ops.h"
#include "lor.h"


// reads a line from a lor file. BE WARNED: it does not care about line endings,
// if your file does not always have the same number of columns it will break.
// expects a lor file with: c_x, c_y, c_z, d_x, d_y, d_z, sigma_l, sigma_t
lor* read_lor(FILE* input) {
	double center_x;
	double center_y;
	double center_z;
	double dir_x;
	double dir_y;
	double dir_z;
	double longitudinal;
	double transverse;

	int worked;

	worked = fscanf(input, "%lf", &center_x);
	worked = fscanf(input, "%lf", &center_y);
	worked = fscanf(input, "%lf", &center_z);
	worked = fscanf(input, "%lf", &dir_x);
	worked = fscanf(input, "%lf", &dir_y);
	worked = fscanf(input, "%lf", &dir_z);
	worked = fscanf(input, "%lf", &longitudinal);
	worked = fscanf(input, "%lf", &transverse);

	if (worked == EOF) {
		return NULL;
	}
	
	// make a new event to be passed out
	lor* new = (lor*)malloc(sizeof(lor));
	if (new == NULL) {
		return NULL;
	}
	
	new->center = three_vec(center_x, center_y, center_z);
	vec3* cleanup = three_vec(dir_x, dir_y, dir_z);
	new->dir = vec_norm(cleanup);
	free(cleanup);
	new->long_uncert = longitudinal;
	new->transverse_uncert = transverse;


	return new;
}



















int main(int argc, char const *argv[])
{
	/* code */
	return 0;
}
