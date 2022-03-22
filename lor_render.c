#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "vector_ops.h"
#include "lor.h"


// reads a line from a lor file. BE WARNED: it does not care about line endings,
// if your file does not always have the same number of columns it will break.
// expects a lor file with: hist_num, c_x, c_y, c_z, d_x, d_y, d_z, sigma_l, sigma_t
lor* read_lor(FILE* input) {
	if (input == NULL) {
		return NULL;
	}

	double center_x;
	double center_y;
	double center_z;
	double dir_x;
	double dir_y;
	double dir_z;
	double longitudinal;
	double transverse;

	int worked;

	fscanf(input, "%*i,%lf,", &center_x);
	fscanf(input, "%lf,", &center_y);
	fscanf(input, "%lf,", &center_z);
	fscanf(input, "%lf,", &dir_x);
	fscanf(input, "%lf,", &dir_y);
	fscanf(input, "%lf,", &dir_z);
	fscanf(input, "%lf", &longitudinal);
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

render* read_render_def(FILE* input) {
	if (input == NULL) {
		fprintf(stderr, "read_render_def: null pointer given to render definitions");
		return NULL;
	}
	double x, y, z;

	fscanf(input, "%lf,%lf,%lf,", &x, &y, &z); // render cube corner 1
	vec3* corner1 = three_vec(x, y, z);
	fscanf(input, "%lf,%lf,%lf", &x, &y, &z); // render cube corner 2
	vec3* corner2 = three_vec(x, y, z);
	vec3* low_corner = three_vec(0,0,0);
	vec3* high_corner = three_vec(0,0,0);
	if (corner1->x < corner2->x) {
		low_corner->x = corner1->x;
		high_corner->x = corner2->x;
	} else {
		low_corner->x = corner2->x;
		high_corner->x = corner1->x;
	}
	if (corner1->y < corner2->y) {
		low_corner->y = corner1->y;
		high_corner->y = corner2->y;
	} else {
		low_corner->y = corner2->y;
		high_corner->y = corner1->y;
	}	
	if (corner1->z < corner2->z) {
		low_corner->z = corner1->z;
		high_corner->z = corner2->z;
	} else {
		low_corner->z = corner2->z;
		high_corner->z = corner1->z;
	}
	free(corner1);
	free(corner2);
	int a, b, c;
	fscanf(input, "%i,%i,%i", &a, &b, &c); // gets the number of voxels in each dimension

	// define the conversion from length space to voxel space
	double x_convert = a / (high_corner->x - low_corner->x);
	double y_convert = b / (high_corner->y - low_corner->y);
	double z_convert = c / (high_corner->z - low_corner->z);


	render* new = (render*)malloc(sizeof(render));
	double* volume = (double*)malloc(sizeof(double) * a * b * c);
	// allocate the volume for all of our rendering

	new->volume = volume;
	new->dimensions[0] = a;
	new->dimensions[1] = b;
	new->dimensions[2] = c;
	new->conversion = three_vec(x_convert, y_convert, z_convert);
	new->least_corner = low_corner;
	new->max_corner = high_corner;


	// add stuff to handle various choices of conversion here when needed.

	new->combiner == NULL;
	return new;
}


vec3* space_to_vox(vec3* a, render* universe) {
	if (a = NULL || universe == NULL) {
		return NULL;
	}
	// first handle the offset
	vec3* offsetless = vec_sub(a, universe->least_corner);
	// now multiply by conversion
	offsetless->x = offsetless->x * universe->conversion->x;
	offsetless->y = offsetless->y * universe->conversion->y;
	offsetless->z = offsetless->z * universe->conversion->z;

	return offsetless;
}



















int main(int argc, char const *argv[])
{
	// expected inputs: lor file name, output file name, rendering def file name
	if (argc != 4) {
		if (strcmp(argv[1], "-h") || strcmp(argv[1], "-H")) {
			printf("Kepler's lor renderer: \nThis is designed to be used with");
			printf(" the reverse kinematics code running on a TOPAS simulation.");
			printf("\nThe code expects three files: a lor file, an output file ");
			printf(" name, and a definitions file.\n");
			printf("The lor file should have one line of response for each history.");
			printf(" Each lor should have the struture:\n\thistory_number,");
			printf(" center_x, center_y, center_z, direction_x, dir_y, dir_z,");
			printf(" sigma_longitudianal, sigma_transverse\n\n");
			printf("The output file is just the name you wish for it to have.\n");
			printf("The defintion file is a set of values:\n\tcorner_1_x, ");
			printf("corner_1_y, corner_1_z, corner_2_x, corner_2_y, corner_2_z\n\t");
			printf("voxels_in_x, voxels_in_y, voxels_in_z");
			// add information on selecting combination method and output type.
		}
		printf("Expected 3 inputs, a lor file name, an output file name and");
		printf(" a rendering defintion file name.");
	}
	return 0;
}
