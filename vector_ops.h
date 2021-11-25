#ifndef vector_ops_h
#define vector_ops_h

#include <stdio.h>

typedef struct vec3_ {
	double x;
	double y;
	double z;
} vec3;

// defines a new 3-vector
vec3* three_vec(double x, double y, double z);

// returns the magnitude of a 3-vector
double vec_mag(vec3* vector);

// returns the dot product of two 3-vectors
double vec_dot(vec3* a, vec3* b);

// returns the addition of two 3-vectors as a new 3-vector
vec3 *vec_add(vec3 *a, vec3* b);

// returns the subtraction of the second vector from the first
vec3 *vec_sub(vec3* a, vec3* b);

// returns the cross product of the two vectors a X b
vec3 *vec_cross(vec3* a, vec3* b);

// calculates the angle between two 3 vectors
double vec_angle(vec3* a, vec3* b);

// determines the distance between two vectors
double vec_dist(vec3* a, vec3* b);

// makes a new copy of a vector
vec3* vec_copy(vec3* a);

// prints a vector as the three values
void vec_print(vec3* a, FILE* output);

#endif