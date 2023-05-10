#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "vector_ops.h"

// defines a new 3-vector
vec3d three_vec(double x, double y, double z) {
	vec3d new;
	new.x = x;
	new.y = y;
	new.z = z;
	return new;
}

// determines the distance between two vectors
double vec_dist(vec3d a, vec3d b) {
	// subtract b from a
	vec3d new = vec_sub(a, b);
	// get the distance of the subtracted vector
	double distance = vec_mag(new);
	return distance;
}

// returns the magnitude of a 3-vector
double vec_mag(vec3d vector) {
	double x2 = vector.x * vector.x;
	double y2 = vector.y * vector.y;
	double z2 = vector.z * vector.z;
	return sqrt(x2 + y2 + z2);
}

// returns the dot product of two 3-vectors
double vec_dot(vec3d a, vec3d b) {
	double first = a.x * b.x;
	double second = a.y * b.y;
	double third = a.z * b.z;
	return first + second + third;
}

// returns the addition of two 3-vectors as a new 3-vector
vec3d vec_add(vec3d a, vec3d b) {
	vec3d sum;
	sum.x = a.x + b.x;
	sum.y = a.y + b.y;
	sum.z = a.z + b.z;
	return sum;
}

// returns the subtraction of the second vector from the first
vec3d vec_sub(vec3d a, vec3d b) {
	vec3d sum;
	sum.x = a.x - b.x;
	sum.y = a.y - b.y;
	sum.z = a.z - b.z;
	return sum;
}

// returns the cross product of the two vectors a X b
vec3d vec_cross(vec3d a, vec3d b) {
	vec3d cross;
	cross.x = (a.y * b.z) - (a.z * b.y);
	cross.y = (a.x * b.z) - (a.z * b.x);
	cross.z = (a.x * b.y) - (a.y * b.x);
	return cross;
}

// calculates the angle between two 3 vectors in radians
double vec_angle(vec3d a, vec3d b) {
	// arccos((a * b) / (||a|| ||b||)) = theta
	double a_dot_b = vec_dot(a, b);
	double norms = vec_mag(a) * (vec_mag(b));
	return acos(a_dot_b / norms);
}

// // makes a new copy of a vector
// vec3d vec_copy(vec3d a) {
// 	return three_vec(a.x, a->y, a->z);
// }
// this is no longer needed now that we are
// just using the true structure rather than a pointer

// prints a vector as the three values
void vec_print(vec3d a, FILE* output) {
	if (output == NULL) {
		return;
	}
	fprintf(output, "%f, %f, %f", a.x, a.y, a.z);
}

// normalizes the given vector, returns as a new vector structure
vec3d vec_norm(vec3d a) {
	double mag = vec_mag(a);
	return three_vec(a.x / mag, a.y / mag, a.z / mag);
}

// multiplies a vector by a scalar
vec3d vec_scaler(vec3d a, double b) {
	return three_vec(b * a.x, b * a.y, b * a.z);
}

// projects a onto b
vec3d vec_projection(vec3d a, vec3d b) {
	return vec_scaler(b, (vec_dot(a,b)/vec_dot(b,b)));
}

// rejects b from a (i.e. projects a onto the plane perpendicular to b)
vec3d vec_rejection(vec3d a, vec3d b) {
	// find the projection of a onto b
	vec3d projection = vec_projection(a, b);
	return vec_sub(a, projection);
}

