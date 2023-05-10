#include "vector_ops.h"
#include "ray_trace.h"
#include <stdlib.h>
#include <math.h>


/* 
 * Simple ray tracing code to allow for attenuation correction in the LOR
 * rendering operations.
 * Allows for the building of a combined shape using subshapes. The current
 * allowed shapes are the sphere, cylinder, and rectangular prism. The shapes
 * are defined by their centerpoint, their dimentions, and in the case of the
 * cylinder: what axis the circle cross-section is normal to. Also includes
 * information about the attenuation coefficent, however it will produce
 * undefined behavior if the shape has different attenuation coefficents at
 * different locations.
 * 
 * Allows for overlapping geometry, as the method of ray tracing is taking a
 * ray, checking if it is inside the geometry, then calculating the next
 * intersection point in the direction of travel. When it exits the geometry
 * object it is currently in, it checks if it is inside of any other objects.
 * If it is, then the distance to the next exit is found and added to the
 * running total.
 * 
 * If a ray does not intersect a shape in the given direction then a NULL
 * pointer is returned by the intersection function of the shape.
 * 
 * IN OPERATION, ALL RAY DIRECTION VECTORS MUST BE UNIT VECTORS
 * 
 */

// builds a ray structure
ray* ray_build(vec3d pos, vec3d dir) {
    ray* new = (ray*)malloc(sizeof(ray));
    if (new == NULL) {
        fprintf(stderr, "make_ray: malloc failed, likely fatal\n");
        return NULL;
    }
    new->pos = pos;
    new->dir = dir;
    return new;
}

// copy a ray structure
ray* ray_copy(ray* src) {
    if (src == NULL) {
        return NULL;
    }
    return ray_build(src->pos, src->dir);
}

// frees a ray structure
void* ray_free(ray* src) {
    free(src);
    return NULL;
}

traversal* traversal_build(vec3d intersect, double t) {
    traversal* new = (traversal*)malloc(sizeof(traversal));
    if (new == NULL) {
        fprintf(stderr, "traversal_build: malloc failed, likely fatal\n");
        return NULL;
    }
    new->intersection = intersect;
    new->t = t;
    return new;
}

traversal* traversal_copy(traversal* src) {
    return traversal_build(src->intersection, src->t);
}

void traversal_free(traversal* src) {
    if (src != NULL) {        
        free(src);
    }
}

shape* shape_build(int type, float* pos, float* dim, int axis, float attenuation) {
    if ((pos == NULL) || (dim == NULL)) {
        return NULL;
    }
    shape* new = (shape*)malloc(sizeof(shape));
    new->type = type;
    new->axis = axis;
    new->atten = attenuation;
    for (int i = 0; i < 3; i++) {
        new->pos[i] = pos[i];
    }
    if (type == REC_PRISM) {
        for (int i = 0; i < 3; i++) {
            new->dim[i] = dim[i];
        }
    } else if (type == CYLINDER) {
        for (int i = 0; i < 2; i++) {
            new->dim[i] = dim[i];
        }
    } else if (type == SPHERE) {
        new->dim[0] = dim[0];
    } else {
        fprintf(stderr, "shape_build: passed non-valid shape\n");
        free(new);
        return NULL;
    }
    return new;
}

void* shape_free(void* a) {
    if (a == NULL) {
        return NULL;
    }
    free(a);
    return(NULL);
}

geometry* geometry_build(shape** geos, uint size) {
    geometry* new = malloc(sizeof(geometry));
    new->geo = geos;
    new->size = size;
    return new;
}

void* geometry_free(geometry* a) {
    if (a != NULL) {
        for (int i = 0; i < a->size; i++) {
            free(a->geo[i]);
        }
        free(a->geo);
        free(a);
    }
    return NULL;
}

/*
 * Takes the vector and adjusts the values to transfer them to the given axis.
 * The change is moving the z axis to the given location. 
 * (1 = x, 2 = y, 3=z(do nothing)). Negative values are also accepted. A
 * negative value is interpreted as a rotation in the opposite direction. If it 
 * was able to do the transformation a 1 is returned. If it failed in any way it
 * returns a 0.
 */
int coord_swap(vec3d* coord, int axis) {
    if ((coord == NULL) || (axis > 3) || (axis < -3) || (axis == 0)) {
        return 0;
    }
    if ((axis == 1) || (axis == -1)) {
        // move the z axis onto the x axis 
        double x_old = coord->x;
        double z_old = coord->z;
        coord->x = ((double)axis) * z_old;
        coord->z = -((double)axis) * x_old;
        return 1;
    } else if ((axis == 2) || (axis == -2)) {
        // move the z axis onto the y axis
        double y_old = coord->y;
        double z_old = coord->z;
        coord->y = ((double)(axis/abs(axis))) * z_old;
        coord->z = -((double)(axis/abs(axis))) * y_old;
        return 1;
    } else {
        // axis is 3, for 3 we keep the z where it is, same for -3
    }
    return 1;
}

/*
 * finds if the given ray (propagate) intersects a given plane. The plane is
 * defined by the normal axis (axis = 1,2,3 = x,y,z normal) and the size in
 * floats. The normal axis is to be treated as the z axis, with the size 
 * defining the size along the x (size[0]) and y (size[1]) directions. This plane
 * is centered at the location given by center.
 */
traversal* plane_intersect_rec(float size[2], int axis, vec3d center_src, ray* propagate_src) {
    if ((size == NULL) || (propagate_src == NULL)) {
        return NULL;
    }
    // rotate everything so that it is all lined up
    vec3d center = center_src;
    ray* propagate = ray_copy(propagate_src);
    coord_swap(&center, axis);
    coord_swap(&(propagate->pos), axis);
    coord_swap(&(propagate->dir), axis);
    // find when the intersection should happen (when ray crosses z of plane)
    double t = (center.z - propagate->pos.z) / propagate->dir.z;
    double x_cross = propagate->pos.x + (propagate->dir.x * t);
    double x_hl = 0.5 * size[0];
    if ((x_cross > (x_hl + center.x)) || (x_cross < (-x_hl + center.x))) {
        // we missed the x edge of the plane, no intersection
        ray_free(propagate);
        return NULL;
    }
    double y_cross = propagate->pos.y + (propagate->dir.y * t);
    double y_hl = 0.5 * size[1];
    if ((y_cross > (y_hl + center.y)) || (y_cross < (-y_hl + center.y))) {
        ray_free(propagate);
        return NULL;
    }
    vec3d intersection = three_vec(x_cross, y_cross, center.z);
    // now that we have the intersection we need to rotate everything back to
    // its original position (negative axis value)
    coord_swap(&intersection, -axis);
    ray_free(propagate);
    traversal* output = traversal_build(intersection, t);
    return output;
}

/*
 * finds if the given ray (propagate) intersects a given plane. The plane is
 * defined by the normal axis (axis = 1,2,3 = x,y,z normal) and the size in
 * floats. The normal axis is to be treated as the z axis, with the size 
 * defining the radius as size. This plane
 * is centered at the location given by center.
 */

traversal* plane_intersect_circle(float size, int axis, vec3d center_src, ray* propagate_src) {
    if (propagate_src == NULL) {
        return NULL;
    }
    // rotate everything so that it is all lined up (as if the normal was zhat)
    vec3d center = center_src;
    ray* propagate = ray_copy(propagate_src);
    coord_swap(&center, axis);
    coord_swap(&(propagate->pos), axis);
    coord_swap(&(propagate->dir), axis);
    // find when the intersection should happen (when ray crosses z of plane)
    double t = (center.z - propagate->pos.z) / propagate->dir.z;
    double x_cross = propagate->pos.x + (propagate->dir.x * t);
    double y_cross = propagate->pos.y + (propagate->dir.y * t);
    double rad2 = size * size;
    double x_disp = x_cross - center.x;
    double y_disp = y_cross - center.y;
    double disp2 = (x_disp * x_disp) + (y_disp * y_disp);
    if (rad2 < disp2) {
        // intersection not within the given radius
        ray_free(propagate);
        return NULL;
    }
    vec3d intersection = three_vec(x_cross, y_cross, center.z);
    // now that we have the intersection we need to rotate everything back to
    // its original position (negative axis value)
    coord_swap(&intersection, -axis);
    ray_free(propagate);
    traversal* output = traversal_build(intersection, t);
    return output;
}

/*
 * Takes a ray and prism and finds how far the ray traveled (in the direction of
 * travel) through the prism. Returns the exit location and distance travled.
 * If the traversal started outside the box (and therefore there are two crossings
 * in the direction of travel) the total distance inside the box is returned,
 * with the farther exit location and full_crossing is set to true.
 */
traversal* exit_rectangular_prism(ray* path, shape* prism, int* full_crossing) {
    if ((path == NULL) || (prism == NULL) || (prism->type != REC_PRISM)) {
        return NULL;
    }
    if (full_crossing != NULL) {
        full_crossing[0] = 0; // the interaction did not start outside of the box
    }
    double long_path = -1.0;
    double short_path = -1.0;
    traversal* best = NULL;
    // check the distance for all six faces.
    for (int i = 0; i < 3; i++) {
        // check the face in the + axis direction and the -axis direction
        // each face has the dimensions of the two axis not using i
        float size[2];
        if (i == 0) {
            // rotation from +x to +z
            size[0] = prism->dim[2]; // x = dim[z]
            size[1] = prism->dim[1]; // y = dim[y]
        } else if (i == 1) {
            // rotation from +y to +z
            size[0] = prism->dim[0]; // x = dim[x]
            size[1] = prism->dim[2]; // y = dim[z]
        } else /*i == 2*/ {
            // rotation from +z to +z
            size[0] = prism->dim[0]; // x = dim[x]
            size[1] = prism->dim[1]; // y = dim[y]
        }
        // get the offset of the plane we will check from the center of the shape
        double offset = prism->dim[i] * 0.5;
        // centerpoints of the planes
        vec3d direction_p = three_vec(prism->pos[0] + (offset * (i == 0)),
                                        prism->pos[1] + (offset * (i == 1)),
                                        prism->pos[2] + (offset * (i == 2)));
        vec3d direction_m = three_vec(prism->pos[0] - (offset * (i == 0)),
                                        prism->pos[1] - (offset * (i == 1)),
                                        prism->pos[2] - (offset * (i == 2)));
        traversal* exit1 = plane_intersect_rec(size, i + 1, direction_p, path);
        traversal* exit2 = plane_intersect_rec(size, i + 1, direction_m, path);

        // now check if either has an intersection at all,
        if ((exit1 != NULL) && (exit1->t > 0)) {
            if (long_path < 0) {
                long_path = exit1->t;
                best = exit1;
            } else {
                if (exit1->t > long_path) {
                    traversal_free(best);
                    short_path = long_path;
                    long_path = exit1->t;
                    best = exit1;
                } else {
                    short_path = exit1->t;
                    traversal_free(exit1);
                }
            }
        } else {
            traversal_free(exit1);
        }
        if ((exit2 != NULL) && (exit2->t > 0)) {
            if (long_path < 0) {
                long_path = exit2->t;
                best = exit2;
            } else {
                if (exit2->t > long_path) {
                    traversal_free(best);
                    short_path = long_path;
                    long_path = exit2->t;
                    best = exit2;
                } else {
                    short_path = exit2->t;
                    traversal_free(exit2);
                }
            }
        } else {
            traversal_free(exit2);
        }
    }
    // we have now checked all faces. If there is only one traversal with a
    // positive value (the original ray was inside of the box) then the distance
    // to its exit is given. If two returned a value (the ray started outside, 
    // and passed through), then the difference between them is returned
    if (best == NULL) {
        // no intersection found!
        return NULL;
    }
    if (short_path > 0) {
        // passed completely through
        if (full_crossing != NULL) {
            full_crossing[0] = 1;
        }
        best->t = long_path - short_path; // set traversal to dist inside
        return best;
    }
    return best;
}

// returns the two crossing t of a sphere (or circle if you project to 2d)
// returns null if it is not 
double* sphere_crossing(ray* path, vec3d center, double r) {
    if (path == NULL) {
        return NULL;
    }
    vec3d st_ls_sph = vec_sub(path->pos, center);
    double u_dot_st_ls_sph = vec_dot(path->dir, st_ls_sph);
    double determinator = (u_dot_st_ls_sph * u_dot_st_ls_sph)
                        - (vec_dot(st_ls_sph, st_ls_sph) - (r * r));
    if (determinator <= 0) {
        // no distance spent inside of the sphere
        return NULL;
    }
    double determined = sqrt(determinator);
    double t_high = -u_dot_st_ls_sph + determined;
    double t_low  = -u_dot_st_ls_sph - determined;
    double* crossings = (double*)malloc(2 * sizeof(double));
    crossings[0] = t_low;
    crossings[1] = t_high;
    return crossings;
}

/*
 * finds the distance to exiting the sphere for the given ray. If the ray
 * starts outside the sphere and travels through, the distance is given as the
 * travel inside of the sphere and full_crossing is set to true. Works using
 * a line crossing sphere formula on wikipedia
 * 
 * Assumes that the ray is a unit vector
 */
traversal* exit_sphere(ray* path, shape* sphere, int* full_crossing){
    if ((path == NULL) || (sphere == NULL) || (sphere->type != SPHERE)) {
        return NULL;
    }
    if (full_crossing != NULL) {
        full_crossing[0] = 0;
    }
    vec3d sphere_center = three_vec(sphere->pos[0], sphere->pos[1], sphere->pos[2]);
    double* crossings = sphere_crossing(path, sphere_center, sphere->dim[0]);
    if (crossings == NULL) {
        return NULL;
    }
    double t_high = crossings[1];
    double t_low  = crossings[0];
    free(crossings);
    if (t_low > 0) {
        vec3d dist = vec_scaler(path->dir, t_high);
        traversal* exit = traversal_build(vec_add(dist, path->pos), t_high - t_low);
        if (full_crossing != NULL) {
            full_crossing[0] = 1;
        }
        return exit;
    } else if (t_high > 0) {
        vec3d dist = vec_scaler(path->dir, t_high);
        traversal* exit = traversal_build(vec_add(dist, path->pos), t_high);
        return exit;
    }
    return NULL;
}

/* 
 * Finds the distance to exiting a cylinder from the current ray. If the ray
 * starts outside and goes through it returns the distance traveled inside and
 * sets the full crossing flag to true. Works by checking for intersections with
 * the two flat ends, then looking for ones along the cylinder. The cylinder
 * check reuses the code for the sphere, just projected onto the xy plane.
 */
traversal* exit_cyl(ray* path, shape* cyl, int* full_crossing) {
    if ((path == NULL) || (cyl == NULL) || (cyl->type != CYLINDER)) {
        return NULL;
    }
    if (full_crossing != NULL) {
        full_crossing[0] = 0;
    }
    vec3d center = three_vec(cyl->pos[0], cyl->pos[1], cyl->pos[2]);
    double half_height = cyl->dim[1] * 0.5;
    vec3d offset = three_vec(half_height * (cyl->axis == 1),
                                half_height * (cyl->axis == 2),
                                half_height * (cyl->axis == 3));
    vec3d plane1 = vec_add(center, offset);
    vec3d plane2 = vec_sub(center, offset);
    traversal* exit1 = plane_intersect_circle(cyl->dim[0], cyl->axis, plane1, path);
    traversal* exit2 = plane_intersect_circle(cyl->dim[0], cyl->axis, plane2, path);
    traversal* plane_intersect = NULL;
    if ((exit1 != NULL) && (exit1->t > 0)) {
        plane_intersect = exit1;
        if ((exit2 != NULL) && (exit2->t > 0)) {
            if (full_crossing != NULL) {
                full_crossing[0] = 1;
            }
            if (exit1->t > exit2->t) {
                double dist = exit1->t - exit2->t;
                exit1->t = dist;
                traversal_free(exit2);
                return exit1;
            } else {
                double dist = exit2->t - exit1->t;
                exit2->t = dist;
                traversal_free(exit1);
                return exit2;
            }
        }
        traversal_free(exit2);
    } else if ((exit2 != NULL) && (exit2->t > 0)) {
        plane_intersect = exit2;
        traversal_free(exit1);
    } else {
        traversal_free(exit1);
        traversal_free(exit2);
    }
    // now for projected cylinder/sphere work. First rotate everything so that
    // the cylinder points in the z axis
    vec3d ray_dir = path->dir;
    vec3d ray_pos = path->pos;
    coord_swap(&center, cyl->axis);
    coord_swap(&ray_dir, cyl->axis);
    coord_swap(&ray_pos, cyl->axis);
    // project onto xy plane (but keep the z avaliable for later)

    double z_dir = ray_dir.z;
    ray_dir.z = 0.0;
    double z_pos = ray_pos.z;
    ray_pos.z = 0.0;
    double z_center = center.z;
    center.z = 0.0;
    double flat_mag = vec_mag(ray_dir);
    double inverse_mag = 1.0 / flat_mag;
    ray_dir = vec_scaler(ray_dir, inverse_mag);
    ray* new_path = ray_build(ray_pos, ray_dir);
    double* flat_cyl = sphere_crossing(new_path, center, cyl->dim[0]);
    if (flat_cyl == NULL) {
        // no sphere crossing occured
        ray_free(new_path);
        return plane_intersect;
    } else {
        // calculate the z intersection points
        // first correct the lengths given by the normalization used
        flat_cyl[0] *= inverse_mag;
        flat_cyl[1] *= inverse_mag;
        // re-add z component of position and direction for endpoint calcs
        new_path->pos.z = z_pos;
        new_path->dir.z = z_dir;
        center.z = z_center;
        // adjust x and y directions back to unit values.
        new_path->dir.x *= flat_mag;
        new_path->dir.y *= flat_mag;
        vec3d ends[2];
        char ends_exist[2] = {0,0};
        int num_end = 0;
        for (int i = 0; i < 2; i++) {
            if (flat_cyl[i] > 0) {
                // now to find the intersection point. We are still in the case
                // where we can say the cylinder points in the z direction
                vec3d travel = vec_scaler(new_path->dir, flat_cyl[i]);
                ends[i] = vec_add(travel, new_path->pos);

                if (fabs(ends[i].z - center.z) > half_height) {
                    ends[i] = three_vec(NAN,NAN,NAN);
                    ends_exist[i] = 0;
                    // removes the end if it is not within the cylinder height
                } else {
                    num_end++;
                    ends_exist[i] = 1;
                }
            } else {
                ends[i] = three_vec(NAN,NAN,NAN);
                ends_exist[i] = 0;
            }
        }
        ray_free(new_path);
        coord_swap(&(ends[0]), -(cyl->axis));
        coord_swap(&(ends[1]), -(cyl->axis));
        // we now have the available ends: plane_intersect, ends[0], ends[1]
        // no more than 2 of 3 can exist
        if (num_end == 0) {
            // no interaction with the walls occured in positive t.
            free(flat_cyl);
            return plane_intersect;
        }
        if (num_end == 1) {
            // one interaction occured, could be 0 or 1
            traversal* out;
            if (ends_exist[0] != 0) {
                out = traversal_build(ends[0], flat_cyl[0]);
            } else {
                out = traversal_build(ends[1], flat_cyl[1]);
            }
            free(flat_cyl);
            if ((plane_intersect != NULL) && (plane_intersect->t > 0)) {
                if (full_crossing != NULL) {
                    full_crossing[0] = 1;
                }
                if (plane_intersect->t > out->t) {
                    plane_intersect->t = plane_intersect->t - out->t;
                    traversal_free(out);
                    return plane_intersect;
                } else {
                    out->t = out->t - plane_intersect->t;
                    traversal_free(plane_intersect);
                    return out;
                }
            }
            traversal_free(plane_intersect);
            return out;
        } else {
            // we passed in one round side, out the other
            if (full_crossing != NULL) {
                full_crossing[0] = 1;
            }
            traversal_free(plane_intersect); // should never run, but just in case
            double dist = abs(flat_cyl[0] - flat_cyl[1]);
            if (flat_cyl[0] > flat_cyl[1]) {
                free(flat_cyl);
                traversal* out = traversal_build(ends[0], dist);
                return out;
            } else {
                free(flat_cyl);
                traversal* out = traversal_build(ends[1], dist);
                return out;
            }
        }
    }
}


/*
 * takes a ray and finds how far it spends inside of the given geometry. To do
 * this it steps through each object in the geometry. It checks to see if any
 * of the geometry has the path starting within the object. If it does, it marks
 * that geometry as crossed and starts again with the same requirements. If
 * the ray does not start in any of the geometry then the distance to all exits
 * is compared. The shortest distance from the ray start to the ray end less the
 * ray travel inside of the object is used as the next distance, and the
 * procedure continues from the exitpoint.
 */
double propagate(ray* path_src, geometry* all) {
    if ((path_src == NULL) || (all == NULL)) {
        // fprintf(stderr, "propagate: given bad geometry or ray\n");
        return -1; // errorcode
    }
    ray* path = ray_copy(path_src);
    char* mask = (char*)calloc(sizeof(char), all->size);
    // will act as a bitmask to show what geometry componets have already been
    // searched
    traversal** crossings = (traversal**)calloc(all->size, sizeof(traversal*));
    // will store the distances in the case where no start inside of geometry
    // was found. Starts with all NULL to simplify coding
    double distance = 0;
    // search all of the geometry
    double closest = INFINITY;
    int best_find = -1;
    for (int i = 0; i < all->size; i++) {
        // have we already looked at this geometry?
        if (!mask[i]) {
            // run the check for the next crossing
            shape* next_shape = all->geo[i];
            int full;
            if (next_shape->type == REC_PRISM) {
                // we have a box of some sort
                crossings[i] = exit_rectangular_prism(path, next_shape, &full);
            } else if (next_shape->type == SPHERE) {
                crossings[i] = exit_sphere(path, next_shape, &full);
            } else if (next_shape->type == CYLINDER) {
                crossings[i] = exit_cyl(path, next_shape, &full);
            } else {
                fprintf(stderr, "propagate: unknown shape type\n");
                crossings[i] = NULL;
            }
            if ((crossings[i] != NULL) && (full == 0)) {
                // we had an inside start, we can just use the given information
                distance += crossings[i]->t * next_shape->atten;
                mask[i] = 1; // mask this geometry off
                // move the ray
                path->pos = crossings[i]->intersection;
                // clean up the crossings array
                for (int j = 0; j <= i; j++) {
                    traversal_free(crossings[j]);
                    crossings[j] = NULL;
                }
                // restart the search at the beginning
                i = 0;
                closest = INFINITY;
                best_find = -1;
            } else if (crossings[i] != NULL) {
                // we had a full crossing, is the entrance closer than any
                // previous one
                double entry_dist = vec_dist(path->pos, crossings[i]->intersection) - crossings[i]->t;
                if (entry_dist < closest) {
                    closest = entry_dist;
                    best_find = i;
                }
            } else {
                // no crossing occured in the direction of travel of the array,
                // so this geometry component will not factor into our future
                // work.
                mask[i] = 1;
            }
        }
        if ((i + 1 == all->size) && (best_find >= 0)) {
            // we have reached the end of the loop, but there is a crossing
            // the closest geometry entry happened with best_find.
            mask[best_find] = 1;
            distance += crossings[best_find]->t * (all->geo[best_find])->atten;
            path->pos = crossings[best_find]->intersection;
            // clean up the crossings array
            for (int j = 0; j <= i; j++) {
                traversal_free(crossings[j]);
                crossings[j] = NULL;
            }
            // restart the search at the beginning
            i = 0;
            closest = INFINITY;
            best_find = -1;
        }
    }
    for (int j = 0; j < all->size; j++) {
        traversal_free(crossings[j]);
        crossings[j] = NULL;
    }
    free(crossings);
    free(mask);
    ray_free(path);
    return distance;
}

// now for a group of test cases and similar functions

// checks rectangular prisms with a full crossing
int test_prism_1() {
    // build the prism for testing
    float pos[3] = {1.0, -1.0, 1.0};
    float dim[3] = {2.0, 4.0, 2.0};
    shape* box = shape_build(REC_PRISM, pos, dim, 0, 1.0);
    vec3d start = three_vec(1.0, 3.0, -2.0);
    vec3d point = three_vec(0.0, -2.0, 1.0);
    vec3d unit_point = vec_norm(point);

    ray* path = ray_build(start, unit_point);
    geometry* world = geometry_build(&box, 1);
    double dist = propagate(path, world);
    int full;
    traversal* test = exit_rectangular_prism(path, box, &full);

    free(box);
    free(world);
    ray_free(path);

    int fail = 0;
    if (!full) {
        fprintf(stderr, "test_prism_1: full crossing not reported\n");
        fail++;
    }
    double correct = sqrt(5);
    double dir_test = fabs(test->t - correct);
    if (dir_test > 0.001) {
        fprintf(stderr, "test_prism_1: exit_rectangular_prism gives %lf, should be %lf\n",test->t, correct);
        fail++;
    }
    if (fabs(dist - correct) > 0.001) {
        fprintf(stderr, "test_prism_1: propagate gives wrong distance %lf, should be %lf\n", dist, correct);
        fail++;
    }
    if (fail) {
        fprintf(stderr, "test_prism_1: exiting with %i errors\n", fail);
    } else {
        printf("test_prism_1: passed tests\n");
    }
    traversal_free(test);
    return !fail;
}


// checks rectangualr prisms with a partial crossing
int test_prism_2() {
    // build the prism for testing
    float pos[3] = {1.0, -1.0, 1.0};
    float dim[3] = {22.0, 22.0, .1};
    shape* box = shape_build(REC_PRISM, pos, dim, 0, 1.0);
    vec3d start = three_vec(1.0, -7.0, 1.0);
    vec3d point = three_vec(0.0, 2.0, 0.0);
    vec3d unit_point = vec_norm(point);

    ray* path = ray_build(start, unit_point);
    geometry* world = geometry_build(&box, 1);
    double dist = propagate(path, world);
    int full;
    traversal* test = exit_rectangular_prism(path, box, &full);

    free(box);
    free(world);
    ray_free(path);

    int fail = 0;
    if (full) {
        fprintf(stderr, "test_prism_2: full crossing reported\n");
        fail++;
    }
    double correct = 17;
    double dir_test = fabs(test->t - correct);
    if (dir_test > 0.001) {
        fprintf(stderr, "test_prism_2: exit_rectangular_prism gives %lf, should be %lf\n",test->t, correct);
        fail++;
    }
    if (fabs(dist - correct) > 0.001) {
        fprintf(stderr, "test_prism_2: propagate gives wrong distance %lf, should be %lf\n", dist, correct);
        fail++;
    }
    if (fail) {
        fprintf(stderr, "test_prism_2: exiting with %i errors\n", fail);
    } else {
        printf("test_prism_2: passed tests\n");
    }
    traversal_free(test);
    return !fail;
}

// checks rectangular prisms with a full crossing
int test_prism_3() {
    // build the prism for testing
    float pos[3] = {0.0, 0.0, 0.0};
    float dim[3] = {22.0, 22.0, .1};
    shape* box = shape_build(REC_PRISM, pos, dim, 0, 1.0);
    vec3d start = three_vec(0.0, 1.0, 0.0);
    vec3d point = three_vec(1.0, 1.0, 0.0);
    vec3d unit_point = vec_norm(point);

    ray* path = ray_build(start, unit_point);
    geometry* world = geometry_build(&box, 1);
    double dist = propagate(path, world);
    int full;
    traversal* test = exit_rectangular_prism(path, box, &full);

    free(box);
    free(world);
    ray_free(path);

    int fail = 0;
    if (test == NULL) {
        fprintf(stderr, "test_prism_3: no intersection reported!\n");
        fail++;
        return 0;
    }

    if (full) {
        fprintf(stderr, "test_prism_3: full crossing reported\n");
        fail++;
    }
    double correct = sqrt((10 * 10) + (10 * 10));
    double dir_test = fabs(test->t - correct);
    if (dir_test > 0.001) {
        fprintf(stderr, "test_prism_3: exit_rectangular_prism gives %lf, should be %lf\n",test->t, correct);
        fail++;
    }
    if (fabs(dist - correct) > 0.001) {
        fprintf(stderr, "test_prism_3: propagate gives wrong distance %lf, should be %lf\n", dist, correct);
        fail++;
    }
    if (fail) {
        fprintf(stderr, "test_prism_3: exiting with %i errors\n", fail);
    } else {
        printf("test_prism_3: passed tests\n");
    }
    traversal_free(test);
    return !fail;
}



// checks spheres with a full crossing
int test_sphere_1() {
    // build the sphere for testing
    float pos[3] = {1.0, -1.0, 1.0};
    float dim[3] = {2.0, 4.0, 2.0};
    shape* sphere = shape_build(SPHERE, pos, dim, 0, 1.0);
    vec3d start = three_vec(0.0, 2.0, 1.0);
    vec3d point = three_vec(1.0, -1.0, 0.0);
    vec3d unit_point = vec_norm(point);

    ray* path = ray_build(start, unit_point);
    geometry* world = geometry_build(&sphere, 1);
    double dist = propagate(path, world);
    int full;
    traversal* test = exit_sphere(path, sphere, &full);

    free(sphere);
    free(world);
    ray_free(path);

    int fail = 0;
    if (!full) {
        fprintf(stderr, "test_sphere_1: full crossing not reported\n");
        fail++;
    }
    double correct = sqrt(8);
    double dir_test = fabs(test->t - correct);
    if (dir_test > 0.001) {
        fprintf(stderr, "test_sphere_1: exit_sphere gives %lf, should be %lf\n",test->t, correct);
        fail++;
    }
    if (fabs(dist - correct) > 0.001) {
        fprintf(stderr, "test_sphere_1: propagate gives wrong distance %lf, should be %lf\n", dist, correct);
        fail++;
    }
    if (fail) {
        fprintf(stderr, "test_sphere_1: exiting with %i errors\n", fail);
    } else {
        printf("test_sphere_1: passed tests\n");
    }
    traversal_free(test);
    return !fail;
}

// checks spheres with a partial crossing


// checks cylinder with a full crossing (in round, out plane)
int test_cyl_1() {
    // build the cylinder for testing
    float pos[3] = {0.0, 0.0, 0.0};
    float dim[3] = {1.0, 2.0, 2.0}; // radius 1, height 2
    shape* cyl = shape_build(CYLINDER, pos, dim, 1, 1.0);
    vec3d start = three_vec(-1.0, -2.0, 0.0);
    vec3d point = three_vec(1.0, 1.0, 0.0);
    vec3d unit_point = vec_norm(point);

    ray* path = ray_build(start, unit_point);
    geometry* world = geometry_build(&cyl, 1);
    double dist = propagate(path, world);
    int full;
    traversal* test = exit_cyl(path, cyl, &full);


    free(cyl);
    free(world);
    ray_free(path);

    if (test == NULL) {
        fprintf(stderr, "test_cyl_1: exit_sphere failed to give result\n");
        return 0;
    }

    int fail = 0;
    if (!full) {
        fprintf(stderr, "test_cyl_1: full crossing not reported\n");
        fail++;
    }
    double correct = sqrt(2);
    double dir_test = fabs(test->t - correct);
    if (dir_test > 0.001) {
        fprintf(stderr, "test_cyl_1: exit_sphere gives %lf, should be %lf\n",test->t, correct);
        fail++;
    }
    if (fabs(dist - correct) > 0.001) {
        fprintf(stderr, "test_cyl_1: propagate gives wrong distance %lf, should be %lf\n", dist, correct);
        fail++;
    }
    if (fail) {
        fprintf(stderr, "test_cyl_1: exiting with %i errors\n", fail);
    } else {
        printf("test_cyl_1: passed tests\n");
    }
    traversal_free(test);
    return !fail;
}

// checks cylinder with a partial crossing (out plane)
int test_cyl_2() {
    // build the cylinder for testing
    float pos[3] = {1.0, -1.0, 1.0};
    float dim[3] = {1.0, 2.0, 2.0}; // radius 1, height 2
    shape* cyl = shape_build(CYLINDER, pos, dim, 3, 1.0);
    vec3d start = three_vec(1.0, -1.0, 1.0);
    vec3d point = three_vec(0.0, 0.0, 1.0);
    vec3d unit_point = vec_norm(point);

    ray* path = ray_build(start, unit_point);
    geometry* world = geometry_build(&cyl, 1);
    double dist = propagate(path, world);
    int full;
    traversal* test = exit_cyl(path, cyl, &full);


    free(cyl);
    free(world);
    ray_free(path);

    if (test == NULL) {
        fprintf(stderr, "test_cyl_2: exit_cyl failed to give result\n");
        return 0;
    }

    int fail = 0;
    if (full) {
        fprintf(stderr, "test_cyl_2: full crossing reported incorrectly\n");
        fail++;
    }
    double correct = 1.0;
    double dir_test = fabs(test->t - correct);
    if (dir_test > 0.001) {
        fprintf(stderr, "test_cyl_2: exit_cyl gives %lf, should be %lf\n",test->t, correct);
        fail++;
    }
    if (fabs(dist - correct) > 0.001) {
        fprintf(stderr, "test_cyl_2: propagate gives wrong distance %lf, should be %lf\n", dist, correct);
        fail++;
    }
    if (fail) {
        fprintf(stderr, "test_cyl_2: exiting with %i errors\n", fail);
    } else {
        printf("test_cyl_2: passed tests\n");
    }
    traversal_free(test);
    return !fail;
}

// checks cylinder with a partial crossing (out round)



// the full suite of tests for ray tracing
int geometry_full_tests() {
    int num_of_tests = 6;
    int a = test_prism_1();
    a += test_prism_2();
    a += test_prism_3();
    a += test_sphere_1();
    a += test_cyl_1();
    a += test_cyl_2();
    printf("ray tracing passed %i out of %i tests\n", a, num_of_tests);
    return (a == num_of_tests);
}

void print_geometry(FILE* output, geometry* a) {
    uint i = a->size;
    fprintf(output, "%u shapes in geometry\n",i);
    for (uint j = 0; j < i; j++) {
        shape* cur = a->geo[j];
        fprintf(output, "\tShape type: ");
        if (cur->type == REC_PRISM) {
            fprintf(output, "rectangular prism\n");
            fprintf(output, "\t\tlx = %f, ly = %f, lz = %f\n", cur->dim[0],cur->dim[1],cur->dim[2]);
        } else if (cur->type == CYLINDER) {
            fprintf(output, "cylinder\n");
            fprintf(output, "\t\theight = %f, radius = %f\n", cur->dim[1], cur->dim[0]);
            char axis_array[4] = {'\0', 'x', 'y', 'z'};
            fprintf(output, "\t\tmain axis of shape: %c\n", axis_array[cur->axis]);
        } else if (cur->type == SPHERE) {
            fprintf(output, "sphere\n");
            fprintf(output, "\t\tradius = %lf\n", cur->dim[0]);
        } else {
            fprintf(output, "UNKNOWN SHAPE ERROR!\n");
        }
        fprintf(output, "\t\tx = %f,  y = %f,  z = %f\n", cur->pos[0],cur->pos[1],cur->pos[2]);
        fprintf(output, "\t\tattenuation: %f\n", cur->atten);
    }
}