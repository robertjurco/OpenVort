/**
    Code for initialization and injecion of vorticies.
    @file vortex_injection.c
*/

#include <math.h>

#include "vortex_injection.h"
#include "vortex_constants.h"
#include "boundaries.h"
#include "tangle.h"

/**
    Creates a perpendicular vector to the given direction vector. This can be used to create a base of perpendicular vectors to the given direction vector.
    @param dir: Vector to witch we want to find the perpendiclar one.
    @return Vector perpendiulr to the given direction vector.
*/
struct vec3 perpendicular(const struct vec3 *dir)
{
    // Definitions
    struct vec3 a = vec3(1,0,0);
    struct vec3 b = vec3(0,1,0);
    struct vec3 res;

    // Perpendicular vector to *dir from the plane given by *dir and (1,0,0).
    vec3_mul(&res, dir, -vec3_dot(&a, dir));
    vec3_add(&res, &res, &a);
    // If *dir is not colinear with (1,0,0).
    if(vec3_len(&res) > 1e-8)
    {
        vec3_normalize(&res);
        return res;
    }

    // If *dir is colinear with (1,0,0) we use vector perpedicular to *dir from the plane given by *dir and (0,1,0).
    vec3_mul(&res, dir, -vec3_dot(&b, dir));
    vec3_add(&res, &res, &b);
    vec3_normalize(&res);
    return res;
}

/**
    Add one loop vortex into the tangle at the point center, direction dir and with radius r.
    @param tangle: Tangle structure holding all the information about the vorices.
    @param center: Center of the loop vortex given by vector from the origin.
    @param dir: Direction of the loop vortex given by vector.
    @param r: Radius of the circle in centimeters.
*/
void add_loop(struct tangle_state *tangle, struct vec3 *center, struct vec3 *dir, double r)
{
    // Calculate number of points.
    int Npoints = (int) (4 * M_PI * r / (global_dl_max + global_dl_min));
    // Is there enough of free points in the tangle?
    while (number_of_empty_nodes(tangle) < Npoints) expand_tangle(tangle, 2 * tangle->N);

    // Declare and define variables.
    struct vec3 u, v, p, ptmp;
    int curr_point, last_point, first_point;

    // Create a base for circle vectors.
    struct vec3 zdir = perpendicular(dir);
    vec3_cross(&u, &zdir, dir);
    vec3_cross(&v, dir, &u);
    vec3_normalize(&u);
    vec3_normalize(&v);

    first_point = -1;
	curr_point = -1; // To suppres warnings.

    // Loop through the points of the loop.
    for(int k = 0; k < Npoints; ++k)
    {
        // Find position of the k-th point, store it into vector p.
        double phi = 2*M_PI/Npoints * k;
        double c = cos(phi);
        double s = sin(phi);
        vec3_mul(&p, &u, r*c);
        vec3_mul(&ptmp, &v, r*s);

        vec3_add(&p, &p, &ptmp);
        vec3_add(&p, &p, center);
        curr_point = get_next_empty_node(tangle);

        // Check if it is the first point.
        if (first_point < 0) first_point = curr_point;

        // Create the k-th point.
        tangle->status[curr_point].status = FREE;
        tangle->vnodes[curr_point] = p;

        // Connect k-th point with the previous point, exept if it is the first point.
        if(curr_point != first_point)
	    {
            tangle->connections[curr_point].reverse = last_point;
            tangle->connections[last_point].forward = curr_point;
        }
        // Rememeber the last point.
        last_point = curr_point;
    }
    // Connect last and first points.
    tangle->connections[curr_point].forward  = first_point;
    tangle->connections[first_point].reverse = curr_point;
}

/**
    Add one loop vortex into the tangle at the point center, direction dir, with radius r and with helical osilation of given wavelength and amplitude.
    @param tangle: Tangle structure holding all the information about the vorices.
    @param center: Center of the loop vortex given by vector from the origin.
    @param dir: Direction of the loop vortex given by vector.
    @param r: Radius of the circle in centimeters.
    @param amplitude: Amplitude of the helical oscilation.
    @param wavelength: Wavelength of the helical oscilation.
*/
void add_loop_with_oscilation(struct tangle_state* tangle, struct vec3* center, struct vec3* dir, double r, double amplitude, double wavelength)
{
    // Calculate number of points.
    int Npoints = (int)(4 * M_PI * r / (global_dl_max + global_dl_min));
    // Is there enough of free points in the tangle?
    while (number_of_empty_nodes(tangle) < Npoints) expand_tangle(tangle, 2 * tangle->N);

    // Declare and define variables.
    struct vec3 u, v, p, ptmp, paral, perp;
    int curr_point, last_point, first_point;
    struct vec3 zdir = perpendicular(dir);

    vec3_cross(&u, &zdir, dir);
    vec3_cross(&v, dir, &u);

    vec3_normalize(&u);
    vec3_normalize(&v);

    first_point = -1;
	curr_point = -1; // To suppres warnings.

    // Loop through the points of the loop.
    for (int k = 0; k < Npoints; ++k)
    {
        // Find position of the k-th point, store it into vector p.
        double phi = 2 * M_PI / Npoints * k;
        double c = cos(phi);
        double s = sin(phi);
        vec3_mul(&p, &u, r * c);
        vec3_mul(&ptmp, &v, r * s);

        vec3_add(&p, &p, &ptmp);

        // Add helical oscilation.
        paral = p;
        perp = *dir;
        vec3_normalize(&paral);
        vec3_normalize(&perp);

        vec3_mul(&paral, &paral, amplitude * cos(2 * M_PI * phi * r / wavelength));
        vec3_mul(&perp, &perp, amplitude * sin(2 * M_PI * phi * r / wavelength));

        vec3_add(&p, &p, &paral);
        vec3_add(&p, &p, &perp);

        // Move loop into its posiion.
        vec3_add(&p, &p, center);

        // Find next free node inside the tangle.
        curr_point = get_next_empty_node(tangle);

        // Check if it is the first point.
        if (first_point < 0) first_point = curr_point;

        // Create the k-th point.
        tangle->status[curr_point].status = FREE;
        tangle->vnodes[curr_point] = p;

        // Connect k-th point with the previous point, exept if it is the first point.
        if (curr_point != first_point)
        {
            tangle->connections[curr_point].reverse = last_point;
            tangle->connections[last_point].forward = curr_point;
        }
        // Rememeber the last point.
        last_point = curr_point;
    }
    // Connect last and first points.
    tangle->connections[curr_point].forward = first_point;
    tangle->connections[first_point].reverse = curr_point;
}