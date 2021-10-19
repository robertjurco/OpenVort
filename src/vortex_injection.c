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

/**
    Add N random loop vortices into the tangle box with radius from interval [avg_radius*(1-variation), avg_radius*(1+variation)] and maximal distance from origin of tangle box given by dist.
    @param tangle: Tangle structure holding all the information about the vorices.
    @param N: Number of the generated loop vortices.
    @param avg_radius: Avarage radius of the generated loop vortices.
    @param variation: Gives the lower and upper bound of radius variation from its mean in fraction of avg_radius.
    @param dist: Maximal distance of centers of generated loop from origin of the tangle box.
*/
void add_random_loops_centered(struct tangle_state *tangle, int N, double avg_radius, double variation, double dist)
{
  for(int k = 0; k < N; ++k) {
      // Randomise direction of k-th loop.
      struct vec3 dir = vec3(0,0,0);
      while (dir.p[0] * dir.p[0] + dir.p[1] * dir.p[1] + dir.p[2] * dir.p[2] < 1e-8)
      {
          dir = vec3(drand48() - 0.5, drand48() - 0.5, drand48() - 0.5);
      }
      vec3_normalize(&dir);

      // Randomise center of k-th loop.
      struct vec3 c = vec3(0, 0, 0);
      while (c.p[0] * c.p[0] + c.p[1] * c.p[1] + c.p[2] * c.p[2] < 1e-8)
      {
          c = vec3(drand48() - 0.5, drand48() - 0.5, drand48() - 0.5);
      }
      vec3_normalize(&c);
      double d1 = dist * drand48();
      vec3_mul(&c,&c, d1);

      // Randomise radius of k-th loop.
      double d2 = 2 * (drand48() - 0.5);
      double r = avg_radius + avg_radius * variation * d2;

      // Add k-th loop.
      add_loop(tangle, &c, &dir, r);
    }
}

/**
    Add N random loop vortices into the tangle box with radius from interval [avg_radius*(1-variation), avg_radius*(1+variation)].
    @param tangle: Tangle structure holding all the information about the vorices.
    @param N: Number of the generated loop vortices.
    @param avg_radius: Avarage radius of the generated loop vortices.
    @param variation: Gives the lower and upper bound of radius variation from its mean in fraction of avg_radius.
 */
void add_random_loops(struct tangle_state *tangle, int N, double avg_radius, double variation)
{
  for(int k = 0; k<N; ++k) {
      // Randomise direction of k-th loop.
      struct vec3 dir = vec3(0, 0, 0);
      while (dir.p[0] * dir.p[0] + dir.p[1] * dir.p[1] + dir.p[2] * dir.p[2] < 1e-8)
      {
          dir = vec3(drand48() - 0.5, drand48() - 0.5, drand48() - 0.5);
      }
      vec3_normalize(&dir);
      
      // Randomise center of k-th loop.
      struct vec3 c;
      for (int i = 0; i < 3; ++i)
      {
          double blb = tangle->box.bottom_left_back.p[i] + avg_radius * (1 + variation);
          double trf = tangle->box.top_right_front.p[i] - avg_radius * (1 + variation);
          c.p[i] = blb + drand48() * (trf - blb);
      }

      // Randomise radius of k-th loop.
      double d2 = 2 * (drand48() - 0.5);
      double r = avg_radius + avg_radius * variation * d2;

      // Add k-th loop.
      add_loop(tangle, &c, &dir, r);
    }
}

/**
    Add straight line vortex in z direction at position x,y and direction up (+1) or down (-1).
    @param tangle: Tangle structure holding all the information about the vorices.
    @param x: Position of the points in the x direction.
    @param y: Position of the points in the y direction.
    @param direction: Direction of the line vortex: +1 for up, -1 for down.
 */
void add_line_in_zdir(struct tangle_state *tangle, double x, double y, int direction)
{
    // Boundary points.
    double zmin = tangle->box.bottom_left_back.p[2];
    double zmax = tangle->box.top_right_front.p[2];

    // Number of points.
    int points = (int) (2 * (zmax-zmin) / (global_dl_max - global_dl_min));
    while (number_of_empty_nodes(tangle) < points) expand_tangle(tangle, 2 * tangle->N);

    // Distance between points.
    double dz = (zmax - zmin) / points;

    // Orientation of the line vortex.
    double zstart = direction > 0 ? zmin : zmax;
    double zend = direction > 0 ? zmax : zmin;

    // Add first point.
    struct vec3 s = vec3(x, y, zstart);

    int new_pt = get_next_empty_node(tangle);
    int last_pt;
    tangle->vnodes[new_pt] = s;
    tangle->status[new_pt].status = pin_mode;
    tangle->status[new_pt].pin_wall = direction > 0 ? DOWN : UP;
    tangle->connections[new_pt].reverse = -1;

    // Loop through all the points, expet first and last, and add them.
    for(int k=1; k<points-1; ++k)
    {
        last_pt = new_pt;
        new_pt = get_next_empty_node(tangle);
        s.p[2] = zstart + direction*k*dz;

        tangle->vnodes[new_pt] = s;
        tangle->status[new_pt].status = FREE;
        tangle->status[new_pt].pin_wall = NOT_A_FACE;
        tangle->connections[new_pt].reverse = last_pt;
        tangle->connections[last_pt].forward = new_pt;
    }

    // Add last point.
    s.p[2] = zend;
    last_pt = new_pt;
    new_pt = get_next_empty_node(tangle);
    tangle->vnodes[new_pt] = s;
    tangle->status[new_pt].status = pin_mode;
    tangle->status[new_pt].pin_wall = direction > 0 ? UP : DOWN;
    tangle->connections[new_pt].reverse = last_pt;
    tangle->connections[new_pt].forward = -1;
    tangle->connections[last_pt].forward = new_pt;
}

/**
    Add straight line vortex given by segment (two vectors giving the endpoints of line).
    @param tangle: Tangle structure holding all the information about the vorices.
    @param seg: Segment specifying the line vortex to add.
    @param pin_wall_1: Wall on which the starting point (seg->r1) is pinned on.
    @param pin_wall_2: Wall on which the ending point (seg->r2) is pinned on.
 */
void add_line(struct tangle_state *tangle, const struct segment *seg, int pin_wall_1, int pin_wall_2)
{
    // Number of points.
    int points = 2 * segment_len(seg) / (global_dl_min + global_dl_max);
    while (number_of_empty_nodes(tangle) < points) expand_tangle(tangle, 2 * tangle->N);

    // Distance between points.
    double dx = (global_dl_min + global_dl_max) / 2;

    // Orientation of the line vortex.
    struct vec3 shift, dir;
    dir = segment_to_vec(seg);
    vec3_normalize(&dir);

    // Add first point.
    int new_pt = get_next_empty_node(tangle);
    int last_pt;
    vec3_add(tangle->vnodes+new_pt, &seg->r1, &VEC_NULL);
    tangle->status[new_pt].status = pin_mode;
    tangle->status[new_pt].pin_wall = pin_wall_1;
    tangle->connections[new_pt].reverse = -1;

    // Loop through all the points, expet first and last, and add them.
    for(int k = 1; k < points-1; ++k)
    {
        last_pt = new_pt;
        new_pt = get_next_empty_node(tangle);
        vec3_mul(&shift, &dir, dx*k);
        vec3_add(tangle->vnodes+new_pt, &seg->r1, &shift);
        tangle->status[new_pt].status = FREE;
        tangle->connections[new_pt].reverse = last_pt;
        tangle->connections[last_pt].forward = new_pt;
    }

    // Add last point.
    last_pt = new_pt;
    new_pt = get_next_empty_node(tangle);
    vec3_add(tangle->vnodes+new_pt, &seg->r2, &VEC_NULL);
    tangle->status[new_pt].status = pin_mode;
    tangle->status[new_pt].pin_wall = pin_wall_2;
    tangle->connections[new_pt].reverse = last_pt;
    tangle->connections[last_pt].forward = new_pt;
    tangle->connections[new_pt].forward = -1;
}

/**
    Add straight line vortex given by segment (two vectors giving the endpoints of line) and with helical osilation of given wavelength and amplitude.
    @param tangle: Tangle structure holding all the information about the vorices.
    @param seg: Segment specifying the line vortex to add.
    @param pin_wall_1: Wall on which the starting point (seg->r1) is pinned on.
    @param pin_wall_2: Wall on which the ending point (seg->r2) is pinned on.
    @param amplitude: Amplitude of the helical oscilation.
    @param wavelength: Wavelength of the helical oscilation.
 */
void add_line_with_oscilation(struct tangle_state *tangle, const struct segment *seg, int pin_wall_1, int pin_wall_2, double amplitude, double wavelength)
{
    // Number of points.
    int points = 2 * segment_len(seg) / (global_dl_min + global_dl_max);
    while (number_of_empty_nodes(tangle) < points) expand_tangle(tangle, 2 * tangle->N);

    // Distance between points.
    double dx = (global_dl_min + global_dl_max) / 2;

    // Orientation of the line vortex.
    struct vec3 shift, dir, u, v, deviation_u, deviation_v;
    struct segment position_seg;
    position_seg.r1 = seg->r1;
    dir = segment_to_vec(seg);
    vec3_normalize(&dir);

    // Create a base for helical oscilation vectors.
    struct vec3 pdir = perpendicular(&dir);
    vec3_cross(&u, &pdir, &dir);
    vec3_cross(&v, &dir, &u);
    vec3_normalize(&u);
    vec3_normalize(&v);

    // Prepare variables.
    int new_pt = get_next_empty_node(tangle);
    int last_pt;

    // Add first point into the tangle.
    tangle->vnodes[new_pt] = seg->r1;

    // Add oscilation to the first point.
    double phi = 0;
    double c = cos(phi);
    double s = sin(phi);

    vec3_mul(&deviation_u, &u, amplitude * c);
    vec3_mul(&deviation_v, &v, amplitude * s);

    vec3_add(tangle->vnodes + new_pt, tangle->vnodes + new_pt, &deviation_u);
    vec3_add(tangle->vnodes + new_pt, tangle->vnodes + new_pt, &deviation_v);

    // Connect first point.
    tangle->status[new_pt].status = pin_mode;
    tangle->status[new_pt].pin_wall = pin_wall_1;
    tangle->connections[new_pt].reverse = -1;

    // Loop through all the points, expet first and last, and add them.
    for(int k = 1; k < points-1; ++k)
    {
        // Add k-th point into the tangle.
        last_pt = new_pt;
        new_pt = get_next_empty_node(tangle);
        vec3_mul(&shift, &dir, dx*k);
        vec3_add(tangle->vnodes+new_pt, &seg->r1, &shift);

        // Add oscilation to the k-th point.
        position_seg.r2 = tangle->vnodes[new_pt];

        phi = 2 * M_PI * segment_len(&position_seg) / wavelength;
        c = cos(phi);
        s = sin(phi);

        vec3_mul(&deviation_u, &u, amplitude * c);
        vec3_mul(&deviation_v, &v, amplitude * s);

        vec3_add(tangle->vnodes + new_pt, tangle->vnodes + new_pt, &deviation_u);
        vec3_add(tangle->vnodes + new_pt, tangle->vnodes + new_pt, &deviation_v);

        // Connect k-th point.
        tangle->status[new_pt].status = FREE;
        tangle->connections[new_pt].reverse = last_pt;
        tangle->connections[last_pt].forward = new_pt;
    }

    // Add last point into the tangle.
    last_pt = new_pt;
    new_pt = get_next_empty_node(tangle);
    vec3_add(tangle->vnodes+new_pt, &seg->r2, &VEC_NULL);

    // Add oscilation to the last point.
    phi = 2 * M_PI * segment_len(seg) / wavelength;
    c = cos(phi);
    s = sin(phi);

    vec3_mul(&deviation_u, &u, amplitude * c);
    vec3_mul(&deviation_v, &v, amplitude * s);

    vec3_add(tangle->vnodes + new_pt, tangle->vnodes + new_pt, &deviation_u);
    vec3_add(tangle->vnodes + new_pt, tangle->vnodes + new_pt, &deviation_v);

    // Connect last point.
    tangle->status[new_pt].status = pin_mode;
    tangle->status[new_pt].pin_wall = pin_wall_2;
    tangle->connections[new_pt].reverse = last_pt;
    tangle->connections[last_pt].forward = new_pt;
    tangle->connections[new_pt].forward = -1;
}


/**
    Add straight line vortex given by segment (two vectors giving the endpoints of line) and with two helical osilations of given wavelength and amplitude.
    @param tangle: Tangle structure holding all the information about the vorices.
    @param seg: Segment specifying the line vortex to add.
    @param pin_wall_1: Wall on which the starting point (seg->r1) is pinned on.
    @param pin_wall_2: Wall on which the ending point (seg->r2) is pinned on.
    @param amplitude_1: Amplitude of the first helical oscilation.
    @param wavelength_1: Wavelength of the first helical oscilation.
    @param amplitude_2: Amplitude of the second helical oscilation.
    @param wavelength_2: Wavelength of the second helical oscilation.
 */
void add_line_with_two_oscilations(struct tangle_state* tangle, const struct segment* seg, int pin_wall_1, int pin_wall_2, double amplitude_1, double wavelength_1, double amplitude_2, double wavelength_2)
{
    // Number of points.
    int points = 2 * segment_len(seg) / (global_dl_min + global_dl_max);
    while (number_of_empty_nodes(tangle) < points) expand_tangle(tangle, 2 * tangle->N);

    // Distance between points.
    double dx = (global_dl_min + global_dl_max) / 2;

    // Orientation of the line vortex.
    struct vec3 shift, dir, u, v, deviation_u, deviation_v;
    struct segment position_seg;
    position_seg.r1 = seg->r1;
    dir = segment_to_vec(seg);
    vec3_normalize(&dir);

    // Create a base for helical oscilation vectors.
    struct vec3 pdir = perpendicular(&dir);
    vec3_cross(&u, &pdir, &dir);
    vec3_cross(&v, &dir, &u);
    vec3_normalize(&u);
    vec3_normalize(&v);

    // Prepare variables.
    int new_pt = get_next_empty_node(tangle);
    int last_pt;

    // Add first point into the tangle.
    tangle->vnodes[new_pt] = seg->r1;

    // Add oscilations to the first point.
    double phi_1 = 0;
    double phi_2 = 0;
    double c = cos(phi_1) * amplitude_1 + cos(phi_2) * amplitude_2;
    double s = sin(phi_1) * amplitude_1 + sin(phi_2) * amplitude_2;

    vec3_mul(&deviation_u, &u, c);
    vec3_mul(&deviation_v, &v, s);

    vec3_add(tangle->vnodes + new_pt, tangle->vnodes + new_pt, &deviation_u);
    vec3_add(tangle->vnodes + new_pt, tangle->vnodes + new_pt, &deviation_v);

    // Connect first point.
    tangle->status[new_pt].status = pin_mode;
    tangle->status[new_pt].pin_wall = pin_wall_1;
    tangle->connections[new_pt].reverse = -1;

    // Loop through all the points, expet first and last, and add them.
    for (int k = 1; k < points - 1; ++k)
    {
        // Add k-th point into the tangle.
        last_pt = new_pt;
        new_pt = get_next_empty_node(tangle);
        vec3_mul(&shift, &dir, dx * k);
        vec3_add(tangle->vnodes + new_pt, &seg->r1, &shift);

        // Add oscilations to the k-th point.
        position_seg.r2 = tangle->vnodes[new_pt];

        phi_1 = 2 * M_PI * segment_len(&position_seg) / wavelength_1;
        phi_2 = 2 * M_PI * segment_len(&position_seg) / wavelength_2;

        c = cos(phi_1) * amplitude_1 + cos(phi_2) * amplitude_2;
        s = sin(phi_1) * amplitude_1 + sin(phi_2) * amplitude_2;

        vec3_mul(&deviation_u, &u, c);
        vec3_mul(&deviation_v, &v, s);

        vec3_add(tangle->vnodes + new_pt, tangle->vnodes + new_pt, &deviation_u);
        vec3_add(tangle->vnodes + new_pt, tangle->vnodes + new_pt, &deviation_v);

        // Connect k-th point.
        tangle->status[new_pt].status = FREE;
        tangle->connections[new_pt].reverse = last_pt;
        tangle->connections[last_pt].forward = new_pt;
    }

    // Add last point into the tangle.
    last_pt = new_pt;
    new_pt = get_next_empty_node(tangle);
    vec3_add(tangle->vnodes + new_pt, &seg->r2, &VEC_NULL);

    // Add oscilations to the last point.
    phi_1 = 2 * M_PI * segment_len(seg) / wavelength_1;
    phi_2 = 2 * M_PI * segment_len(seg) / wavelength_2;

    c = cos(phi_1) * amplitude_1 + cos(phi_2) * amplitude_2;
    s = sin(phi_1) * amplitude_1 + sin(phi_2) * amplitude_2;

    vec3_mul(&deviation_u, &u, c);
    vec3_mul(&deviation_v, &v, s);

    vec3_add(tangle->vnodes + new_pt, tangle->vnodes + new_pt, &deviation_u);
    vec3_add(tangle->vnodes + new_pt, tangle->vnodes + new_pt, &deviation_v);

    // Connect last point.
    tangle->status[new_pt].status = pin_mode;
    tangle->status[new_pt].pin_wall = pin_wall_2;
    tangle->connections[new_pt].reverse = last_pt;
    tangle->connections[last_pt].forward = new_pt;
    tangle->connections[new_pt].forward = -1;
}

/**
    Add straight line vortex vortex in z direction at position x,y and direction up (+1) or down (-1), 
    and with two helical osilations of given wavelength and amplitude. Added line is not clipped to the wall, but rather is periodic and used in periodic boundaries.
    @param tangle: Tangle structure holding all the information about the vorices.
    @param x: Position of the points in the x direction.
    @param y: Position of the points in the y direction.
    @param direction: Direction of the line vortex: +1 for up, -1 for down.
    @param amplitude_1: Amplitude of the first helical oscilation.
    @param wavelength_1: Wavelength of the first helical oscilation.
    @param amplitude_2: Amplitude of the second helical oscilation.
    @param wavelength_2: Wavelength of the second helical oscilation.
 */
void add_line_in_zdir_with_two_oscilations_with_periodic_boundary(struct tangle_state *tangle, double x, double y, int direction, double amplitude_1, int number_of_waves_1, double amplitude_2, int number_of_waves_2)
{
    // Boundary points.
    double zmin = tangle->box.bottom_left_back.p[2];
    double zmax = tangle->box.top_right_front.p[2];

    // Number of points.
    int points = (int)(2 * (zmax - zmin) / (global_dl_max - global_dl_min));
    while (number_of_empty_nodes(tangle) < points) expand_tangle(tangle, 2 * tangle->N);

    // Distance between points.
    double dz = (zmax - zmin) / points;

    // Orientation of the line vortex.
    double zstart = direction > 0 ? zmin : zmax;

    // Define variables.
    struct vec3 shift, dir, start, deviation_1, deviation_2;
    dir = vec3(0.0,0.0,direction);
    start = vec3(x, y, zstart);

    int first_pt = get_next_empty_node(tangle);
    int last_pt;
    int new_pt = first_pt;

    // Add first point and its oscillations.
    tangle->vnodes[new_pt] = start;
    deviation_1 = vec3(amplitude_1, 0, 0);
    deviation_2 = vec3(amplitude_2, 0, 0);

    vec3_add(tangle->vnodes+new_pt, tangle->vnodes+new_pt, &deviation_1);
    vec3_add(tangle->vnodes+new_pt, tangle->vnodes+new_pt, &deviation_2);
    tangle->status[new_pt].status = FREE;

    // Loop through all the points, expet first, and add them.
    for(int k = 1; k < points-1; ++k)
    {
        // Add k-th point.
        last_pt = new_pt;
        new_pt = get_next_empty_node(tangle);
        vec3_mul(&shift, &dir, dz*k);
        vec3_add(&shift, &start, &shift);
        tangle->vnodes[new_pt] = shift;

        // Add oscilations to k-th point.
        deviation_1 = vec3(amplitude_1 * cos(2 * M_PI * k * number_of_waves_1 / points), -amplitude_1 * sin(2 * M_PI * k * number_of_waves_1 / points), 0);
        deviation_2 = vec3(amplitude_2 * cos(2 * M_PI * k * number_of_waves_2 / points), -amplitude_2 * sin(2 * M_PI * k * number_of_waves_2 / points), 0);

        vec3_add(tangle->vnodes+new_pt, tangle->vnodes+new_pt, &deviation_1);
        vec3_add(tangle->vnodes+new_pt, tangle->vnodes+new_pt, &deviation_2);

        // Connect k-th point.
        tangle->status[new_pt].status = FREE;
        tangle->connections[new_pt].reverse = last_pt;
        tangle->connections[last_pt].forward = new_pt;
    }

    // Connect first and last points.
    tangle->connections[new_pt].forward = first_pt;
    tangle->connections[first_pt].reverse = new_pt;
}

/**
    Add N random straight line vortices in z direction and random orientation.
    @param tangle: Tangle structure holding all the information about the vorices.
    @param N: Number of line vortices to add.
*/
void random_straight_lines_in_zdir(struct tangle_state *tangle, int N)
{
    // Finds the boundaries of the box to randomise poitions of the straight vortex lines.
    double xmin = tangle->box.bottom_left_back.p[0];
    double xmax = tangle->box.top_right_front.p[0];
    double ymin = tangle->box.bottom_left_back.p[1];
    double ymax = tangle->box.top_right_front.p[1];

    // Loop over all vortex lines.
    for(int k=0; k < N; ++k)
    {
        // Randomise poitions of the straight vortex line.
        double x1 = xmin + (xmax - xmin)*drand48();
        double y1 = ymin + (ymax - ymin)*drand48();
        double x2 = xmin + (xmax - xmin)*drand48();
        double y2 = ymin + (ymax - ymin)*drand48();

        // Add line vortex into the tangle.
        add_line_in_zdir(tangle, x1, y1, +1);
        add_line_in_zdir(tangle, x2, y2, -1);
    }
}

/**
    Main function for injecting vortices. Parmeters of injected vortices are set up in initial_conditions.c.
    This function calls inject_loop(tangle, t, loop_injection_frequency) and inject_line_pairs(tangle, t, line_injection_frequency, line_injection_n, line_injection_polarized).
    @param tangle: Tangle structure holding all the information about the vorices.
    @param t: Time of the simulation evolution.
*/
void inject_vortices(struct tangle_state *tangle, double t)
{
    // Inject vortices.
    if(loop_injection) inject_loop(tangle, t, loop_injection_frequency);
    if(line_injection) inject_line_pairs(tangle, t, line_injection_frequency, line_injection_n, line_injection_polarized);
}

/**
    Injects random loops placed randomly in the top plane of the domain box with given frequency.
    @param tangle: Tangle structure holding all the information about the vorices.
    @param t: Time of the simulation evolution.
    @param frequency: Frequency of injection of the loop vortices.
*/
void inject_loop(struct tangle_state *tangle, double t, double frequency)
{
    static double last_injection = 0;

    // Direction and center of the injected loop inject in random direction pointing down-ish.
    struct vec3 dir = vec3(2*(drand48() - 0.5), 2*(drand48()-0.5), -drand48());
    vec3_normalize(&dir);
    struct vec3 cent;

    // Is it time yet to inject a loop?
    if(t - last_injection < 1/frequency) return;
    last_injection = t;

    // Boundaries of the box, where to inject loops.
    double left, right, back, front;
    left = tangle->box.bottom_left_back.p[0];
    right = tangle->box.top_right_front.p[0];

    back = tangle->box.bottom_left_back.p[1];
    front = tangle->box.top_right_front.p[1];

    // The injected loop is placed randomly in the top plane of the domain box.
    cent.p[2] = tangle->box.top_right_front.p[2];
    cent.p[0] = left + (right - left)*drand48();
    cent.p[1] = back + (front - back)*drand48();

    // Randomize radius of the injected loop.
    double D = fabs(right - left);
    double rmin = 0.05*D;
    double rmax = 0.25*D;
	
	double r = rmin + (rmax - rmin)*drand48();

    // Add loop.
    add_loop(tangle, &cent, &dir, r);
}

/**
Injects a pair of straight vortex lines of opposite circulations oriented along z-axis.
   * Only makes sense for the 2-wall and periodic boundaries.	@param tangle: Tangle structure holding all the information about the vorices.
	@param t: Time of the simulation evolution.
	@param frequency: Frequency of injection of the loop vortices.
	@param line_injection_n: Number of line injected.
	@param polarized: Polarization 
*/
void inject_line_pairs(struct tangle_state *tangle, double t, double frequency, int line_injection_n, int polarized)
{
	static double last_injection = 0;

	// Is it time yet to inject the lines?
	if(t - last_injection < 1/frequency) return;
	last_injection = t;
	
	// Loop through the all line vortex pairs to inject.
	for(int k = 0; k < line_injection_n; ++k)
    {
		// Boundaries of the box, where to inject line vortices.
		double left = tangle->box.bottom_left_back.p[0];
		double right = tangle->box.top_right_front.p[0];

		double back = tangle->box.bottom_left_back.p[1];
		double front = tangle->box.top_right_front.p[1];

		// Positions of the vortex line pair.
		double x1, y1, x2, y2;

		x1 = left + drand48()*(right-left);
		x2 = left + drand48()*(right-left);

		// Get random positions for both lines.
		if(polarized)
		{
			y1 = back + (1 + drand48())*(front-back)/2; // Positive vortex in (0.5, 1).
			y2 = back + drand48()*(front-back)/2; // Negative vortex in (0, 0.5).
		}
		else // Unpolarized, fully random positions.
		{
			y1 = back + drand48()*(front-back);
			y2 = back + drand48()*(front-back);
		}

		// Add both lines.
		add_line_in_zdir(tangle, x1, y1, +1);
		add_line_in_zdir(tangle, x2, y2, -1);
    }
}
