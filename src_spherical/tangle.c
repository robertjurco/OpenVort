#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <assert.h>
#include <stdio.h>

#include <external_velocity.h>
#include "tangle.h"
#include "vortex_constants.h"
#include "vec3_math.h"
#include "boundaries.h"
#include "octree.h"

#define _GNU_SOURCE
#include <fenv.h>

/*******************************************************************************
 *************************** TANGLE BASIC FUNCTIONS ****************************
 ******************************************************************************/

 /*
	Initialize the tangle structue with the given size.
	@param tangle: Tangle to initialize.
	@param n: Total number of nodes inside the tangle.
 */
void create_tangle(struct tangle_state *tangle, size_t n)
{
	tangle->vnodes      = (struct vec3*)malloc(sizeof(struct vec3)*n);
	tangle->vels        = (struct vec3*)malloc(sizeof(struct vec3)*n);
	tangle->vs          = (struct vec3*)malloc(sizeof(struct vec3)*n);
	tangle->tangents    = (struct vec3*)malloc(sizeof(struct vec3)*n);
	tangle->normals     = (struct vec3*)malloc(sizeof(struct vec3)*n);

	tangle->recalculate = (int*)malloc(sizeof(int)*n);
	tangle->status      = (struct node_status*)malloc(sizeof(struct node_status)*n);

	tangle->connections = (struct neighbour_t*)malloc(sizeof(struct neighbour_t)*n);

	// All the nodes in the tangle are EMPTY.
	for(size_t k=0; k<n; ++k)
    {
		tangle->status[k].status = EMPTY;
		tangle->status[k].pin_surface = -1;
		tangle->connections[k].forward = -1;
		tangle->connections[k].reverse = -1;
    }

	// Default initialisation is to open bounadry conditions.
	struct domain domain_section = {
		.inner_radius = 0.02,
		.outer_radius = 0.1,
		.azimut_angle = DEG2RAD(15),
		.polar_angle = DEG2RAD(15)
	};

	tangle->domain_section = domain_section;

	tangle->N = n;
	tangle->next_empty = 0;
	tangle->total_empty = n;
}

/*
	Expant the tangle structue into the given size.
	@param tangle: Tangle to expant.
	@param n: Total number of nodes inside the expanded tangle.
 */
void expand_tangle(struct tangle_state *tangle, size_t n)
{
    size_t k;
    size_t old_n = tangle->N;
	
	tangle->vnodes      = (struct vec3*)realloc(tangle->vnodes, sizeof(struct vec3)*n);
    tangle->vels        = (struct vec3*)realloc(tangle->vels, sizeof(struct vec3)*n);
    tangle->vs          = (struct vec3*)realloc(tangle->vs, sizeof(struct vec3)*n);
    tangle->tangents    = (struct vec3*)realloc(tangle->tangents, sizeof(struct vec3)*n);
    tangle->normals     = (struct vec3*)realloc(tangle->normals, sizeof(struct vec3)*n);

    tangle->recalculate = (int*)realloc(tangle->recalculate, sizeof(int)*n);
    tangle->status      = (struct node_status*)realloc(tangle->status, sizeof(struct node_status)*n);

    tangle->connections = (struct neighbour_t*)realloc(tangle->connections, n*sizeof(struct neighbour_t));

	// All the new nodes added into the tangle are EMPTY.
    if(old_n < n)
    {
        for(k = old_n; k < n; ++k)
        {
            tangle->status[k].status = EMPTY;
            tangle->status[k].pin_surface = -1;
            tangle->connections[k].forward = -1;
            tangle->connections[k].reverse = -1;
        }
    }

	tangle->N = n;
	tangle->total_empty += n - old_n;

	// This has to go last, as get_tangle_next_free(tangle) is using new tange informations AND can use expand_tangle.
    tangle->next_empty = get_next_empty_node(tangle);
}

/*
	Frees the tangle structue.
	@param tangle: Tangle to free.
 */
void free_tangle(struct tangle_state *tangle)
{
	free(tangle->vnodes);
	free(tangle->vels);
	free(tangle->vs);
	free(tangle->tangents);
	free(tangle->normals);
	free(tangle->connections);
	free(tangle->recalculate);
	free(tangle->status);
}

/*******************************************************************************
 ************************** TANGLE POINTS FUNCTIONS ****************************
 ******************************************************************************/

 /*
	Returns the number of unused (EMPTY) nodes inside the tangle.
	@param tangle: Tangle structure holding all the informations about the vortices.
	@returns Returns the number of unused (EMPTY) nodes inside the tangle.
 */
int number_of_empty_nodes(struct tangle_state *tangle)
{
	int sum = 0;
	for (int k = 0; k < tangle->N; ++k) if (tangle->status[k].status == EMPTY) sum++;
	return sum;
}

/*
   Returns the number of used (FREE, PINNED or PINNED_SLIP) nodes inside the tangle.
   @param tangle: Tangle structure holding all the informations about the vortices.
   @returns Returns the number of used (FREE, PINNED or PINNED_SLIP) nodes inside the tangle.
*/
int number_of_used_nodes(struct tangle_state *tangle)
{
	int sum = 0;
	for(int k = 0; k < tangle->N; ++k) if (tangle->status[k].status != EMPTY) sum++;
	return sum;
}

/*
   Returns the position of the node given number of the steps away from the given node.
   In the case of wall, we continue stepping throught the mirror image on the other side.
   @param tangle: Tangle structure holding all the informations about the vortices.
   @param i: Id of the point we want to make number of steps away from.
   @param where: Number of steps we wwant to make. Sign is determining the direction, +1 is forward and -1 is reverse.
   @returns Returns the position of the node given number of the steps away from the given node.
*/
struct vec3 step_node(const struct tangle_state *tangle, int i, int where)
{
	// i-th node is not empty.
	assert(tangle->status[i].status != EMPTY);

	// If where = 0 returns the current position.
	if (where == 0) return tangle->vnodes[i];

	// If node is free, recursively calls itself.
	if (tangle->status[i].status == FREE)
	{
		if (where > 0) return step_node(tangle, tangle->connections[i].forward, where - 1);
		else if (where < 0)	return step_node(tangle, tangle->connections[i].reverse, where + 1);
	}
	// If the node is pinned, walk the rest of the steps back.
	else if(tangle->status[i].status == PINNED || tangle->status[i].status == PINNED_SLIP)
	{
		struct vec3 out;
		if(where > 0)
		{
			if(tangle->connections[i].forward < 0) out = step_node(tangle, i, -where);
			else return step_node(tangle, tangle->connections[i].forward, where-1); 
		}
  		else if(where < 0)
		{
	  		if(tangle->connections[i].reverse < 0) out = step_node(tangle, i, -where);
	  		else return step_node(tangle, tangle->connections[i].reverse, where+1);
		}

		// We remembered to inverse cause we walked through the surface.
		if (tangle->status[i].pin_surface == INNER) return spherical_inversion(&out, tangle->domain_section.inner_radius);
		if (tangle->status[i].pin_surface == OUTER) return spherical_inversion(&out, tangle->domain_section.outer_radius);
	}
	else
	{
		printf("Walking across empty node."); // We should never get here.
		fflush(stdout);
		return vec3(0, 0, 0);
	}

	// To suppress warning.
	printf("ERROR in step_node().\n");
	printf("dfs %d, %d, %d\n", where, i, tangle->status[i].status);

	fflush(stdout);
	return vec3(0,0,0);
}

/*
   Returns ID of the next empty node. In the case of no empty nodes inside the tangle, returns -1. This function is called by get_next_empty_node(tangle).
   @param tangle: Tangle structure holding all the informations about the vortices.
   @returns Returns ID of the next empty node. Returns -1 in the case of no empty nodes inside the tangle.
*/
static inline int search_next_empty_node(struct tangle_state *tangle)
{
	// If something is already available, just return that.
	if (tangle->next_empty > 0) {
		if (tangle->next_empty < tangle->N && tangle->status[tangle->next_empty].status == EMPTY)
		{
			int return_id = tangle->next_empty;

			// If next node is empty, remember it. Otherwise set tangle->next_empty to -1;
			if (tangle->next_empty + 1 < tangle->N && tangle->status[tangle->next_empty + 1].status == EMPTY) ++tangle->next_empty;
			else tangle->next_empty = -1;

			// Returns ID of the empty node.
			return return_id;
		}
		// If something is wrong, if there is a non-empty node in tangle->next_free.
		else printf("ERROR: There is a non-empty point in tangle->next_free.");
	}
	
	// Otherwise we have to search for it.
	for (int k = 0; k < tangle->N; ++k) {
		if(tangle->status[k].status == EMPTY)
		{
			// If next node is empty, remember it. Otherwise set tangle->next_empty to -1;
			if (tangle->status[k + 1].status == EMPTY) tangle->next_empty = k + 1;
			else tangle->next_empty = -1;

			return k;
		}
	}

	// We haven't found anything, return -1.
	return -1;
}

/*
   Returns ID of the next empty node. Expands the tangle if needed. Calls search_next_empty_node(tangle).
   @param tangle: Tangle structure holding all the informations about the vortices.
   @returns Returns ID of the next empty node.
*/
int get_next_empty_node(struct tangle_state *tangle)
{
	int idx = search_next_empty_node(tangle);

	// Expands the tangle if needed.
	while(idx < 0)
	{
		expand_tangle(tangle, 2*tangle->N);
		idx = search_next_empty_node(tangle);
	}

	return idx;
}

// Needed for add point.
struct vec3 normal_protate(struct tangle_state* tangle, int p, int next, const struct domain* domain_section)
{
	// Angles of the domain.
	double domain_azimut_angle = 2 * domain_section->azimut_angle;
	double domain_polar_angle = 2 * domain_section->polar_angle;

	// Vectors and normals
	struct vec3 v1 = tangle->vnodes[p];
	struct vec3 v2 = tangle->vnodes[next];
	struct vec3 n2 = tangle->normals[next];

	// Don't forget that it may be needed to wrap one point multiple times periodically.

	// Check azimutal angles.
	double azimut_v1 = azimut_angle(&v1);
	double azimut_v2 = azimut_angle(&v2);
	double azimut_n2 = azimut_angle(&n2);
	while (azimut_v2 - azimut_v1 > domain_azimut_angle)
	{
		azimut_v2 -= domain_azimut_angle;
		azimut_n2 -= domain_azimut_angle;
	}
	while (azimut_v2 - azimut_v1 < -domain_azimut_angle)
	{
		azimut_v2 += domain_azimut_angle;
		azimut_n2 += domain_azimut_angle;
	}

	// Check polar angles.
	double polar_v1 = polar_angle(&v1);
	double polar_v2 = polar_angle(&v2);
	double polar_n2 = polar_angle(&n2);
	while (polar_v2 - polar_v1 > domain_polar_angle)
	{
		polar_v2 -= domain_polar_angle;
		polar_n2 -= domain_polar_angle;
	}
	while (polar_v2 - polar_v1 < -domain_polar_angle)
	{
		polar_v2 += domain_polar_angle;
		polar_n2 += domain_polar_angle;
	}

	double rad_n2 = radius(&n2);

	return spherical_to_vector(rad_n2, azimut_n2, polar_n2);
}

/*
	Add a point between p and p+1 (p+1 in the sense of connections) into the tangle.
	@param tangle: Tangle structure to which we want the point to add.
	@param p: Id of the point after which we want the point to add.
	@returns Returns the id of the added point.
*/
int add_point(struct tangle_state *tangle, int p)
{
	// Find and initialize new point.
	int new_pt = get_next_empty_node(tangle);
	tangle->status[new_pt].status = FREE;
	tangle->status[new_pt].pin_surface = NOT_A_SURFACE;

	// Find next (p+1) point.
	int next = tangle->connections[p].forward;

	// Update tangents and normals of points p and p+1.
	update_tangent_normal(tangle, p);
	update_tangent_normal(tangle, next);
	// Also rotate normals.
	struct vec3 s0pp = tangle->normals[p];
	struct vec3 s1pp = normal_protate(tangle, p, next, &tangle->domain_section);

	// Find positions of points p and p+1, pwraps it.
	struct vec3 s0 = tangle->vnodes[p];
	struct vec3 s1 = tangle->vnodes[next];
	struct segment seg = seg_pwrap(&s0, &s1, &tangle->domain_section);
	s1 = seg.r2;

	// And their distance.
	double l = vec3_dist(&s0, &s1);

	// Define variables.
	struct vec3 a, b, new, n;

	// Normal vector is zero if there is no curvature, hence we find out the curvature at the new point by avarage.
	vec3_add(&n, &s0pp, &s1pp);
	vec3_mul(&n, &n, 0.5);
	
	if(vec3_len(&n) > 1e-5) // n will be identically 0 for a straight vortex.
	{
		// The idea is to fit the line between p and p+1 by a circle with the right curvature.
		// dR is the squared distance from centre of the circle into the mid-point A of the straight line connecting p and p+1.
		// delta is the distance from A to the circle (delta = R - sqrt(dR)).
		// New point will sit in the middle of the circle arc bounded by p and p+1.

		double R = 1/vec3_len(&n);
		double dR = R*R - l*l/4;
		// dR can become < 0 for sharp cusps.
		// Simplest way to deal with it is to treat the s0 and s1 as sitting on opposite ends of a circle, for which dR = 0.
		// This does not preserve curvature, but this is below our resolution anyway.
		double delta = dR > 0 ? R - sqrt(dR) : R;

		// Normal is pointing into the circle, we needs it to point out.
		vec3_normalize(&n);
		vec3_mul(&n, &n, -1);

		// The mid-point A of the straight line connecting p and p+1.
		vec3_add(&a, &s0, &s1);
		vec3_mul(&a, &a, 0.5);

		// Vector from A to the circle (delta = R - sqrt(dR)).
		vec3_mul(&b, &n, delta);

		// New point sits in the middle of the circle arc bounded by p and p+1.
		vec3_add(&new, &a, &b);
	}
	else // We basically have a straight vortex, just average s0 and s1.
	{
		vec3_add(&new, &s0, &s1);
		vec3_mul(&new, &new, 0.5);
	}

	tangle->vnodes[new_pt] = new;
	tangle->connections[new_pt].reverse = p;
	tangle->connections[new_pt].forward = next;
	tangle->connections[p].forward = new_pt;
	tangle->connections[next].reverse = new_pt;

	return new_pt;
}

/*
	Removes the point from the tangle.
	@param tangle: Tangle structure from which we want the point to remove.
	@param point_idx: Id of the point to remove.
*/
void remove_point(struct tangle_state *tangle, int point_idx)
{
	// Finds id of previous and next nodes.
	int prev = tangle->connections[point_idx].reverse;
	int next = tangle->connections[point_idx].forward;

	// If possible, connects the previous node with the next one.
	if(prev >= 0) tangle->connections[prev].forward = next;
	if(next >= 0) tangle->connections[next].reverse = prev;

	// Set the node EMPTY.
	tangle->connections[point_idx].reverse = tangle->connections[point_idx].forward = -1;
	tangle->status[point_idx].status = EMPTY;
}

/*******************************************************************************
 ********************** TANGLE FUNCTIONS FOR SIMULATION ************************
 ******************************************************************************/

 /*
	Updates the tangent and normal of k-th node of the tangle. Calculation for k-th node is done from 5 nodes (k-2, k-1, k, k+1, k+2).
	ADD THE NAME OF THE METHOD (I can't remember it now).
	@param tangle: Tangle structure holding all informations about vortices we want to update.
	@param k: ID of the node for which we want to update tangent and normal.
*/
void update_tangent_normal(struct tangle_state *tangle, size_t k)
{
	// If EMPTY return.
	if (tangle->status[k].status == EMPTY)
		return;

	// Define variables.
	struct vec3 s0, s1, sm1;
	struct vec3 s2, sm2;
	size_t i;

	// Vector differences.
	struct vec3 ds[4];
	struct segment dseg[4];
	struct segment dseg_12;
	struct segment dseg_m12;

	s0  = tangle->vnodes[k];
	s1 = step_node(tangle, k, 1);
	s2 = step_node(tangle, k, 2);
	sm1 = step_node(tangle, k, -1);
	sm2 = step_node(tangle, k, -2);

	dseg[0] = seg_pwrap(&s0, &s2, &tangle->domain_section);
	dseg[1] = seg_pwrap(&s0, &s1, &tangle->domain_section);
	dseg[2] = seg_pwrap(&s0, &sm1, &tangle->domain_section);
	dseg[3] = seg_pwrap(&s0, &sm2, &tangle->domain_section);
	dseg_12 = seg_pwrap(&s1, &s2, &tangle->domain_section);
	dseg_m12 = seg_pwrap(&sm1, &sm2, &tangle->domain_section);

	for(int j = 0; j<4; ++j) ds[j] = segment_to_vec(&dseg[j]);

	// Lengths.
	double d1 = segment_len(&dseg[1]);
	double d2 = d1 + segment_len(&dseg_12);
	double dm1 = segment_len(&dseg[2]);
	double dm2 = dm1 + segment_len(&dseg_m12);

	// Four point coefficients, denominators.
	double d_s_diff[] = {
		d2*(d2 - d1)*(dm1 + d2)*(dm2 + d2),
		d1*(d2 - d1)*(dm1 + d1)*(dm2 + d1),
		dm1*(dm1 + d1)*(dm1 + d2)*(dm2 - dm1),
		dm2*(dm2 + d1)*(dm2 + d2)*(dm2 - dm1)
	};

  // First derivative, four point coefficients, nominators, O(d^4).
	double s_1_cf[] = {
		-d1*dm1*dm2,
		d2*dm1*dm2,
		-d1*d2*dm2,
		d1*d2*dm1
	};

	// Second derivative, four point coefficients, nominators, O(d^3).
	double s_2_cf[] = {
		2*((dm1 - d1)*dm2 - d1*dm1),
		-2*((dm1 - d2)*dm2 - d2*dm1),
		2*((d2  + d1)*dm2 - d1*d2),
		-2*((d2  + d1)*dm1 - d1*d2)
	};

	// Solve for tangent and normal.
	for(i=0; i<3; ++i)
	{
		tangle->tangents[k].p[i] = 0;
		tangle->normals[k].p[i]  = 0;
		for(int z = 0; z<4; ++z)
		{
	  		tangle->tangents[k].p[i] += s_1_cf[z]/d_s_diff[z]*ds[z].p[i];
	  		tangle->normals[k].p[i]  += s_2_cf[z]/d_s_diff[z]*ds[z].p[i];
		}
	}
}

/*
	Updates the tangents and normals of all nodes of the tangle. Calculation for k-th node is done from 5 nodes (k-2, k-1, k, k+1, k+2).
	ADD THE NAME OF THE METHOD (I can't remember it now).
	@param tangle: Tangle structure holding all informations about vortices we want to update.
*/
void update_tangents_normals(struct tangle_state* tangle)
{
	int i;
	#pragma omp parallel private(i) num_threads(global_num_threads)
	{
		#pragma omp for
		for (i = 0; i < tangle->N; ++i) update_tangent_normal(tangle, i);
	}
}


/*
	Calculates the velocity field of both segments in the vicinity of point i (at that point) using Local Induction Approximation.
	@param tangle: Tangle structure holding all informations about vortices.
	@param i: ID of node where we want the velocity field generated by neighbour segments to evaluate.
	@returns Returns the velocity field of both neighbour segments given by nodes i-1, i and i+1 at the position of the i-th node.
 */
static inline struct vec3 lia_velocity(const struct tangle_state *tangle, int i)
{
	// Find the position of neighbour nodes.
	const struct vec3 *p = tangle->vnodes + i;
	struct vec3 next = step_node(tangle, i, +1);
	struct vec3 prev = step_node(tangle, i, -1);

	// Pwrap those positions.
	struct segment sf = seg_pwrap(p, &next, &tangle->domain_section);
	struct segment sr = seg_pwrap(&prev, p, &tangle->domain_section);

	// Calculate the distance to the neighbour nodes.
	double l_next = segment_len(&sf);
	double l_prev = segment_len(&sr);

	// Equation for LIA approximation.
	double f = KAPPA*log(2*sqrt(l_next*l_prev)/sqrt(M_E)/VORTEX_WIDTH)/4/M_PI;

	// Calculate the direction of the velocity field at the i-th node induced by both its neighbourhood segments.
	struct vec3 vv;
	vec3_cross(&vv, &tangle->tangents[i], &tangle->normals[i]);
	vec3_mul(&vv, &vv, f);

	return vv;
}

struct vec3 shifted(const struct vec3* r, double shift_azimut_angle, double shift_polar_angle)
{
	double rad = radius(r);
	double azim = azimut_angle(r);
	double polar = polar_angle(r);

	return spherical_to_vector(rad, azim + shift_azimut_angle, polar + shift_polar_angle);
}

/*
	Updates the velocities of k-th node of the tangle.
	@param tangle: Tangle structure holding all informations about vortices we want to update.
	@param t: Time of the simulation when we are updating the tangle.
	@param tree: Pointer to the BH tree. If we are not using BH algortihm *tree = NULL;
*/
void update_velocity(struct tangle_state *tangle, struct tangle_state* inner_img, struct tangle_state* outer_img, int k, double t, struct octree *tree, struct octree* inner_img_tree, struct octree* outer_img_tree)
{
	// If the node is empty not needed to update it.
    if(tangle->status[k].status == EMPTY) return;

	// If the node is pinned, it has the velocity of the boundary.
    if(tangle->status[k].status == PINNED)
    {
		// Node by itself is not moving.
        tangle->vs[k] = vec3(0,0,0);
        // Get the boundary velocity, by default non-moving.
        get_vb(&tangle->vnodes[k], t, &tangle->vels[k]);
        return;
    }

	// Calculate the velocity induced by the vicinity of the node using LIA.
    tangle->vs[k] = lia_velocity(tangle, k);

	// Add the external velocity of the superfluid component.
    struct vec3 evs;
    get_vs(&tangle->vnodes[k], t, &evs);
    vec3_add(&tangle->vs[k], &tangle->vs[k], &evs);

	// Define variables.
    struct vec3 shift_r, v_shift, v_safe;
    struct vec3 v_shift_total = VEC_NULL;

	// Use the tree approximation to the full Biot-Savart.
	struct vec3 v_tree;
	octree_get_vs(tangle, tree, &tangle->vnodes[k], BH_resolution, &v_tree);
	vec3_add(&tangle->vs[k], &tangle->vs[k], &v_tree);
		
	// Calculate the velocity due to boundary images.
	
	// Angles of the domain.
	double domain_azimut_angle = 2 * tangle->domain_section.azimut_angle;
	double domain_polar_angle = 2 * tangle->domain_section.polar_angle;

	// Velocity from the periodic boundary image, do not forget to shift velocity back on the right place.

	shift_r = shifted(&tangle->vnodes[k], domain_azimut_angle, 0.0);
	octree_get_vs(tangle, tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, -domain_azimut_angle, 0.0);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], -domain_azimut_angle, 0.0);
	octree_get_vs(tangle, tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, domain_azimut_angle, 0.0);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], 0.0, domain_polar_angle);
	octree_get_vs(tangle, tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, 0.0, -domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], 0.0, -domain_polar_angle);
	octree_get_vs(tangle, tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, 0.0, domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	// Velocity from corner boundary images.

	shift_r = shifted(&tangle->vnodes[k], domain_azimut_angle, domain_polar_angle);
	octree_get_vs(tangle, tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, -domain_azimut_angle, -domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], -domain_azimut_angle, domain_polar_angle);
	octree_get_vs(tangle, tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, domain_azimut_angle, -domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], domain_azimut_angle, -domain_polar_angle);
	octree_get_vs(tangle, tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, -domain_azimut_angle, domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], -domain_azimut_angle, -domain_polar_angle);
	octree_get_vs(tangle, tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, domain_azimut_angle, domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	// The case of the mirror wall, no need to mirror reflet since we have already sperical inversion.
	octree_get_vs(inner_img, inner_img_tree, &tangle->vnodes[k], BH_resolution, &v_shift);
    vec3_add(&v_shift_total, &v_shift_total, &v_shift);
	octree_get_vs(outer_img, outer_img_tree, &tangle->vnodes[k], BH_resolution, &v_shift);
	vec3_add(&v_shift_total, &v_shift_total, &v_shift);

	// Inner wall sides.
	shift_r = shifted(&tangle->vnodes[k], domain_azimut_angle, 0.0);
	octree_get_vs(inner_img, inner_img_tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, -domain_azimut_angle, 0.0);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], -domain_azimut_angle, 0.0);
	octree_get_vs(inner_img, inner_img_tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, domain_azimut_angle, 0.0);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], 0.0, domain_polar_angle);
	octree_get_vs(inner_img, inner_img_tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, 0.0, -domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], 0.0, -domain_polar_angle);
	octree_get_vs(inner_img, inner_img_tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, 0.0, domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	// Inner wall corners.
	shift_r = shifted(&tangle->vnodes[k], domain_azimut_angle, domain_polar_angle);
	octree_get_vs(inner_img, inner_img_tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, -domain_azimut_angle, -domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], -domain_azimut_angle, domain_polar_angle);
	octree_get_vs(inner_img, inner_img_tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, domain_azimut_angle, -domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], domain_azimut_angle, -domain_polar_angle);
	octree_get_vs(inner_img, inner_img_tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, -domain_azimut_angle, domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], -domain_azimut_angle, -domain_polar_angle);
	octree_get_vs(inner_img, inner_img_tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, domain_azimut_angle, domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	// Outer wall sides.
	shift_r = shifted(&tangle->vnodes[k], domain_azimut_angle, 0.0);
	octree_get_vs(outer_img, outer_img_tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, -domain_azimut_angle, 0.0);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], -domain_azimut_angle, 0.0);
	octree_get_vs(outer_img, outer_img_tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, domain_azimut_angle, 0.0);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], 0.0, domain_polar_angle);
	octree_get_vs(outer_img, outer_img_tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, 0.0, -domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], 0.0, -domain_polar_angle);
	octree_get_vs(outer_img, outer_img_tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, 0.0, domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	// Outer wall corners.
	shift_r = shifted(&tangle->vnodes[k], domain_azimut_angle, domain_polar_angle);
	octree_get_vs(outer_img, outer_img_tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, -domain_azimut_angle, -domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], -domain_azimut_angle, domain_polar_angle);
	octree_get_vs(outer_img, outer_img_tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, domain_azimut_angle, -domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], domain_azimut_angle, -domain_polar_angle);
	octree_get_vs(outer_img, outer_img_tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, -domain_azimut_angle, domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);

	shift_r = shifted(&tangle->vnodes[k], -domain_azimut_angle, -domain_polar_angle);
	octree_get_vs(outer_img, outer_img_tree, &shift_r, BH_resolution, &v_shift);
	v_safe = shifted(&v_shift, domain_azimut_angle, domain_polar_angle);
	vec3_add(&v_shift_total, &v_shift_total, &v_safe);
	
    // Add the velocity due to boundary images to the result.
    vec3_add(&tangle->vs[k], &tangle->vs[k], &v_shift_total);

	// We have the superfluid velocity field at every node. In case of zero temperature, vortices are moving together with superfuid component.
	// In case of non zero temperatures, we have to find the action of normal component on motion of vortices.
    tangle->vels[k] = tangle->vs[k];

	// In non-zero temperature use also the mutual friction.
    if(use_mutual_friction)
    {
		// Define variables.
        struct vec3 tmp, dv;

        // The velocity difference. Here enters the velocity of the normal part.
        get_vn(&tangle->vnodes[k], t, &dv);
        vec3_sub(&dv, &dv, &tangle->vs[k]);

        // The dissipative term.
        vec3_cross(&tmp, &tangle->tangents[k], &dv);
        vec3_mul(&tmp, &tmp, alpha);
        vec3_add(&tangle->vels[k], &tangle->vels[k], &tmp);

        // The non-dissipative term.
        vec3_cross(&tmp, &tangle->tangents[k], &dv);
        vec3_cross(&tmp, &tangle->tangents[k], &tmp);
        vec3_mul(&tmp, &tmp, -alpha_p);
        vec3_add(&tangle->vels[k], &tangle->vels[k], &tmp);
    }

	// In the case of PINNED_SLIP subtract the normal component.
    if(tangle->status[k].status == PINNED_SLIP)
    {
        struct vec3 n = surface_normal(tangle->vnodes + k, tangle->status[k].pin_surface);
        double normal_velocity = vec3_dot(&n, &tangle->vels[k]);
        vec3_mul(&n, &n, normal_velocity);
        vec3_sub(&tangle->vels[k], &tangle->vels[k], &n);
    }
}

void update_one_velocity(struct tangle_state* tangle, int k, double t)
{
	// Create image tangles.
	struct tangle_state* inner_img = (struct tangle_state*)malloc(sizeof(struct tangle_state));
	struct tangle_state* outer_img = (struct tangle_state*)malloc(sizeof(struct tangle_state));
	create_tangle(inner_img, tangle->N);
	create_tangle(outer_img, tangle->N);

	// Only copy the stuff we actualy need. Here we need to reclculate domains.
	double domain_azimut_angle = tangle->domain_section.azimut_angle;
	double domain_polar_angle = tangle->domain_section.polar_angle;
	double domain_inner_radius = tangle->domain_section.inner_radius;
	double domain_outer_radius = tangle->domain_section.outer_radius;

	inner_img->domain_section.azimut_angle = domain_azimut_angle;
	inner_img->domain_section.polar_angle = domain_polar_angle;
	inner_img->domain_section.outer_radius = domain_inner_radius;
	inner_img->domain_section.inner_radius = domain_inner_radius * domain_inner_radius / domain_outer_radius;

	outer_img->domain_section.azimut_angle = domain_azimut_angle;
	outer_img->domain_section.polar_angle = domain_polar_angle;
	outer_img->domain_section.outer_radius = domain_outer_radius * domain_outer_radius / domain_inner_radius;
	outer_img->domain_section.inner_radius = domain_outer_radius;

	for (int kk = 0; kk < tangle->N; ++kk)
	{
		inner_img->status[kk] = tangle->status[kk];
		outer_img->status[kk] = tangle->status[kk];
		inner_img->connections[kk] = tangle->connections[kk];
		outer_img->connections[kk] = tangle->connections[kk];

		if (tangle->status[kk].status == EMPTY) continue;

		inner_img->vnodes[kk] = spherical_inversion(tangle->vnodes + kk, domain_inner_radius);
		outer_img->vnodes[kk] = spherical_inversion(tangle->vnodes + kk, domain_outer_radius);
	}

	// Build BH trees.
	struct octree* tree, * inner_img_tree, * outer_img_tree = NULL;
	tree = octree_build(tangle);
	inner_img_tree = octree_build(inner_img);
	outer_img_tree = octree_build(outer_img);

	// Velocity of node.
	update_velocity(tangle, inner_img, outer_img, k, t, tree, inner_img_tree, outer_img_tree);

	// Free BH trees.
	octree_destroy(tree);
	octree_destroy(inner_img_tree);
	octree_destroy(outer_img_tree);

	free_tangle(inner_img);
	free_tangle(outer_img);

	free(inner_img);
	free(outer_img);
}

/*
	Updates the velocities of all nodes of the tangle.
	@param tangle: Tangle structure holding all informations about vortices we want to update.
	@param t: Time of the simulation when we are updating the tangle.
*/
void update_velocities(struct tangle_state *tangle, double t)
{
	// Create image tangles.
	struct tangle_state* inner_img = (struct tangle_state*)malloc(sizeof(struct tangle_state));
	struct tangle_state* outer_img = (struct tangle_state*)malloc(sizeof(struct tangle_state));
	create_tangle(inner_img, tangle->N);
	create_tangle(outer_img, tangle->N);

	// Only copy the stuff we actualy need. Here we need to reclculate domains.
	double domain_azimut_angle = tangle->domain_section.azimut_angle;
	double domain_polar_angle = tangle->domain_section.polar_angle;
	double domain_inner_radius = tangle->domain_section.inner_radius;
	double domain_outer_radius = tangle->domain_section.outer_radius;


	inner_img->domain_section.azimut_angle = domain_azimut_angle;
	inner_img->domain_section.polar_angle = domain_polar_angle;
	inner_img->domain_section.outer_radius = domain_inner_radius;
	inner_img->domain_section.inner_radius = domain_inner_radius * domain_inner_radius / domain_outer_radius;

	outer_img->domain_section.azimut_angle = domain_azimut_angle;
	outer_img->domain_section.polar_angle = domain_polar_angle;
	outer_img->domain_section.outer_radius = domain_outer_radius * domain_outer_radius / domain_inner_radius;
	outer_img->domain_section.inner_radius = domain_outer_radius;

	for (int kk = 0; kk < tangle->N; ++kk)
	{
		inner_img->status[kk] = tangle->status[kk];
		outer_img->status[kk] = tangle->status[kk];
		inner_img->connections[kk] = tangle->connections[kk];
		outer_img->connections[kk] = tangle->connections[kk];

		if (tangle->status[kk].status == EMPTY) continue;

		inner_img->vnodes[kk] = spherical_inversion(tangle->vnodes + kk, domain_inner_radius);
		outer_img->vnodes[kk] = spherical_inversion(tangle->vnodes + kk, domain_outer_radius);
	}

	// Build BH trees.
	struct octree* tree, * inner_img_tree, * outer_img_tree = NULL;

	tree = octree_build(tangle);
	inner_img_tree = octree_build(inner_img);
	outer_img_tree = octree_build(outer_img);


	// Velocity of every node is updated separately, parallelize it.
    int i;
	//#pragma omp parallel private(i) num_threads(global_num_threads)
    {
		//#pragma omp for
        for(i = 0; i < tangle->N; ++i) update_velocity(tangle, inner_img, outer_img, i, t, tree, inner_img_tree, outer_img_tree);
    }

	// Free BH trees.
    octree_destroy(tree);
	octree_destroy(inner_img_tree);
	octree_destroy(outer_img_tree);

	free_tangle(inner_img);
	free_tangle(outer_img);

	free(inner_img);
	free(outer_img);
}

/*******************************************************************************
 *********** TANGLE FUNCTIONS FOR BOUNDARIES, MALL LOOPS REMESH ****************
 ******************************************************************************/

/*
	Enforce periodc bounadries, move all points inside the domain section.
	@param tangle: Tangle structure that holds all information about the nodes and domain box.
*/
void enforce_periodic_boundaries(struct tangle_state *tangle)
{
	// Loop thourgh all vortices.
	for(int k=0; k < tangle->N; ++k)
    {
		// Continue if empty.
		if(tangle->status[k].status == EMPTY) continue;

		// Angles of the domain.
		double domain_azimut_angle = tangle->domain_section.azimut_angle;
		double domain_polar_angle = tangle->domain_section.polar_angle;

		// Check azimutal angles.
		double azimut_node = azimut_angle(tangle->vnodes + k);
		double azimut_tangent = azimut_angle(tangle->tangents + k);
		double azimut_normal = azimut_angle(tangle->normals + k);
		double azimut_vels = azimut_angle(tangle->vels + k);
		double azimut_vs = azimut_angle(tangle->vs + k);

		while (azimut_node > domain_azimut_angle)
		{
			azimut_node -= 2 * domain_azimut_angle;
			azimut_tangent -= 2 * domain_azimut_angle;
			azimut_normal -= 2 * domain_azimut_angle;
			azimut_vels -= 2 * domain_azimut_angle;
			azimut_vs -= 2 * domain_azimut_angle;
		}
		while (azimut_node < -domain_azimut_angle)
		{
			azimut_node += 2 * domain_azimut_angle;
			azimut_tangent += 2 * domain_azimut_angle;
			azimut_normal += 2 * domain_azimut_angle;
			azimut_vels += 2 * domain_azimut_angle;
			azimut_vs += 2 * domain_azimut_angle;
		}

		// Check polar angles.
		double polar_node = polar_angle(tangle->vnodes + k);
		double polar_tangent = polar_angle(tangle->tangents + k);
		double polar_normal = polar_angle(tangle->normals + k);
		double polar_vels = polar_angle(tangle->vels + k);
		double polar_vs = polar_angle(tangle->vs + k);

		while (polar_node > domain_polar_angle)
		{
			polar_node -= 2 * domain_polar_angle;
			polar_tangent -= 2 * domain_polar_angle;
			polar_normal -= 2 * domain_polar_angle;
			polar_vels -= 2 * domain_polar_angle;
			polar_vs -= 2 * domain_polar_angle;
		}
		while (polar_node < -domain_polar_angle)
		{
			polar_node += 2 * domain_polar_angle;
			polar_tangent += 2 * domain_polar_angle;
			polar_normal += 2 * domain_polar_angle;
			polar_vels += 2 * domain_polar_angle;
			polar_vs += 2 * domain_polar_angle;
		}

		// Reverse coordinates.
		double r_node = radius(tangle->vnodes + k);
		double r_tangent = radius(tangle->tangents + k);
		double r_normal = radius(tangle->normals + k);
		double r_vels = radius(tangle->vels + k);
		double r_vs = radius(tangle->vs + k);

		tangle->vnodes[k] = spherical_to_vector(r_node, azimut_node, polar_node);
		tangle->tangents[k] = spherical_to_vector(r_tangent, azimut_tangent, polar_tangent);
		tangle->normals[k] = spherical_to_vector(r_normal, azimut_normal, polar_normal);
		tangle->vels[k] = spherical_to_vector(r_vels, polar_vels, polar_vels);
		tangle->vs[k] = spherical_to_vector(r_vs, polar_vs, polar_vs);
    }
}

// Functions needed for enfrce wall boundaries.

void remesh_segment(struct tangle_state* tangle, int pt) 
{
	int next = tangle->connections[pt].forward;

	struct segment seg = seg_pwrap(&tangle->vnodes[pt], &tangle->vnodes[next], &tangle->domain_section);

	//do we need an extra point?
	if (segment_len(&seg) > global_dl_max)
	{
		int new_pt = add_point(tangle, pt);

		// Update tangets and normals where needed.
		update_tangent_normal(tangle, new_pt);
		update_tangent_normal(tangle, pt);
		update_tangent_normal(tangle, next);

		if (tangle->connections[pt].reverse > -1) update_tangent_normal(tangle, tangle->connections[pt].reverse);
		if (tangle->connections[next].forward > -1) update_tangent_normal(tangle, tangle->connections[next].forward);

		// Do we need more extra points?
		remesh_segment(tangle, pt);
		remesh_segment(tangle, new_pt);
	}
}

void move_pinned_points_onto_the_surface(struct tangle_state* tangle)
{
	// Looks if any PINNED_SLIP point havent moved away from wall. Yes it is hard to make radial motion on the wall other way.
	for (int k = 0; k < tangle->N; ++k)
	{
		if (tangle->status[k].status == PINNED_SLIP || tangle->status[k].status == PINNED)
		{
			vec3_normalize(tangle->vnodes + k);

			if (tangle->status[k].pin_surface == INNER)
			{
				vec3_mul(tangle->vnodes + k, tangle->vnodes + k, tangle->domain_section.inner_radius);
				if (tangle->connections[k].forward == -1) remesh_segment(tangle, tangle->connections[k].reverse);
				if (tangle->connections[k].reverse == -1) remesh_segment(tangle, k);
			}
			if (tangle->status[k].pin_surface == OUTER)
			{
				vec3_mul(tangle->vnodes + k, tangle->vnodes + k, tangle->domain_section.outer_radius);
				if (tangle->connections[k].forward == -1) remesh_segment(tangle, tangle->connections[k].reverse);
				if (tangle->connections[k].reverse == -1) remesh_segment(tangle, k);
			}
		}
	}
}

//
void remove_and_pin_neighbours(struct tangle_state* tangle, int point_idx, int wall)
{
	// Finds id of previous and next nodes.
	int prev = tangle->connections[point_idx].reverse;
	int next = tangle->connections[point_idx].forward;

	if (prev >= 0)
	{
		tangle->connections[prev].forward = -1;
		tangle->status[prev].pin_surface = wall;
		tangle->status[prev].status = pin_mode;
	}
	if (next >= 0)
	{
		tangle->connections[next].reverse = -1;
		tangle->status[next].pin_surface = wall;
		tangle->status[next].status = pin_mode;
	}

	// Set the node EMPTY.
	tangle->connections[point_idx].reverse = -1;
	tangle->connections[point_idx].forward = -1;
	tangle->status[point_idx].status = EMPTY;
	tangle->status[point_idx].pin_surface = NOT_A_SURFACE;
}

void remove_all_behind_walls(struct tangle_state* tangle)
{
	for (int k = 0; k < tangle->N; ++k)
	{
		if (tangle->status[k].status == EMPTY ||
			tangle->status[k].status == PINNED ||
			tangle->status[k].status == PINNED_SLIP)
			continue;

		double r = radius(tangle->vnodes + k);

		if (r < tangle->domain_section.inner_radius) remove_and_pin_neighbours(tangle, k, INNER);
		if (r > tangle->domain_section.outer_radius) remove_and_pin_neighbours(tangle, k, OUTER);
	}
}

//
void remove_with_neighbours_in_rec_dist(struct tangle_state* tangle, int k, int wall)
{
	int next = tangle->connections[k].forward;
	double next_r = 0;
	if (next > -1) next_r = radius(tangle->vnodes + next);
	while (next > -1 && next_r < tangle->domain_section.inner_radius + rec_dist_with_walls)
	{
		int current = next;
		next = tangle->connections[k].forward;
		next_r = 0;
		if (next > -1) next_r = radius(tangle->vnodes + next);
		remove_and_pin_neighbours(tangle, current, wall);
	}

	int prev = tangle->connections[k].reverse;
	double prev_r = 0;
	if (prev > -1) prev_r = radius(tangle->vnodes + prev);
	while (prev > -1 && prev_r > tangle->domain_section.inner_radius - rec_dist_with_walls)
	{
		int current = prev;
		prev = tangle->connections[k].reverse;
		prev_r = 0;
		if (prev > -1) prev_r = radius(tangle->vnodes + prev);
		remove_and_pin_neighbours(tangle, current, wall);
	}

	remove_and_pin_neighbours(tangle, k, wall);
}

void pin_vortices_in_rec_distance(struct tangle_state* tangle)
{
	for (int k = 0; k < tangle->N; ++k)
	{
		if (tangle->status[k].status == EMPTY ||
			tangle->status[k].status == PINNED ||
			tangle->status[k].status == PINNED_SLIP)
			continue;

		double r = radius(tangle->vnodes + k);

		int is_parallel = 0;
		int is_going_away = 0;

		// Do not remove parallel vortices if they are near the wall.

		if (r < tangle->domain_section.inner_radius + rec_dist_with_walls)
		{
			struct vec3 normal = surface_normal(tangle->vnodes + k, INNER);
			if (abs(vec3_ndot(&tangle->tangents[k], &normal)) < sin(reconnection_angle_cutoff)) is_parallel = 1;
			if (vec3_ndot(&tangle->vels[k], &normal) > 0) is_going_away = 1;

			if (is_going_away && is_parallel) remove_with_neighbours_in_rec_dist(tangle, k, INNER);
		}
		if (r > tangle->domain_section.outer_radius - rec_dist_with_walls)
		{
			struct vec3 normal = surface_normal(tangle->vnodes + k, OUTER);
			if (abs(vec3_ndot(&tangle->tangents[k], &normal)) < sin(reconnection_angle_cutoff)) is_parallel = 1;
			if (vec3_ndot(&tangle->vels[k], &normal) > 0) is_going_away = 1;

			if (is_going_away && is_parallel)remove_with_neighbours_in_rec_dist(tangle, k, OUTER);
		}
	}
}

/*
	Enforce wall bounadries, move all points inside the domain box and reconnect points close to a solid wall.
	@param tangle: Tangle structure that holds all information about the nodes and domain box.
	@param t: Time of the simulation, to deretmine the velocity of the point we are clipping to the wall.
*/
void enforce_wall_boundaries(struct tangle_state* tangle)
{	
	move_pinned_points_onto_the_surface(tangle);
	
	remove_all_behind_walls(tangle);
 
	pin_vortices_in_rec_distance(tangle);

	move_pinned_points_onto_the_surface(tangle);
}


// Does it take into account that it may be needed to add more then one point???
void remesh(struct tangle_state *tangle, double min_dist, double max_dist)
{
	int added = 0;
	for(int k=0; k<tangle->N; ++k)
    {
		if(tangle->status[k].status == EMPTY) continue;

		int next = tangle->connections[k].forward;
		int prev = tangle->connections[k].reverse;

		struct segment sf;
		struct segment sr;
		double lf = 0; // To suppres warnings.
		double lr = 0; // To suppres warnings.

		if(next >= 0)
		{
			sf = seg_pwrap(&tangle->vnodes[k], &tangle->vnodes[next], &tangle->domain_section);
			lf = segment_len(&sf);
		}
		if(prev >= 0)
		{
			sr = seg_pwrap(&tangle->vnodes[k], &tangle->vnodes[prev], &tangle->domain_section);
			lr = segment_len(&sr);
		}

		//can we remove point k?
		if(next >= 0 && prev >= 0 && ((lf < min_dist || lr < min_dist) && (lf + lr) < max_dist )) remove_point(tangle, k);

		// Do we need an extra point?
		if(next >= 0 && (lf > max_dist)) //since we are adding between k and next, check only lf
		{
			added++;
			int new_pt = add_point(tangle, k);
			// We may need to update all 4 points around it.
			update_tangent_normal(tangle, prev);
			update_tangent_normal(tangle, k);
			update_tangent_normal(tangle, new_pt);
			update_tangent_normal(tangle, next);
			if (tangle->connections[next].forward > -1) update_tangent_normal(tangle, tangle->connections[next].forward);
		}
		// Do we need an extra point in the segment before? This solves the probem of adding more then one point.
		if (prev >= 0 && (lr > max_dist)) //since we are adding between k and next, check only lf
		{
			added++;
			int new_pt = add_point(tangle, prev);
			// We may need to update all 4 points around it.
			if (tangle->connections[prev].reverse > -1) update_tangent_normal(tangle, tangle->connections[prev].reverse);
			update_tangent_normal(tangle, prev);
			update_tangent_normal(tangle, new_pt);
			update_tangent_normal(tangle, k);
			update_tangent_normal(tangle, next);
		}
	}
	//we could have added points outside of the domain
	if(added){
	    enforce_periodic_boundaries(tangle);
		enforce_wall_boundaries(tangle);
	}
}

void eliminate_small_loops(struct tangle_state *tangle)
{
	// Prepare array to remember visited points.
	int visited[tangle->N];
	for (int i = 0; i < tangle->N; i++) {
		visited[i] = 0;
	}
	
	// Loop through the vortices.
    for(int k=0; k < tangle->N; ++k)
    {
		// Skip empty or already visited points.
        if(tangle->status[k].status == EMPTY || visited[k]) continue;
		
		// Mark point as visited;
		visited[k] = 1;

		// Length of the vortex.
        int loop = 0;

		// Loop over the vortex.
        int here = k;
        int next = tangle->connections[here].forward;
        while(next != k)
        {
			// The vortex is not closed.
            if((tangle->status[here].status == PINNED || tangle->status[here].status == PINNED_SLIP) && next < 0)
            {
                // We hit a wall, turn back from k.
                here = k;
                next = tangle->connections[here].reverse;
                while(next >= 0)
                {
					visited[here] = 1;
                    here = next;
                    next = tangle->connections[here].reverse;
                    loop++;
                }
                // 'here' now points to a node with its back to the wall.
                break; // Exit the outer loop.
            }

			// The vortex is closed.
			visited[here] = 1;
            here = next;
            next = tangle->connections[here].forward;
            loop++;
        }
		// The loop is short, delete it.
        if(loop < small_loop_cutoff)
        {
            // For loops, the starting point doesn't matter, but for wall-pinned lines the code bellow only goes
            // forward, so we have to start at the end facing away from the wall.
            next = here;
            while(1)
            {
                int tmp = next;
                next = tangle->connections[next].forward;
                tangle->connections[tmp].forward = -1;
                tangle->connections[tmp].reverse = -1;
                tangle->status[tmp].status = EMPTY;
                if(next == here || next < 0) break;
            }
        }
    }
}

