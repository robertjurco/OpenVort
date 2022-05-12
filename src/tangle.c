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
		tangle->status[k].pin_wall = -1;
		tangle->connections[k].forward = -1;
		tangle->connections[k].reverse = -1;
    }

	// Default initialisation is to open bounadry conditions.
	struct domain_box box = {
		.bottom_left_back = {{0,0,0}},
		.top_right_front  = {{1,1,1}},
		.wall = {WALL_OPEN, WALL_OPEN, WALL_OPEN, WALL_OPEN, WALL_OPEN, WALL_OPEN}
	};

	tangle->box = box;

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
            tangle->status[k].pin_wall = -1;
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
	if(where == 0) return tangle->vnodes[i];

	// If node is free, recursively calls itself.
	if(tangle->status[i].status == FREE)
	{
		if(where > 0) return step_node(tangle, tangle->connections[i].forward, where-1);
		else if (where < 0)	return step_node(tangle, tangle->connections[i].reverse, where+1);
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
		// If we are here it means we have ran into a wall and we need to flip the node.
		return mirror_shift(&out, &tangle->box, tangle->status[i].pin_wall);
	}
	else
	{
		printf("Walking across empty node."); // We should never get here.
		return vec3(0,0,0);
	}

	// To suppress warning.
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
	while(idx < 0 || tangle->status[idx].status != EMPTY)
	{
		expand_tangle(tangle, 2*tangle->N);
		idx = search_next_empty_node(tangle);
	}

	return idx;
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
	tangle->status[new_pt].pin_wall = NOT_A_FACE;

	// Find next (p+1) point.
	int next = tangle->connections[p].forward;

	// Update tangents and normals of points p and p+1.
	update_tangent_normal(tangle, p);
	update_tangent_normal(tangle, next);
	struct vec3 s0pp = tangle->normals[p];
	struct vec3 s1pp = tangle->normals[next];

	// Find positions of points p and p+1, pwraps it.
	struct vec3 s0 = tangle->vnodes[p];
	struct vec3 s1 = tangle->vnodes[next];
	struct segment seg = seg_pwrap(&s0, &s1, &tangle->box);
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
	tangle->vs[new_pt] = vec3(0.0, 0.0, 0.0);
	tangle->vels[new_pt] = vec3(0.0, 0.0, 0.0);
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
	Helper function that shifts the position vector r according the information containd in image_tangle,
	i.e. positive/negative number of shifts in the units of tangle box size with posibillity of mirror reflect.
	Mirror images of only depth 1 are supported (higher depth could make sense for two mirror facing each other).
	@param shift: Positive/negative number of shifts in the units of tangle box size with posibillity of mirror reflect.
	@param tangle: Structure holding all the informations about vortices and the size of domain box.
	@param r: Position vector we want to shift.
	@return Returns shifted position vector.
 */
struct vec3 shifted(const struct image_tangle* shift, const struct tangle_state* tangle, const struct vec3* r)
{
	// Domain box size.
	struct vec3 rs = *r;
	double Ls[3];
	for (int k = 0; k < 3; ++k) Ls[k] = tangle->box.top_right_front.p[k] - tangle->box.bottom_left_back.p[k];

	// First apply eriodic shift.
	for (int k = 0; k < 3; ++k) rs.p[k] -= Ls[k] * shift->shift[k];

	// Check for mirror reflect - solid wall.
	if (shift->number_of_reflects > 0)
	{
		for (int k = 0; k < shift->number_of_reflects; k++)
		{
			switch (shift->reflect[k])
			{
			case LEFT:
				rs.p[0] = 2 * tangle->box.bottom_left_back.p[0] - r->p[0];
				break;
			case RIGHT:
				rs.p[0] = 2 * tangle->box.top_right_front.p[0] - r->p[0];
				break;
			case BACK:
				rs.p[1] = 2 * tangle->box.bottom_left_back.p[1] - r->p[1];
				break;
			case FRONT:
				rs.p[1] = 2 * tangle->box.top_right_front.p[1] - r->p[1];
				break;
			case DOWN:
				rs.p[2] = 2 * tangle->box.bottom_left_back.p[2] - r->p[2];
				break;
			case UP:
				rs.p[2] = 2 * tangle->box.top_right_front.p[2] - r->p[2];
				break;
				// No reflection.
			default:
				break;
			}
		}
	}

	// Return shifter vector.
	return rs;
}


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

	dseg[0] = seg_pwrap(&s0, &s2, &tangle->box);
	dseg[1] = seg_pwrap(&s0, &s1, &tangle->box);
	dseg[2] = seg_pwrap(&s0, &sm1, &tangle->box);
	dseg[3] = seg_pwrap(&s0, &sm2, &tangle->box);
	dseg_12 = seg_pwrap(&s1, &s2, &tangle->box);
	dseg_m12 = seg_pwrap(&sm1, &sm2, &tangle->box);


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

	// WHY IS THIS HERE ???
	//double x = vec3_d(&tangle->tangents[k]);
	//vec3_mul(&tangle->normals[k], &tangle->normals[k], 1/x/x);
	//vec3_normalize(&tangle->tangents[k]);
}

/*
	Updates the tangents and normals of all nodes of the tangle. Calculation for k-th node is done from 5 nodes (k-2, k-1, k, k+1, k+2).
	ADD THE NAME OF THE METHOD (I can't remember it now).
	@param tangle: Tangle structure holding all informations about vortices we want to update.
*/
void update_tangents_normals(struct tangle_state* tangle)
{
	int i;
	for (i = 0; i < tangle->N; ++i)
	{
		update_tangent_normal(tangle, i);
	}
}

/*
	Calculates the velocity field of a straight vortex segment given by nodes i and i+1 at position r.
	@param tangle: Tangle structure holding all informations about vortices.
	@param i: ID of node that together with next one makes segment which velocity field we want calculate.
	@param r: Position where we want the field to evaluate.
	@returns Returns the velocity field of a straight vortex segment given by nodes i and i+1 at position r.
 */
static inline struct vec3 segment_field(const struct tangle_state *tangle, size_t i, struct vec3 r)
{
	// Look up the next node.
	int next = tangle->connections[i].forward;
	
	// This is an edge point on a wall, no segment.
	if(next == -1) return vec3(0,0,0);

	// Create the segment, pwrap it.
	struct segment seg = seg_pwrap(tangle->vnodes + i, tangle->vnodes + next, &tangle->box);

	// Calculate the vectors from r to the endpoints of the segment.
	struct vec3 R, Rp1;
	vec3_sub(&R, &seg.r1, &r);
	vec3_sub(&Rp1, &seg.r2, &r);
	double lR = vec3_len(&R);
	double lRp1 = vec3_len(&Rp1);

	// Calculate the coefficients for velocity field.
	double denom = lR * lRp1 * (lR * lRp1 + vec3_dot(&R, &Rp1));
	double f = KAPPA / 4 / M_PI;

	// This can happen in periodic boundary conditions.
	// TODO: the logic should be moved higher, test if this happen sometimes?
	if (lR < 1e-8 || lRp1 < 1e-8) return vec3(0, 0, 0);

	// If R and Rp1 are colinear, the result is 0, but code below would try to calculate 0/0.
	if (vec3_ndot(&R, &Rp1) < 1e-8) return vec3(0, 0, 0);

	// Calculate the velocity field.
	struct vec3 vv;
	vec3_cross(&vv, &R, &Rp1);
	vec3_mul(&vv, &vv, f * (lR + lRp1) / denom);

	return vv;
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
	struct segment sf = seg_pwrap(p, &next, &tangle->box);
	struct segment sr = seg_pwrap(&prev, p, &tangle->box);

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

/*
	Updates the velocities of k-th node of the tangle.
	@param tangle: Tangle structure holding all informations about vortices we want to update.
	@param t: Time of the simulation when we are updating the tangle.
	@param tree: Pointer to the BH tree. If we are not using BH algortihm *tree = NULL;
*/
void update_velocity(struct tangle_state *tangle, int k, double t, struct octree *tree)
{
    int m;

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
    struct vec3 shift_r, v_shift;
    struct vec3 v_shift_total = VEC_NULL;

	// Case of using BH algorithm.
    if(use_BH && tree)
    {
        // Use the tree approximation to the full Biot-Savart.
        struct vec3 v_tree;
        octree_get_vs(tangle, tree, &tangle->vnodes[k], BH_resolution, &v_tree);
        vec3_add(&tangle->vs[k], &tangle->vs[k], &v_tree);

        // Calculate the velocity due to boundary images.
        for(int j = 0; j < tangle->bimg.n; ++j)
        {
			// Velocity from the boundary image.
            shift_r = shifted(&tangle->bimg.images[j], tangle, &tangle->vnodes[k]);
			octree_get_vs(tangle, tree, &shift_r, BH_resolution, &v_shift);
			// The case of the mirror wall.
			if (tangle->bimg.images[j].number_of_reflects > 0)
				for (int k = 0; k < tangle->bimg.images[j].number_of_reflects; k++)
					v_shift = mirror_vector_reflect(&v_shift, tangle->bimg.images[j].reflect[k]);
			// Add the result.
            vec3_add(&v_shift_total, &v_shift_total, &v_shift);
        }
    }
	// Case of NOT using BH algorithm.
    else
    {
        // Integrate using Biot-Savart law as usual.
        for (m = 0; m < tangle->N; ++m)
        {
			// Continue if the point is on the wall, or the segment is contained in LIA. 
            if(tangle->connections[m].forward == -1 || k == m || k == tangle->connections[m].forward) continue;
			// Otherwise calculate inducet velocity.
            struct vec3 segment_vel = segment_field(tangle, m, tangle->vnodes[k]);
            vec3_add(tangle->vs + k, tangle->vs + k, &segment_vel);
        }

        // Calculate the velocity due to boundary images.
        for(int j = 0; j < tangle->bimg.n; ++j)
        {
			// Move position to right place with respect to boundary image.
            shift_r = shifted(&tangle->bimg.images[j], tangle, &tangle->vnodes[k]);
			v_shift = vec3(0, 0, 0);

			// Integrate velocity from the boundary image using Biot-Savart law.
			for (int m = 0; m < tangle->N; ++m)
			{
				// Continue if the point is on the wall.
				if (tangle->connections[m].forward == -1) continue;
				// Otherwise calculate inducet velocity.
				struct vec3 segment_vel = segment_field(tangle, m, shift_r);
				vec3_add(tangle->vs + k, tangle->vs + k, &segment_vel);
			}
			// The case of the mirror wall.
			if (tangle->bimg.images[j].number_of_reflects > 0)
				for (int k = 0; k < tangle->bimg.images[j].number_of_reflects; k++)
					v_shift = mirror_vector_reflect(&v_shift, tangle->bimg.images[j].reflect[k]);
			// Add the result.
            vec3_add(&v_shift_total, &v_shift_total, &v_shift);
        }
    }

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
        struct vec3 n = boundary_normals[tangle->status[k].pin_wall];
        double normal_velocity = vec3_dot(&n, &tangle->vels[k]);
        vec3_mul(&n, &n, normal_velocity);
        vec3_sub(&tangle->vels[k], &tangle->vels[k], &n);
    }
}

/*
	Updates the velocities of all nodes of the tangle.
	@param tangle: Tangle structure holding all informations about vortices we want to update.
	@param t: Time of the simulation when we are updating the tangle.
*/
void update_velocities(struct tangle_state *tangle, double t)
{
	// Build BH tree if needed.
    struct octree *tree = NULL;
    if (use_BH) tree = octree_build(tangle);

	// Velocity of every node is updated separately, parallelize it.
    int i;
    #pragma omp parallel private(i) num_threads(global_num_threads)
    {
        #pragma omp for
        for(i = 0; i < tangle->N; ++i) update_velocity(tangle, i, t, tree);
    }

	// Free BH tree.
    octree_destroy(tree);
}

/*******************************************************************************
 *********** TANGLE FUNCTIONS FOR BOUNDARIES, MALL LOOPS REMESH ****************
 ******************************************************************************/

/*
	Checks if the position is inside the domain box.
	@param tangle: Tangle structure that holds all information about the nodes and domain box.
	@param where: The vector giving us the position we want to check.
	@returns Returns -1 if point is inside the domain box. Otherwise returns the face behind which the point is.
*/
static inline int out_of_box(const struct tangle_state *tangle, const struct vec3 *where)
{
	double x = where->p[0];
	double y = where->p[1];
	double z = where->p[2];
	double bounds[] =
	{
		tangle->box.bottom_left_back.p[0],
		tangle->box.top_right_front.p[0],
		tangle->box.bottom_left_back.p[1],
		tangle->box.top_right_front.p[1],
		tangle->box.bottom_left_back.p[2],
		tangle->box.top_right_front.p[2]
     };

	int face = -1;
	if(x < bounds[LEFT]) face = LEFT;
	if(x > bounds[RIGHT]) face = RIGHT;
	if(y < bounds[BACK]) face = BACK;
	if(y > bounds[FRONT]) face = FRONT;
	if(z < bounds[DOWN]) face = DOWN;
	if(z > bounds[UP]) face = UP;

	if(face > -1 && tangle->box.wall[face] != WALL_PERIODIC) face = -1;

	return face;
}

/*
	Enforce periodc bounadries, move all points inside the domain box.
	@param tangle: Tangle structure that holds all information about the nodes and domain box.
*/
void enforce_periodic_boundaries(struct tangle_state *tangle)
{
	// Loop thourgh all vortices.
	for(int k=0; k < tangle->N; ++k)
    {
		// Continue if empty.
		if(tangle->status[k].status == EMPTY) continue;


		int face = out_of_box(tangle, &tangle->vnodes[k]);
		if(face >= 0)
		{
			// This should be only possible with periodic faces.
			if(tangle->status[k].status == PINNED || tangle->status[k].status == PINNED_SLIP)
				printf("Pinned node outside of the box, this should have been caught with enforce_wall_boundaries.");

			// Shift nodes inside the box.
			while((face=out_of_box(tangle, &tangle->vnodes[k])) >= 0)
				tangle->vnodes[k] = periodic_shift(&tangle->vnodes[k], &tangle->box, face);
		}
    }
}

/*
	Checks whether a node k is closer than rec_dist to the wall.
	@param tangle: Tangle structure that holds all information about the nodes and domain box.
	@param k: ID of the node we want to check.
	@param wall: Which wall are we checking.
	@return Return 1 if the node is closer than rec_dist to the wall, otherwise return 0.
*/
int check_wall(struct tangle_state* tangle, int k, int wall)
{
	int idx[6];
	idx[LEFT] = idx[RIGHT] = 0;
	idx[BACK] = idx[FRONT] = 1;
	idx[DOWN] = idx[UP]	   = 2;
	switch (wall)
	{
		case LEFT: case BACK: case DOWN:
			return tangle->vnodes[k].p[idx[wall]] < (rec_dist + tangle->box.bottom_left_back.p[idx[wall]]);
	
		case RIGHT: case FRONT:	case UP:
			return tangle->vnodes[k].p[idx[wall]] > (tangle->box.top_right_front.p[idx[wall]] - rec_dist);
		default:
			break;
	}
	return 0;
}

/*
	Returns the distance of point k to the specified wall.
	@param tangle : Tangle structure that holds all information about the nodesand domain box.
	@param k : ID of the node we want to check.
	@param wall : Which wall are we checking.
	@return Return 1 if the node is closer than rec_dist to the wall, otherwise return 0.
*/
double wall_dist(const struct tangle_state* tangle, int k, boundary_faces wall)
{
	int idx[6];
	idx[LEFT] = idx[RIGHT] = 0;
	idx[BACK] = idx[FRONT] = 1;
	idx[DOWN] = idx[UP]    = 2;
	switch (wall)
	{
		case LEFT: case BACK: case DOWN:
			return tangle->vnodes[k].p[idx[wall]] - tangle->box.bottom_left_back.p[idx[wall]];

		case RIGHT: case FRONT: case UP:
			return tangle->box.top_right_front.p[idx[wall]] - tangle->vnodes[k].p[idx[wall]];

		default:
			printf("wall_dist: unknown wall index %d\n", wall);
	}
	return -1;
}

/*
	Reconnect point close to a solid wall.
	@param tangle : Tangle structure that holds all information about the nodes and domain box.
	@param k : ID of the node we want to check.
	@param wall : Which wall are we checking.
	@param t: Time of the simulation, to deretmine the velocity of the point.
	@return Return 1 if the node is succesfully clipped, otherwise return 0.
*/
int connect_to_wall(struct tangle_state* tangle, int k, int wall, double time)
{
	// Update point.
	update_tangent_normal(tangle, k);
	update_velocity(tangle, k, time, NULL);

	// Find its distance from the wall.
	double d0 = wall_dist(tangle, k, wall);

	// Check that the node is actually getting closer to the wall, if not do not clip it.
	// Boundary_normals are inward-facing unit normals.
	if (d0 > 0 && vec3_dot(&tangle->vels[k], &boundary_normals[wall]) > 0)
		return 0;

	// Check that the vortex is not parallel to the boundary.
	if (abs(vec3_ndot(&tangle->tangents[k], &boundary_normals[wall])) > cos(reconnection_angle_cutoff))
		return 0;

	// Find all the info about the neighbour nodes.
	int next = tangle->connections[k].forward;
	int prev = tangle->connections[k].reverse;
	double d1 = wall_dist(tangle, next, wall);
	double dm1 = wall_dist(tangle, prev, wall);

	struct segment s1 = seg_pwrap(&tangle->vnodes[k], &tangle->vnodes[next], &tangle->box);
	struct segment sm1 = seg_pwrap(&tangle->vnodes[k], &tangle->vnodes[prev], &tangle->box);
	double vd1 = segment_len(&s1);
	double vdm1 = segment_len(&sm1);

	// Check that the total length does not increase.
	if (vd1 < d0 + d1 && vdm1 < d0 + dm1) return 0;

	struct vec3 tmp;

	// New point 1.
	int new_pt = get_next_empty_node(tangle);
	tangle->status[new_pt].status = pin_mode;
	tangle->status[new_pt].pin_wall = wall;

	// New point 2.
	int new_pt2 = get_next_empty_node(tangle);
	tangle->status[new_pt2].status = pin_mode;
	tangle->status[new_pt2].pin_wall = wall;

	// Clip the new point 1.
	vec3_mul(&tmp, &boundary_normals[wall], -d0);
	vec3_add(&tangle->vnodes[new_pt], &tangle->vnodes[k], &tmp);

	// Flag the changed points to avoid overcrowding reconnections.
	tangle->recalculate[k]++;
	tangle->recalculate[new_pt]++;
	tangle->recalculate[new_pt2]++;
	tangle->recalculate[next]++;
	tangle->recalculate[prev]++;

	// Clip the new point 2 and finish clipping 1.
	if (d1 < dm1)
	{
		vec3_mul(&tmp, &boundary_normals[wall], -d1);
		vec3_add(&tangle->vnodes[new_pt2], &tangle->vnodes[next], &tmp);

		tangle->connections[k].forward = new_pt;
		tangle->connections[new_pt].forward = -1;
		tangle->connections[new_pt].reverse = k;
		tangle->connections[new_pt2].forward = next;
		tangle->connections[new_pt2].reverse = -1;
		tangle->connections[next].reverse = new_pt2;
	}
	else
	{
		vec3_mul(&tmp, &boundary_normals[wall], -dm1);
		vec3_add(&tangle->vnodes[new_pt2], &tangle->vnodes[prev], &tmp);

		tangle->connections[k].reverse = new_pt;
		tangle->connections[new_pt].forward = k;
		tangle->connections[new_pt].reverse = -1;
		tangle->connections[new_pt2].forward = -1;
		tangle->connections[new_pt2].reverse = prev;
		tangle->connections[prev].forward = new_pt2;
	}

	return 1;
}

/*
	Enforce wall bounadries, move all points inside the domain box and reconnect points close to a solid wall.
	@param tangle: Tangle structure that holds all information about the nodes and domain box.
	@param t: Time of the simulation, to deretmine the velocity of the point we are clipping to the wall.
*/
void enforce_wall_boundaries(struct tangle_state* tangle, double time)
{	
	for (int k = 0; k < tangle->N; ++k) tangle->recalculate[k] = 0;

	// Reconnect with walls.
	for (int wall = 0; wall < 6; ++wall)
	{
		// Skip periodic and open boundaries.
		if (tangle->box.wall[wall] == WALL_OPEN || tangle->box.wall[wall] == WALL_PERIODIC) continue;

		// Loop through all the points.
		for (int k = 0; k < tangle->N; ++k)
		{
			// Skip PINNED, PINNED_SLIP, EMPTY and already calculated nodes.
			if (tangle->status[k].status == PINNED ||
				tangle->status[k].status == PINNED_SLIP ||
				tangle->status[k].status == EMPTY || 
				tangle->recalculate[k])
				continue;

			// Clip the pointss close to the wall.
			if (check_wall(tangle, k, wall))
				connect_to_wall(tangle, k, wall, time);
		}
	}

	// Eliminate all points that are active and outside the walls. This can potentially happen
	// if the time step is too long and more than one discretisation point gets outside the domain.
	for (int wall = 0; wall < 6; ++wall)
	{
		// Skip periodic and open boundaries.
		if (tangle->box.wall[wall] == WALL_OPEN || tangle->box.wall[wall] == WALL_PERIODIC)
			continue;

		// Loop through all the points.
		for (int k = 0; k < tangle->N; ++k)
		{
			// Skip PINNED, PINNED_SLIP and EMPTY nodes.
			if (tangle->status[k].status == PINNED ||
				tangle->status[k].status == PINNED_SLIP ||
				tangle->status[k].status == EMPTY)
				continue;

			// Remove points behind the wall.
			if (wall_dist(tangle, k, wall) < 0)
			{
				remove_point(tangle, k);
			}
		}
	}
}

void remesh(struct tangle_state *tangle, double min_dist, double max_dist, double time)
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
			sf = seg_pwrap(&tangle->vnodes[k], &tangle->vnodes[next], &tangle->box);
			lf = segment_len(&sf);
		}
		if(prev >= 0)
		{
			sr = seg_pwrap(&tangle->vnodes[k], &tangle->vnodes[prev], &tangle->box);
			lr = segment_len(&sr);
		}

		//can we remove point k?
		if(next >= 0 && prev >= 0 && ((lf < min_dist || lr < min_dist) && (lf + lr) < max_dist )) remove_point(tangle, k);

		//do we need an extra point?
		if(next >= 0 && (lf > max_dist)) //since we are adding between k and next, check only lf
		{
			added++;
			int new_pt = add_point(tangle, k);
			update_tangent_normal(tangle, new_pt);
		}
	}
	//we could have added points outside of the domain
	if(added){
	    enforce_periodic_boundaries(tangle);
		enforce_wall_boundaries(tangle, time);
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

