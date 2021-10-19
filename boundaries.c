#include <math.h>
#include <stdio.h>

#include "vec3_math.h"
#include "boundaries.h"


/*
	Inward-facing normals of the box boundary face walls.
 */
const struct vec3 boundary_normals[] = {
	{{1, 0, 0}}, // LEFT
    {{-1, 0, 0}}, // RIGHT
    {{0, 1, 0}},  // BACK
    {{0, -1, 0}}, // FRONT
    {{0, 0, 1}},  // DOWN
    {{0, 0, -1}}  // UP
};


/*
	Find out if the given vector vec is inside the box.
	@param box: Box where the vector should be located.
	@param vec: Vector we want to check if its inside the box.
	@returns Returns 1 if vec is inside the box, otherwise 0.
 */
int in_box(const struct domain_box *box, const struct vec3 *vec)
{
	double vx = vec->p[0];
	double vy = vec->p[1];
	double vz = vec->p[2];

	int in_box;
	in_box =
		vx > box->bottom_left_back.p[0] && vx < box->top_right_front.p[0] &&
		vy > box->bottom_left_back.p[1] && vy < box->top_right_front.p[1] &&
		vz > box->bottom_left_back.p[2] && vz < box->top_right_front.p[2];
	
	return in_box;
}

/*******************************************************************************
 ********************** PERIODIC AND MIRROR GEOMETRIES *************************
 ******************************************************************************/

 /*
	Periodically wraps the point r2 to follows after the point r1 inside the box.
	In the other words r2 is moved so, that it is in the same period of box as r1.
	@param r1: The reference point.
	@param r2: The next point we want to wrap periodically on the right place.
	@param box: Domain box determining the one period.
	@returns Returns segmen from r1 to r2 that is smaller then one perod of domain box.
 */
struct segment seg_pwrap(const struct vec3 *r1, const struct vec3 *r2, const struct domain_box *box)
{
    struct segment seg = {
		.r1 = *r1,
        .r2 = *r2
    };

	double Lx = box->top_right_front.p[0] - box->bottom_left_back.p[0];
    double Ly = box->top_right_front.p[1] - box->bottom_left_back.p[1];
    double Lz = box->top_right_front.p[2] - box->bottom_left_back.p[2];

	// Don't forget that it may be needed to wrap one point multiple times periodically.
    if(box->wall[LEFT] == WALL_PERIODIC)
    {
        while ((seg.r2.p[0] - seg.r1.p[0]) > Lx/2) seg.r2.p[0] -= Lx;
        while ((seg.r2.p[0] - seg.r1.p[0]) < -Lx/2) seg.r2.p[0] += Lx;
    }
    if(box->wall[BACK] == WALL_PERIODIC)
    {
        while ((seg.r2.p[1] - seg.r1.p[1]) > Ly/2) seg.r2.p[1] -= Ly;
        while ((seg.r2.p[1] - seg.r1.p[1]) < -Ly/2) seg.r2.p[1] += Ly;
    }
    if(box->wall[DOWN] == WALL_PERIODIC)
    {
        while ((seg.r2.p[2] - seg.r1.p[2]) > Lz/2) seg.r2.p[2] -= Lz;
        while ((seg.r2.p[2] - seg.r1.p[2]) < -Lz/2) seg.r2.p[2] += Lz;
    }

    return seg;
}

/*
   Periodically shifts the point by the given number of domain boxes in the given direction.
   @param v: The point we want to shift.
   @param box: Domain box determining the one period.
   @param shift[3]: Number of shifts in every directions in the size of domain box.
   @returns Returns shifted vector.
*/
struct vec3 box_shift(const struct vec3 *v, const struct domain_box *box, int shift[3])
{
	double Lx = box->top_right_front.p[0] - box->bottom_left_back.p[0];
	double Ly = box->top_right_front.p[1] - box->bottom_left_back.p[1];
	double Lz = box->top_right_front.p[2] - box->bottom_left_back.p[2];
	double Ls[] = {Lx, Ly, Lz};

	struct vec3 out = *v;

	for(int k = 0; k < 3; ++k) out.p[k] += Ls[k]*shift[k];

	return out;
}

/*
   Periodically shifts the point into the neighbour domain box given by the shared wall.
   @param v: The point we want to shift.
   @param box: Domain box determining the one period.
   @param wall: Wall though which we want periodically shift the point.
   @returns Returns shifted vector.
*/
struct vec3 periodic_shift(const struct vec3 *v, const struct domain_box *box, boundary_faces wall)
{
	struct vec3 mv = *v;

	int coord;

	double Lx = box->top_right_front.p[0] - box->bottom_left_back.p[0];
	double Ly = box->top_right_front.p[1] - box->bottom_left_back.p[1];
	double Lz = box->top_right_front.p[2] - box->bottom_left_back.p[2];
	double Ls[] = {Lx, Ly, Lz};

	switch(wall)
	{
		//  Which coordinate to shift.
	    case LEFT: case RIGHT: coord=0; break;
	    case BACK: case FRONT: coord=1; break;
	    case DOWN: case UP:    coord=2; break;
	    default:
			printf("Bad wall %d", wall);
			coord=-1;
			break;
	}	
	switch(wall)
	{
		// Shift in possitive or negative direction?
		case LEFT: case BACK: case DOWN:
			mv.p[coord] += Ls[coord];
			break;
		case RIGHT: case FRONT: case UP:
			mv.p[coord] -= Ls[coord];
			break;
		default:
			printf("Bad wall %d", wall);
			break;
	}

	return mv;
}

/*
   Reflects the point into the neighbour domain box given through the shared wall.
   @param v: The point we want to shift.
   @param box: Domain box determining the one period.
   @param wall: Wall though which we want mirror shift the point.
   @returns Returns shifted vector.
*/
struct vec3 mirror_shift(const struct vec3 *v, const struct domain_box *box, boundary_faces wall)
{
	struct vec3 mv = *v;

	int coord;
	double wall_pos;

	switch(wall)
	{
		case LEFT: case RIGHT: coord=0; break;
		case BACK: case FRONT: coord=1; break;
		case DOWN: case UP:    coord=2; break;
		default:
			printf("Bad wall %d", wall);
			coord=-1;
			break;
	}
	switch(wall)
	{
		case LEFT: wall_pos = box->bottom_left_back.p[0]; break;
		case RIGHT: wall_pos = box->top_right_front.p[0]; break;
		case BACK: wall_pos = box->bottom_left_back.p[1]; break;
		case FRONT: wall_pos = box->top_right_front.p[1]; break;
		case DOWN: wall_pos = box->bottom_left_back.p[2]; break;
		case UP: wall_pos = box->top_right_front.p[2]; break;
		default:
			printf("Bad wall %d", wall);
			wall_pos = -1;
			break;
	}

	mv.p[coord] = 2*wall_pos - mv.p[coord];

	return mv;
}

/*
	Reflect a directional vector in a mirror (i.e., the component normal to the mirror reverses).
	@param v: The vector e want to reflect.
	@param wall: The wall though which we want to reflect the vector.
	@returns Returns shifted vector.
*/
struct vec3 mirror_vector_reflect(const struct vec3 *v, boundary_faces wall)
{  
	struct vec3 mv = *v;

	double nc = vec3_dot(&mv, &boundary_normals[wall]);
	struct vec3 tmp;
	vec3_mul(&tmp, &boundary_normals[wall], -2*nc);
	vec3_add(&mv, &mv, &tmp);

	return mv;
}


