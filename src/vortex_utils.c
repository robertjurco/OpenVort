
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vec3_math.h"
#include "vortex_utils.h"
#include "vortex_constants.h"
#include "tangle.h"

/*
 * Debugging utilities
 */

int check_loop(const struct tangle_state *tangle, int *visited, int k);
int check_integrity(const struct tangle_state *tangle)
{
  int *visited = calloc(tangle->N, sizeof(int));
  int errors = 0;

  for(int k=0; k < tangle->N; ++k)
    {
      int next = tangle->connections[k].forward;
      int prev = tangle->connections[k].reverse;
      if(next >= 0 && tangle->status[k].status == FREE)
	{
	  if(k != tangle->connections[next].reverse)
	    printf("Forward connection broken %d %d %d\n", k, next,
		  tangle->connections[next].reverse);
	}
      if(prev >= 0 && tangle->status[k].status == FREE)
	{
	  if( k!= tangle->connections[prev].forward)
	    printf("Reverse connection broken %d %d %d\n", k, prev,
		  tangle->connections[prev].forward);
	}
      if(!visited[k])
	{
	  visited[k] += 1;
	  errors += check_loop(tangle, visited, k);
	}
    }

  free(visited);
  return errors;
}

int is_empty(const struct tangle_state *tangle, int k)
{
  struct neighbour_t nb = tangle->connections[k];

  if(nb.forward == -1 &&
     nb.reverse == -1)
    return 1; //point is empty

  if(nb.forward >= 0 &&
     nb.reverse >= 0)
    return 0; //both links are valid, point is not empty

  //one of the links is >= 0, the other is ==-1
  //point is corrupted
  #ifdef _DEBUG_
  printf("Corrupted links in point %d (%d, %d)\n",
	 k,
	 tangle->connections[k].forward,
	 tangle->connections[k].reverse);
  #endif
  return -1;
}

int check_loop(const struct tangle_state *tangle, int *visited, int k)
{
  int errors = 0;
  int rval;

  if((rval = is_empty(tangle, k)) > 0)
    return 0;
  if(rval < 0)
    return rval;

  int j = tangle->connections[k].forward;
  int total = 0;
  while(j != k && total < tangle->N)
    {
      if(visited[j])
	{
#ifdef _DEBUG_
	  printf("Ran into point %d again.\n", j);
#endif
	  errors--;
	  break;
	}
      visited[j]++;
      if((rval = is_empty(tangle, j)) > 0)
	{
#ifdef _DEBUG_
	  printf("Connected to an empty point, %d %d.\n",
		 k, j);
#endif
	  errors--;
	  return errors;
	}
      if(tangle->connections[j].forward < 0)
	{
#ifdef _DEBUG_
	  printf("Linked point with -1 forward, %d.\n",
		 j);
#endif
	  errors--;
	  return errors;
	}

      j = tangle->connections[j].forward;

      total++;
    }


  if(j != k && total == tangle->N)
    {
#ifdef _DEBUG_
      printf("Ran through the entire tangle from %d.\n", k);
#endif
      errors--;
    }

  return errors;
}

/*
 * Wall-related stuff
 */

double wall_dist(const struct tangle_state *tangle, int k, boundary_faces wall)
{
  /*
   * Returns the distance of point k to the specified wall.
   */
  int idx[6];
  idx[LEFT] = idx[RIGHT] = 0;
  idx[BACK] = idx[FRONT] = 1;
  idx[DOWN] = idx[UP] = 2;
  switch(wall)
  {
    case LEFT:
    case BACK:
    case DOWN:
      return tangle->vnodes[k].p[idx[wall]] - tangle->box.bottom_left_back.p[idx[wall]];

    case RIGHT:
    case FRONT:
    case UP:
      return tangle->box.top_right_front.p[idx[wall]] - tangle->vnodes[k].p[idx[wall]];

    default:
      printf("wall_dist: unknown wall index %d\n", wall);
  }
  return -1;
}

void clip_at_wall(struct tangle_state *tangle)
{
  /*
   * Clips the tangle at the z-walls
   */
  //limits for unconstrained walls
  double llimit = -INFINITY;
  double ulimit = INFINITY;

  if(tangle->box.wall[DOWN] == WALL_MIRROR)
    llimit = tangle->box.bottom_left_back.p[2];
  if(tangle->box.wall[UP] == WALL_MIRROR)
    ulimit = tangle->box.top_right_front.p[2];

  enum POINT_STATES {
        OK=0,
        KILL,         //point to be clipped
        EDGE_FORWARD_L, //last point above the lower wall
        EDGE_REVERSE_L,  //first point above the lower wall
        EDGE_FORWARD_H, //last point below the upper wall
        EDGE_REVERSE_H  //first point below the upper wall
      };

  //recalculate is used as a tag for points to be clipped
  for(int k=0; k<tangle->N; ++k)
    tangle->recalculate[k] = OK;

  //first find all the points to be clipped
  for(int kk=0; kk<tangle->N; ++kk)
    {
      if(tangle->status[kk].status == EMPTY)
	continue;

      double zkk = tangle->vnodes[kk].p[2];
      if(zkk <= llimit || zkk >= ulimit)
	{
	  tangle->recalculate[kk] = KILL;
	  const int forward = tangle->connections[kk].forward;
	  const int reverse = tangle->connections[kk].reverse;

	  double zf = tangle->vnodes[forward].p[2];
	  double zr = tangle->vnodes[reverse].p[2];

	  if(zkk >= ulimit)
	    {
	      if(zf < ulimit)
		tangle->recalculate[forward] = EDGE_REVERSE_H;
	      if(zr < ulimit)
		tangle->recalculate[reverse] = EDGE_FORWARD_H;
	    }
	  if(zkk <= llimit)
	    {
	      if(zf > llimit)
		tangle->recalculate[forward] = EDGE_REVERSE_L;
	      if(zr > llimit)
		tangle->recalculate[reverse] = EDGE_FORWARD_L;
	    }
	}
    }

  //next handle edge points by projecting them on the wall
  //and pinning them on the lower or upper Z-wall
  for(int kk=0; kk < tangle->N; ++kk)
    {
      switch(tangle->recalculate[kk])
      {
	case EDGE_REVERSE_L:
	  tangle->connections[kk].reverse = -1;
	  tangle->vnodes[kk].p[2] = llimit;
	  tangle->status[kk].status = pin_mode;
	  tangle->status[kk].pin_wall = DOWN;
	  break;
	case EDGE_FORWARD_L:
	  tangle->connections[kk].forward = -1;
	  tangle->vnodes[kk].p[2] = llimit;
	  tangle->status[kk].status = pin_mode;
	  tangle->status[kk].pin_wall = DOWN;
	  break;
	case EDGE_REVERSE_H:
	  tangle->connections[kk].reverse = -1;
	  tangle->vnodes[kk].p[2] = ulimit;
	  tangle->status[kk].status = pin_mode;
	  tangle->status[kk].pin_wall = UP;
	  break;
	case EDGE_FORWARD_H:
	  tangle->connections[kk].forward = -1;
	  tangle->vnodes[kk].p[2] = ulimit;
	  tangle->status[kk].status = pin_mode;
	  tangle->status[kk].pin_wall = UP;
	  break;
	case KILL:
	  tangle->status[kk].status = EMPTY;
	  tangle->status[kk].pin_wall = NOT_A_FACE;
	  tangle->connections[kk].forward = -1;
	  tangle->connections[kk].reverse = -1;
	  break;
	default:
	  break;
      }
    }
}

double vortex_length(const struct tangle_state *tangle, int starting_point)
{
    double length = 0;

    struct vec3 res;
    int k = starting_point;
    int l;

    while (tangle->connections[k].reverse != -1 || tangle->connections[k].reverse != starting_point) {
        l = tangle->connections[k].reverse;
        vec3_sub(&res, &tangle->vnodes[k], &tangle->vnodes[l]);
        length = length + vec3_len(&res);
        k = l;
    }

    if (tangle->connections[k].reverse != starting_point) {
        vec3_sub(&res, &tangle->vnodes[k], &tangle->vnodes[starting_point]);
        length = length + vec3_len(&res);
    }

    return length;
}
