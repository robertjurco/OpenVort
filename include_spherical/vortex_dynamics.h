#ifndef VORTEX_DYNAMICS
#define VORTEX_DYNAMICS

#include "tangle.h"

/* 
	Structre holding information about four stepping tangles in RK4.
	TODO: The tangle states here should be reduced, only calculate what we really need for the velocity.
*/
struct rk4_state {
    struct tangle_state *s1, *s2, *s3, *s4;
    double dt;
};

/*
	Update the simulation using Runge-Kutta 4.
	@param result: The tangle to update.
	@param t: Current time of the simulation.
	@param dt: The time step.
*/
void rk4_step(struct tangle_state *tangle, double t, double dt);

/*
	Run through all the pairs of nodes and check their distance if they are not
	immediate neighbors. Reconnect them if they are close enough. 
	This could, and should, in the future include some smarter estimation of possibility
	of a reconnections -- i.e., where in the BSP tree the two nodes are and only check them
	if can, in principle, be close enough.
	@param tangle: Tangle structure holding all the information about the vorices.
	@param t: Current time of the simulation.
*/
int reconnect(struct tangle_state *tangle, double t);

#endif //VORTEX_DYNAMICS
