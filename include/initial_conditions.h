#ifndef INCLUDE_INITIAL_CONDITIONS_H_
#define INCLUDE_INITIAL_CONDITIONS_H_

#include "tangle.h"

#define DEG2RAD(X) ((X)*M_PI/180.0)

/*
	Setup the initial conditions for this simulation. Initial conditions can be set up in initial_conditions.c.
	@param tangle: Tangle holding all the informations about the vortices.
*/
void setup_init_conditions(struct tangle_state *tangle);

#endif /* INCLUDE_INITIAL_CONDITIONS_H_ */
