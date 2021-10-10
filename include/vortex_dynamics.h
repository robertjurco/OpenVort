#ifndef VORTEX_DYNAMICS
#define VORTEX_DYNAMICS

#include "tangle.h"

// Stepping
struct rk4_state {
    //the tangle states here should be reduced
    //only calculate what we really need for the velocity
    struct tangle_state *s1, *s2, *s3, *s4;

    double dt;
};

void rk4_step(struct tangle_state *tangle, double t, double dt);
//reconnections
int reconnect(struct tangle_state *tangle, double t, double rec_dist, double rec_angle);

#endif //VORTEX_DYNAMICS
