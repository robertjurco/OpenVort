#ifndef EXTERNAL_VELOCITY_H
#define EXTERNAL_VELOCITY_H

#include "vec3_math.h"

/*
	Structure holding all the information about the external velocity of the fluid.
	@param type: NO_FLOW, SIMPLE, DECREASING
	@param strength: Strength of constant SIMPLE fluid velocity AND final constant strength for the DECREASING flow type.
	@param time_of_linear_decrease: Time for which the DECREASING flow type is lineary decreasing from its starting strength onto strength.
	@param starting_strength: Starting strength for the DECREASING flow type.
*/
struct ext_vel_param {
	int type;
	double strength;
	double time_of_linear_decrease;
	double starting_strength;
	struct vec3 direction;
};


/*
	Types of external velocities.
*/
typedef enum _ev_t {
	NO_FLOW, // External veloity: no flow.
	SIMPLE, // External veloity: constant in one direction given by vector.
	DECREASING, // External veloity: decreases from given value and in given time to a constant in one direction given by vector.
} external_velocity_t;

extern struct ext_vel_param vn_conf;

extern struct ext_vel_param vs_conf;

extern struct ext_vel_param vb_conf;

int get_vn(const struct vec3 *where, double t, struct vec3 *res);

int get_vs(const struct vec3 *where, double t, struct vec3 *res);

int get_vb(const struct vec3 *where, double t, struct vec3 *res);

#endif /* EXTERNAL_VELOCITY_H */
