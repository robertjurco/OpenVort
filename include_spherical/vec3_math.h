/**
    Header for functions dealing with vector and matrices.
    @file vec3_math.h
*/

#ifndef VEC3_MATH_H
#define VEC3_MATH_H

/*******************************************************************************
 **************************** BASICS STRUCTURES ********************************
 ******************************************************************************/

 /**
	3 dimensional vector. Elements are stored in .p[i].
 */
struct vec3 {
  double p[3];
};

/**
   3x3 matrix. Elements are stored in .m[i][j].
*/
struct mat3 {
  double m[3][3];
};

/**
   Two vectors giving end points of the line. Vectors are stored in .r1 and .r2.
*/
struct segment {
  struct vec3 r1, r2;
};

/*******************************************************************************
 ************************ BASICS CONSTANT STRUCTURES ***************************
 ******************************************************************************/

extern const struct vec3 VEC_NULL;

extern const struct vec3 DIR_X;

extern const struct vec3 DIR_Y;

extern const struct vec3 DIR_Z;

extern const struct mat3 MAT_NULL;

extern const struct mat3 MAT_IDENTITY;

/*******************************************************************************
 **************************** VECTOR OPERSTIONS ********************************
 ******************************************************************************/

struct vec3 vec3(double x, double y, double z);

void vec3_assign(struct vec3 *v, double x, double y, double z);

void vec3_print(const struct vec3 *v);

double vec3_dot(const struct vec3 *u, const struct vec3 *v);

double vec3_ndot(const struct vec3 *u, const struct vec3 *v);

void vec3_cross(struct vec3 *res, const struct vec3 *u, const struct vec3 *v);

void vec3_mul(struct vec3 *res, const struct vec3 *u, double m);

void vec3_sub(struct vec3 *res, const struct vec3 *u, const struct vec3 *v);

void vec3_add(struct vec3 *res, const struct vec3 *u, const struct vec3 *v);

double vec3_len(const struct vec3 *u);

double vec3_dist(const struct vec3 *u, const struct vec3 *v);

void vec3_normalize(struct vec3 *v);

int vec3_equals(struct vec3 *u, struct vec3 *v);

struct vec3 segment_to_vec(const struct segment *seg);

/**
	Calculates the length of the segment.
	@param seg: Segment for which we want to calculate its length.
	@returns Returns the length of the given segment.
 */
static inline double segment_len(const struct segment *seg)
{
  return vec3_dist(&seg->r1, &seg->r2);
}

/*******************************************************************************
 *************************** SPHERICAL COORDINATES *****************************
 ******************************************************************************/

double radius(struct vec3* v);

double azimut_angle(struct vec3* v);

double polar_angle(struct vec3* v);

struct vec3 spherical_to_vector(double radius, double azim_angle, double polar_angle);

struct vec3 spherical_inversion(struct vec3* v, double radius_of_inversion);

void print_spherical(struct vec3* vec);

/*******************************************************************************
 **************************** MATRIX OPERSTIONS ********************************
 ******************************************************************************/

struct mat3 mat3(const struct vec3 *a, const struct vec3 *b, const struct vec3 *c);

void mat3_assign(struct mat3 *res, const struct vec3 *a, const struct vec3 *b, const struct vec3 *c);

void mat3_print(const struct mat3 *res);

void mat3_dmul(struct mat3 *res, const struct mat3 *M, double m);

void mat3_add(struct mat3 *res, const struct mat3 *a, const struct mat3 *b);

void mat3_sub(struct mat3 *res, const struct mat3 *a, const struct mat3 *b);

void mat3_mul(struct mat3 *res, const struct mat3 *a, const struct mat3 *b);

void mat3_vmul(struct vec3 *res, const struct mat3 *a, const struct vec3 *v);

void vec3_outer(struct mat3 *res, const struct vec3 *u, const struct vec3 *v);

#endif//VEC3_MATH_H
