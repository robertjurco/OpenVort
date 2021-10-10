/**
	Functions dealing with vector and matrices.
	@file vec3_math.c
*/

#include <stdio.h>
#include <math.h>

#include "vec3_math.h"

/*******************************************************************************
 **************************** VECTOR OPERSTIONS ********************************
 ******************************************************************************/

 /**
	Assigment of the vector.
	@param x: x value of vector.
	@param x: y value of vector.
	@param x: z value of vector.
	@return Returns vector we assigned the values to.
 */
struct vec3 vec3(double x, double y, double z)
{
    struct vec3 out = {.p = {x, y, z}};
    return out;
}

/**
	Assigment of the vector.
	@param v: Vector to which we want values to assign.
	@param x: x value of vector.
	@param x: y value of vector.
	@param x: z value of vector.
 */
void vec3_assign(struct vec3 *v, double x, double y, double z)
{
    v->p[0] = x;
    v->p[1] = y;
    v->p[2] = z;
}

/**
	Prints the given vector.
	@param v: Vector we want to print.
 */
void vec3_print(const struct vec3 *v)
{
    printf("[%g, %g, %g]\n", v->p[0], v->p[1], v->p[2]);
}

/**
	Returns the dot product of two given vectors.
	@param u: First vector.
	@param v: Second vector.
	@return Returns the dot product of two given vectors.
 */
double vec3_dot(const struct vec3 *u, const struct vec3 *v)
{
    double out = 0;
    for(int k = 0; k < 3; ++k) out += u->p[k]*v->p[k];
    return out;
}

/**
	Returns the normalized dot product (cosine of the angle between them) of two given vectors.
	@param u: First vector.
	@param v: Second vector.
	@return Returns the normalized dot product (cosine of the angle between them) of two given vectors.
 */
double vec3_ndot(const struct vec3 *u, const struct vec3 *v)
{
    double dot = vec3_dot(u, v);
    double d1 = vec3_len(u);
    double d2 = vec3_len(v);

    return dot/d1/d2;
}

/**
	Calculate the cross product of two given vectors.
	@param res: Vector to which we assign the result.
	@param u: First vector.
	@param v: Second vector.
 */
void vec3_cross(struct vec3 *res, const struct vec3 *_u, const struct vec3 *_v)
{
    // We need copies in case res points to the same as _u or _v.
    struct vec3 u = *_u;
    struct vec3 v = *_v;
    res->p[0] = u.p[1] * v.p[2] - u.p[2] * v.p[1];
    res->p[1] = u.p[2] * v.p[0] - u.p[0] * v.p[2];
    res->p[2] = u.p[0] * v.p[1] - u.p[1] * v.p[0];
}

/**
	Multiplies the given vector with the given constant.
	@param res: Vector to which we assign the result.
	@param u: Vector we want to multiply.
	@param m: Constant we want to multiply vector with.
 */
void vec3_mul(struct vec3 *res, const struct vec3 *u, double m)
{
    for(int k=0; k<3; ++k) res->p[k] = u->p[k]*m;
}

/**
	Substracts two given vectors.
	@param res: Vector we assign the difference (result) to.
	@param u: Minuend vector (the one from which we subtract).
	@param v: Subtrahend vector (the one we subtract).
 */
void vec3_sub(struct vec3 *res, const struct vec3 *u, const struct vec3 *v)
{
	for(int k = 0; k < 3; ++k) res->p[k] = u->p[k] - v->p[k];
}

/**
	Adds two given vectors.
	@param res: Vector we assign the sum (result) to.
	@param u: First summand.
	@param v: Second summand.
 */
void vec3_add(struct vec3 *res, const struct vec3 *u, const struct vec3 *v)
{
    for(int k = 0; k < 3; ++k) res->p[k] = u->p[k] + v->p[k];
}

/**
	Returns the length of the given vector.
	@param u: Vector we want to find the length of.
	@returns The length of the given vector.
 */
double vec3_len(const struct vec3 *u)
{
    double out = 0;

    for(int k=0; k<3; ++k) out += u->p[k]*u->p[k];

    return sqrt(out);
}

/**
	Returns the distance between the two points given by two vectors.
	@param u: Position of the first point.
	@param v: Position of the second point.
	@returns The distance between the two points given by two vectors.
 */
double vec3_dist(const struct vec3 *u, const struct vec3 *v)
{
    struct vec3 vv;
    vec3_sub(&vv, u, v);
    return vec3_len(&vv);
}

/**
	Checks if the two given vectors are equal.
	@param u: First vector.
	@param v: Second vector.
	@returns Returns 1 if the given vectors are equal, returns 0 if not.
 */
int vec3_equals(struct vec3 *u, struct vec3 *v)
{
    int check = 0;

    if ((u->p[0] == v->p[0]) && (u->p[1] == v->p[1]) && (u->p[2] == v->p[2])) check = 1;

    return check;
}

/**
	Normalize the given vector.
	@param u: Vector you want to normalize. Result is then assign here.
 */
void vec3_normalize(struct vec3 *v)
{
    if(!vec3_equals(v, &VEC_NULL))
    {
        double d = vec3_len(v);
        vec3_mul(v, v, 1/d);
    } else {
		printf("ERROR: Vector you are trying to normalize is equal to the zero vector.");
	}
}

/**
	Creates the vector out of the segment by subtracting its endpoints.
	@param seg: Segment from which we want to create the vector.
	@returns Returns vector by subtracting seg.r2 - seg.r1.
 */
struct vec3 segment_to_vec(const struct segment *seg)
{
    struct vec3 v;
    vec3_sub(&v, &seg->r2, &seg->r1);

    return v;
}

/*******************************************************************************
 **************************** MATRIX OPERSTIONS ********************************
 ******************************************************************************/

 /**
   Assigment of the matrix.
   @param a: First column of the matrix (m[1][i] = a[i]).
   @param b: Second column of the matrix (m[2][i] = b[i]).
   @param c: Third column of the matrix (m[3][i] = c[i]).
   @return Returns matrix we assigned the values to.
*/
struct mat3 mat3(const struct vec3 *a, const struct vec3 *b, const struct vec3 *c)
{
    struct mat3 out;
    for(int i = 0; i < 3; i++)
    {
        out.m[i][0] = a->p[i];
        out.m[i][1] = b->p[i];
        out.m[i][2] = c->p[i];
    }
    return out;
}

/**
  Assigment of the matrix.
  @param res: Matrix we assign values to.
  @param a: First column of the matrix (m[1][i] = a[i]).
  @param b: Second column of the matrix (m[2][i] = b[i]).
  @param c: Third column of the matrix (m[3][i] = c[i]).
*/
void mat3_assign(struct mat3 *res, const struct vec3 *a, const struct vec3 *b, const struct vec3 *c)
{
    for(int i=0; i<3; i++)
    {
        res->m[i][0] = a->p[i];
        res->m[i][1] = b->p[i];
        res->m[i][2] = c->p[i];
    }
}

/**
	Prints the given matrix.
	@param res: Matrix we want to print.
 */
void mat3_print(const struct mat3 *res)
{
    printf("([%g, %g, %g],\n [%g, %g, %g],\n [%g, %g, %g])\n",
        res->m[0][0], res->m[0][1], res->m[0][2],
        res->m[1][0], res->m[1][1], res->m[1][2],
        res->m[2][0], res->m[2][1], res->m[2][2]);
}

/**
	Multiplies the given matrix with the constant.
	@param res: Result of the multiplication.
	@param m: Matrix we want to multiply.
	@param d: Number we mupltiply the matrix with.
 */
void mat3_dmul(struct mat3 *res, const struct mat3 *m, double d)
{
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            res->m[i][j] = d*m->m[i][j];
}

/**
	Substracts two given matrices.
	@param res: Matrix we assign the difference (result) to.
	@param a: Minuend matrix (the one from which we subtract).
	@param b: Subtrahend matrix (the one we subtract).
 */
void mat3_add(struct mat3 *res, const struct mat3 *a, const struct mat3 *b)
{
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            res->m[i][j] = a->m[i][j] + b->m[i][j];
}

/**
	Adds two given matrices.
	@param res: Matrix we assign the sum (result) to.
	@param a: First summand.
	@param b: Second summand.
 */
void mat3_sub(struct mat3 *res, const struct mat3 *a, const struct mat3 *b)
{
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            res->m[i][j] = a->m[i][j] - b->m[i][j];
}

/**
	Multiplies the two given matrices together.
	@param res: Matrix we assign the product (result) to.
	@param a: Multiplier matrix (the left one in the multiplication).
	@param b: Multiplicand matrix (the right one in the multiplication).
 */
void mat3_mul(struct mat3 *res, const struct mat3 *_a, const struct mat3 *_b)
{
    // We need copies in case res points to the same as _a or _b.
    struct mat3 a = *_a;
    struct mat3 b = *_b;
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            for(int k=0; k<3; ++k)
                res->m[i][j] = a.m[i][k]*b.m[k][j];
}

/**
	Multiplies the matrix by the vector.
	@param res: Matrix we assign the product (result) to.
	@param a: Matrix we multiply the vector with.
	@param b: Vector that is multiplied by the matrixs.
 */
void mat3_vmul(struct vec3 *res, const struct mat3 *a, const struct vec3 *_v)
{
	// We need copy in case res points to the same as _v.
	struct vec3 v = *_v;

    for(int i=0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
            res->p[i] = a->m[i][j] * v.p[j];
}

/**
	Outer product of the two given vectors.
	@param res: Matrix we assign the product (result) to.
	@param a: Multiplier vector (the left one in the multiplication).
	@param b: Multiplicand vector (the right one in the multiplication).
 */
void vec3_outer(struct mat3 *res, const struct vec3 *u, const struct vec3 *v)
{
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            res->m[i][j] = u->p[i]*v->p[j];
}


/*******************************************************************************
 ************************ BASICS CONSTANT STRUCTURES ***************************
 ******************************************************************************/

 /**
	Returns the null matrix.
	@returns Returns the null matrix.
 */
 struct mat3 mat3_NULL()
 {
     struct mat3 out;
     for(int i=0; i<3; i++)
        for (int j = 0; j < 3; j++) {
            out.m[i][j] = 0;
        }
     return out;
 }

 /**
	Null vector.
 */
const struct vec3 VEC_NULL = {{0, 0, 0}};

/**
  Unit vector pointing in the x direction (1,0,0).
*/
const struct vec3 VEC_X = {{1, 0, 0}};

/**
  Unit vector pointing in the y direction (0,1,0).
*/
const struct vec3 VEC_Y ={{0, 1, 0}};

/**
  Unit vector pointing in the z direction (0,0,1).
*/
const struct vec3 VEC_Z = {{0, 0, 1}};

/**
   Null matrix.
*/
const struct mat3 MAT_NULL = {{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};

/**
   Identity matrix.
*/
const struct mat3 MAT_IDENTITY = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
