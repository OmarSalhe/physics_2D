#ifndef PHYS_H
#define PHYS_H

#include "math.h"

/*
* Vectors:
*   A discrete struct comprised of components, with each being respect to an independent dimension. 
*   In this case, only 2D motion will be done, so all vectors will comprised of two components 
*/
typedef struct{
    float x, y;
} Vector;

/*
* Calculates the magnitude of one vectors projection onto another vector
* The dot product is derived from the eq: a • b = (a.x * b.x) + (a.y * b.y)
*/
float dot_product(Vector *a, Vector *b){
    return (a->x * b->x) + (a->y * b->y); 
}

/*
* Scales a vector's magnitude by a given scalar (AKA constant)
* v = k * u = <k * u.x, k * u.y>
* where:
* - v = scaled vector
* - u = initial vector
* - k = scalar
*/
Vector scalar_product(Vector *a, float k){ // increases or decreases the magnitude of a vector by a factor of k, where k is a constant
    return (Vector){a->x * k, a->y*k};
}

/*
* Calcuates the magnitude of a vector using the 2D distance equation between the components and the origin of a 2D coordinate plane
* The 2D distance equation is defined as: d = sqrt((a.x - b.x) ^ 2 - (a.y - b.y) ^ 2)
* where:
* - d = distance
* - a = intial vector
* - b = the reference vector (in this case, the origin vector, i.e., <0,0>) from which the distance of vector 'a' is measured
*/
float mag(Vector *a){
    return sqrt((a->x * a->x) + (a->y * a->y));
}

/*
* Reflects vector 'a' off of a reflective surface defined by vector 'b'
* The reflection is derived from the reflection thrm: r = a - 2 * proj_n(a)
* where:
* - r = reflected vector
* - a = initial vector
* - proj_n(a) = the projection of vector 'a' onto vector 'n' (how much of a aligns with n)
* - n = the normal vector to the reflective surface, derived from b
*/
Vector reflect(Vector *a, Vector *b){
    Vector n = {-b->y, b->x};
    Vector tmp = scalar_product(&n, pow(mag(b), -1)); // intermediary normalization
    Vector c = scalar_product(&tmp, -2 * (dot_product(a, &tmp) / dot_product(&tmp, &tmp)));
    return (Vector){a->x - c.x, a->y - c.y};
}

float angle_between(Vector *a, Vector *b){
    /*
    * Returns the angle formed between two vectors
    * Derived from the cosine law and dot product, the formula: a • b = |a|*|b| * cos(theta)
    * allows the angle to be defined as: theta = arccos(a • b / (|a|*|b|))
    * where:
    * - a, b = given vectors
    * - theta = the angle formed between a, b
    * - || = shorthand for the magnitude of a vector
    */
    return dot_product(a, b) / (mag(a) * mag(b));
}

float cross_product(Vector *a, Vector *b){
    /*
    * Returns the determinant of the matrix formed by two vectors
    * If a = <a.x, a.y>, b = <b.x, b.y>, a
    * 
    * 
    */
}

// Particle
typedef struct {
    Vector pos;
    Vector velocity;
    Vector force;
} Particle;

/*
* Solvers
* Depending on the scenario being simulated, the differential equation describing it
* can either be solved directly, giving an accurate position equation, or,
* through numerical methods like Runge-Kutta 4th order, a close approximation of the particles position
* can be found
*
* These two approaches define the two solvers used: Analytical Solver, Numerical Solver
* One is a direct process from equation to solution, whereas the other is iterative and more computationally intensive
*/

float analytical_solver(float alpha, float gamma, ){
    
}

#endif