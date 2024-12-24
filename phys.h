#ifndef PHYS_H
#define PHYS_H

#include "math.h"

typedef struct{
    float x, y;
} Vector; // definition of a vector

float dot_product(Vector *a, Vector *b){
    /*
    * Calculates the magnitude of one vectors projection onto another vector
    * The dot product is derived from the eq: a * b = (a.x * b.x) + (a.y * b.y)
    */
    return (a->x * b->x) + (a->y * b->y); 
}

Vector scalar_product(Vector *a, float k){ // increases or decreases the magnitude of a vector by a factor of k, where k is a constant
    /*
    * Scales a vector's magnitude by a given constant (AKA scalar)
    * v = k * u = <k * u.x, k * u.y>
    * where:
    * - v = scaled vector
    * - u = initial vector
    * - k = given constant
    */
    return (Vector){a->x * k, a->y*k};
}

float mag(Vector *a){ // the 'size' or 'strength' of a vector (for all intents and purposes), found via finding the distance between the components on a 2D coordinate plane and the origin
    /*
    * Calcuates the magnitude of a vector using the 2D distance equation between the components and the origin of a 2D coordinate plane
    * Distance Eq: d = sqrt((a.x - b.x) ^ 2 - (a.y - b.y) ^ 2)
    * where:
    * - d = distance
    * - a = intial vector
    * - b = the vector we're finding the distance from (in this case the 0 vector)
    */
    return sqrt((a->x * a->x) + (a->y * a->y));
}

Vector reflect(Vector *a, Vector *b){
    /*
    * Reflects vector 'a' off of a reflective surface defined by vector 'b'
    * The reflection is derived from the reflection thrm: r = a - 2 * proj_n(a)
    * where:
    * - r = reflected vector
    * - a = initial vector
    * - proj_n(a) = the projection of vector 'a' onto vector 'n' (how much of a aligns with n)
    * - n = the normal vector to the reflective surface, derived from b
    */
    Vector n = {-b->y, b->x};
    Vector tmp = scalar_product(&n, pow(mag(b), -1)); // intermediary norm calculation
    Vector c = scalar_product(&tmp, -2 * (dot_product(a, &tmp) / dot_product(&tmp, &tmp)));
    return (Vector){a->x - c.x, a->y - c.y};
}

typedef struct {
    Vector pos;
} Particle;
#endif