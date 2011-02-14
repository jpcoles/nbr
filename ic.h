#ifndef IC_H
#define IC_H

double tbr(struct env *env, float M, float Rbin, float Rmin, float Rmax, float slope);
double ball(struct env *env, 
            float x0, float y0, float z0, float R,
            float M, float Rbin, float Rmin, float Rmax, float slope);

double single_mass(struct env *env, float M, float Rmin, float Rmax, float slope);
double wall(struct env *env, float M);
#endif
