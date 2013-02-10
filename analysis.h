#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "nbr.h"

void tprofile(struct env *env, float Rmin, float Rmax, int Nbins, int use_log);
void energy(struct env *env);
void CoM(struct env *env);
void angmom(struct env *env);

#endif

