#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "nbr.h"

void tprofile(struct env *env, float Rmin, float Rmax, int Nbins, int use_log);
void energy(struct env *env);

#endif

