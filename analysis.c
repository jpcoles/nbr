#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "analysis.h"
#include "io.h"

void tprofile(struct env *env, float Rmin, float Rmax, int Nbins, int use_log)
{
    FILE *fp;
    char *fname;

    int i,b;
    long M;

    float r,V;
    float *bins = calloc(Nbins, sizeof(*bins));

    if (use_log)
    {
        Rmin = log10(Rmin);
        Rmax = log10(Rmax);
    }

    for (i=0; i < env->Nt; i++)
    {
        r = sqrt(pow(env->pt[i].r[0],2)
               + pow(env->pt[i].r[1],2)
               + pow(env->pt[i].r[2],2));

        if (use_log)
            r = log10(r);

        if (Rmin < r && r <= Rmax) 
        {
            b = ((r - Rmin) / (Rmax-Rmin)) * Nbins - 1;
            if (b >= Nbins)
                printf("!! %i %f %e %e\n", b, r, Rmin, Rmax);
            assert(b >= 0);
            assert(b < Nbins);

            bins[b]++;
        }

    }

    fname = make_fname(env, "den", NULL, 0);
    fp = fopen(fname, "wt");

    M = 0;
    time_t t = time(NULL);
    fprintf(fp, "# Density profile (R,rho(R)), %s\n", use_log ? "log" : "linear");
    fprintf(fp, "# %s\n", ctime(&t));
    for (i=0; i < Nbins; i++)
    {
        M += bins[i];
        r = (i+1) * (Rmax-Rmin) / Nbins + Rmin;
        if (use_log)
            r = pow(10,r);

        V = 4./3.*M_PI*pow(r,3);

        fprintf(fp, "%e %e %ld\n", r, bins[i]/V, M);
    }

#if 1
    M = 0;
    double r_e=0;
    double rho_e=0;
    double norm = 0;
    for (i=0; i < Nbins; i++)
    {
        r = (i+1) * (Rmax-Rmin) / Nbins + Rmin;
        if (use_log)
            r = pow(10,r);
        V = 4./3.*M_PI*pow(r,3);

        M += bins[i];
        if (M >= env->Nt/2)
        {
            r_e = r;
            rho_e = bins[i] / V;
            break;
        }
    }
    if (rho_e)
    {
        double n=3;
        double dn = 3*n - 1./3 + 0.0079/n;

        double chi2 = 0;

        for (i=0; i < Nbins; i++)
        {
            r = (i+1) * (Rmax-Rmin) / Nbins + Rmin;
            if (use_log)
                r = pow(10,r);
            V = 4./3.*M_PI*pow(r,3);

            double ein = rho_e * exp(-dn*(pow(r/r_e,1./n) - 1));
	    norm += ein;
	}

        for (i=0; i < Nbins; i++)
        {
            r = (i+1) * (Rmax-Rmin) / Nbins + Rmin;
            if (use_log)
                r = pow(10,r);
            V = 4./3.*M_PI*pow(r,3);

            double ein = rho_e * exp(-dn*(pow(r/r_e,1./n) - 1));
            double rho = bins[i] / V;

            chi2 += pow(rho - ein,2) / ein;
        }

        printf("Chi2 is %f\n", chi2 / norm);
    }
#endif

    free(fname);

    free(bins);
}

void energy(struct env *env)
{
    int i,j;
    double E=0;
    double v2;
    double phi;
    double dx, dy, dz, r2,r;

    for (i=0; i < env->Nm; i++)
    {
        v2 = pow(env->pm[i].v[0],2)
           + pow(env->pm[i].v[1],2)
           + pow(env->pm[i].v[2],2);

        E += 0.5*env->pm[i].m * v2 + env->pm[i].phi;
    }

    for (i=0; i < env->Nt; i++)
    {
	phi = 0;
        for (j=0; j < env->Nm; j++)
        {
            dx = env->pt[i].r[0] - env->pm[j].r[0];
            dy = env->pt[i].r[1] - env->pm[j].r[1];
            dz = env->pt[i].r[2] - env->pm[j].r[2];

            r2 = dx*dx + dy*dy + dz*dz + env->eps2;
            r  = sqrt(r2);

            phi -= env->pm[j].m / r;
	}

        v2 = pow(env->pt[i].v[0],2)
           + pow(env->pt[i].v[1],2)
           + pow(env->pt[i].v[2],2);

	env->pt[i].E = 0.5*v2 + phi;
    }

    //printf("ENERGY %f\n", E);
}

