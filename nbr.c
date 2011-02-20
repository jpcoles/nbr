#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nbr.h"
#include "ic.h"
#include "io.h"
#include "analysis.h"

void drift(struct env *env)
{
    int i;
    for (i=0; i < env->Nt; i++)
    {
        env->pt[i].r[0] += env->dt/2 * env->pt[i].v[0];
        env->pt[i].r[1] += env->dt/2 * env->pt[i].v[1];
        env->pt[i].r[2] += env->dt/2 * env->pt[i].v[2];
    }

    for (i=0; i < env->Nm; i++)
    {
        env->pm[i].r[0] += env->dt/2 * env->pm[i].v[0];
        env->pm[i].r[1] += env->dt/2 * env->pm[i].v[1];
        env->pm[i].r[2] += env->dt/2 * env->pm[i].v[2];
    }
}

void kickM(struct env *env)
{
    int i,j;
    double dx, dy, dz;
    double r, r2, r3;
    double vx, vy, vz;

    struct mparticle *p = env->pm;

    for (i=0; i < env->Nm; i++)
        p[i].phi = 0;

    for (i=0; i < env->Nm-1; i++)
    {
        for (j=i+1; j < env->Nm; j++)
        {
            dx = p[i].r[0] - p[j].r[0];
            dy = p[i].r[1] - p[j].r[1];
            dz = p[i].r[2] - p[j].r[2];

            r2 = dx*dx + dy*dy + dz*dz + env->eps2;
            r  = sqrt(r2);

            r3 = r*r2;

            vx = dx / r3 * env->dt;
            vy = dy / r3 * env->dt;
            vz = dz / r3 * env->dt;

            p[i].v[0] -= p[j].m * vx;
            p[i].v[1] -= p[j].m * vy;
            p[i].v[2] -= p[j].m * vz;
            p[i].phi  -= p[i].m * p[j].m / r;

            p[j].v[0] += p[i].m * vx;
            p[j].v[1] += p[i].m * vy;
            p[j].v[2] += p[i].m * vz;
            p[j].phi  -= p[i].m * p[j].m / r;
        }
    }

#if 0
    for (i=0; i < env->Nm; i++)
    {
        dx = p[i].r[0];
        dy = p[i].r[1];
        dz = p[i].r[2];

        r = sqrt(dx*dx + dy*dy + dz*dz);
        printf("%i %f %f %f %f\n", i, p[i].v[0], p[i].v[1], p[i].v[2], r);
    }
#endif
}

double pot(struct env *env)
{
    double t = env->t;

    return 1e11 * 0.5 * (1 - sin(10*t/env->T));
}

void kickT(struct env *env)
{
    int i,j;
    double dx, dy, dz;
    double r, r2, r3;

    struct mparticle *p = env->pm;
    struct tparticle *t = env->pt;

    for (i=0; i < env->Nt; i++)
        t[i].phi = 0;

    #pragma omp parallel for private(j, dx,dy,dz,r2,r,r3) shared(env, t,p)
    for (i=0; i < env->Nt; i++)
    {
        for (j=0; j < env->Nm; j++)
        {
            dx = t[i].r[0] - p[j].r[0];
            dy = t[i].r[1] - p[j].r[1];
            dz = t[i].r[2] - p[j].r[2];

            r2 = dx*dx + dy*dy + dz*dz + env->eps2;
            r  = sqrt(r2);

            r3 = r*r2;

            t[i].v[0] -= p[j].m * dx / r3 * env->dt;
            t[i].v[1] -= p[j].m * dy / r3 * env->dt;
            t[i].v[2] -= p[j].m * dz / r3 * env->dt;
            t[i].phi  -= p[j].m / r;
        }

#if 0
        dx = t[i].r[0];
        dy = t[i].r[1];
        dz = t[i].r[2];

        r2 = dx*dx + dy*dy + dz*dz;
        r  = sqrt(r2);
        r3 = r*r2;

        t[i].v[0] -= pot(env) * dx / r3 * env->dt;
        t[i].v[1] -= pot(env) * dy / r3 * env->dt;
        t[i].v[2] -= pot(env) * dz / r3 * env->dt;
#endif
    }
}

double time_step(struct env *env)
{
    return env->dt;
}

int main(int argc, char **argv)
{
    struct env env;
    double fft;
    int i;

    memset(&env, 0, sizeof(env));


    strcpy(env.tag, "frame");


    env.image.nr = 480;
    env.image.nc = 480;
    env.image.img = malloc(3 * env.image.nr * env.image.nc * sizeof(*env.image.img));

#if 0
    env.Nm = 2;
    env.Nt = 100000;
    float Rbin = 0.5;
    env.eps2 = 2*Rbin / 1.5;
    env.extent = 20;
    fft = tbr(&env, 1e11, 0.5, 2.0, 20.0, -2);
#endif

#if 0
    env.Nm = 2;
    env.Nt = 100000;
    float Rbin = 1.0;
    //env.eps2 = 0.5; //2*Rbin / 1.5;
    env.eps2 = 0.005; //pow(2*M_PI,-4/3.);
    printf("eps2 is %f\n", env.eps2);
    env.extent = Rbin*10;
    //fft = ball(&env, -0.8,0,0, 0.01, 1e11, Rbin, 0.2, 1, -2);
    fft = tbr(&env, 1, Rbin, .01*Rbin, env.extent, -1);
    //fft = single_mass(&env, 1, 1.1*Rbin, env.extent, 0);
#endif

#if 0
    // Fit very well after T = 8fft with einasto n=2
    env.Nm = 1;
    env.Nt = 10000;
    float Rbin = 1.0;
    env.eps2 = 0.005;
    printf("eps2 is %f\n", env.eps2);
    env.extent = Rbin*10;
    fft = single_mass(&env, 1, 1.1*Rbin, env.extent, -2);
#endif

#if 0
    env.Nm = 1;
    env.Nt = 10000;
    env.eps2 = 0.2;
    env.extent = 1;
    fft = single_mass(&env, 1e11, .2, 1.0, 0);
#endif
#if 0
    env.Nm = 1;
    env.Nt = 10000;
    env.eps2 = 0.2;
    env.extent = 1;
    fft = wall(&env, 1e11);
#endif
    env.T = 4*fft;
    env.Nsteps = 200;
    env.dt = 0.01;
    env.step = 0;

#if 1
    env.Nm = 1;
    env.Nt = 1000;
    env.eps2 = 0.005;
    printf("eps2 is %f\n", env.eps2);
    env.extent = 10;
    //fft = ball(&env, -0.8,0,0, 0.01, 1e11, Rbin, 0.2, 1, -2);
    fft = sgrstream(&env);
    fft = 2;
    env.T = 4*fft;
    env.Nsteps = 0;
    env.dt = 0.01;
    env.step = 0;
#endif

    printf("Maximum free-fall time is ~%.3e Gyr\n", fft);
    printf("Simulation time is ~%.3f Gyr\n", env.T);

    double *E0 = malloc(env.Nt * sizeof(*E0));
    double *E1 = malloc(env.Nt * sizeof(*E1));
    double *R0 = malloc(env.Nt * sizeof(*R0));

#if 1

    env.dt_step = env.T / env.Nsteps;
    double tt=0;

    energy(&env);
    for (i=0; i < env.Nt; i++)
    {
        R0[i] = sqrt(pow(env.pt[i].r[0],2)
                   + pow(env.pt[i].r[1],2)
                   + pow(env.pt[i].r[2],2));
        E0[i] = env.pt[i].E;
    }

    tprofile(&env, 0.01-1e-3, .1, 100, 1);

    for (env.t=0; env.t <= env.T; env.t += env.dt)
    {
        env.dt = time_step(&env);

        drift(&env);
        kickM(&env);
        kickT(&env);
        drift(&env);

        tt += env.dt;
        if (tt >= env.dt_step)
        {
            tt -= env.dt_step;

            tprofile(&env, 0.1, 10, 100, 1);

            fprintf(stderr, "%ld %f\n", env.step, env.t);
            //energy(&env);
            capture(&env);
            save_snapshot(&env);

            env.step++;
        }
    }

    if (tt)
    {
        tprofile(&env, 0.1, 10, 100, 1);

        fprintf(stderr, "%ld %f\n", env.step, env.t);
        capture(&env);
        save_snapshot(&env);
    }
#endif

    energy(&env);
    for (i=0; i < env.Nt; i++)
        E1[i] = env.pt[i].E;

    int Nbins=100;
    double *Pu = calloc(Nbins, sizeof(*Pu));
    double *Pd = calloc(Nbins, sizeof(*Pd));
    for (i=0; i < env.Nt; i++)
    {
        int b = R0[i] * Nbins;
        if (E1[i]-E0[i] >= 0)
            Pu[b]++;
        else
            Pd[b]++;
    }
    for (i=0; i < Nbins; i++)
    {
        double S = Pu[i] + Pd[i];
        Pu[i] /= S;
        Pd[i] /= S;
        printf("EBINS %e %e\n", Pu[i], Pd[i]);
    }

    double r;
    for (i=0; i < env.Nt; i++)
    {
        printf("ENERGY %f %f %f %f %f\n", E0[i], E1[i], E1[i]/E0[i], E1[i]-E0[i], R0[i]);
    }

    return 0;
}

