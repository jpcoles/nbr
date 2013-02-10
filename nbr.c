#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nbr.h"
#include "ic.h"
#include "io.h"
#include "analysis.h"

void drift(struct env *env, double dt)
{
    int i;
    for (i=0; i < env->Nt; i++)
    {
        env->pt[i].r[0] += dt * env->pt[i].v[0];
        env->pt[i].r[1] += dt * env->pt[i].v[1];
        env->pt[i].r[2] += dt * env->pt[i].v[2];
    }

    for (i=0; i < env->Nm; i++)
    {
        env->pm[i].r[0] += dt * env->pm[i].v[0];
        env->pm[i].r[1] += dt * env->pm[i].v[1];
        env->pm[i].r[2] += dt * env->pm[i].v[2];
    }
}

void accelM(struct env *env)
{
    int i,j;
    double dx, dy, dz;
    double r, r2, r3;
    double ax, ay, az;
    double rinv;

    struct mparticle *p = env->pm;

    for (i=0; i < env->Nm; i++)
    {
        p[i].a[0] = 0;
        p[i].a[1] = 0;
        p[i].a[2] = 0;
        p[i].phi = 0;
    }

    for (i=0; i < env->Nm-1; i++)
    {
        for (j=i+1; j < env->Nm; j++)
        {
            dx = p[i].r[0] - p[j].r[0];
            dy = p[i].r[1] - p[j].r[1];
            dz = p[i].r[2] - p[j].r[2];

            r2 = dx*dx + dy*dy + dz*dz + env->eps2;
            rinv  = 1./sqrt(r2);

            //double rhatx = dx / r;
            //double rhaty = dy / r;
            //double rhatz = dz / r;

            //r3 = r*r2;

            //ax = (rhatx / r2);
            //ay = (rhaty / r2);
            //az = (rhatz / r2);
            ax = dx * pow(rinv,3);
            ay = dy * pow(rinv,3);
            az = dz * pow(rinv,3);

            //fprintf(stderr, "%e %e %e %e %e\n", vx, vy, vz, p[j].m, p[j].m * vx);

            p[i].a[0] -= p[j].m * ax;
            p[i].a[1] -= p[j].m * ay;
            p[i].a[2] -= p[j].m * az;

            p[j].a[0] += p[i].m * ax;
            p[j].a[1] += p[i].m * ay;
            p[j].a[2] += p[i].m * az;
        }
    }
}

void kick(struct env *env, double dt)
{
    int i;
//  for (i=0; i < env->Nt; i++)
//  {
//      env->pt[i].v[0] += dt * env->pt[i].a[0];
//      env->pt[i].v[1] += dt * env->pt[i].a[1];
//      env->pt[i].v[2] += dt * env->pt[i].a[2];
//  }

    for (i=0; i < env->Nm; i++)
    {
        env->pm[i].v[0] += dt * env->pm[i].a[0];
        env->pm[i].v[1] += dt * env->pm[i].a[1];
        env->pm[i].v[2] += dt * env->pm[i].a[2];
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

double pot(struct env *env)
{
    double t = env->t;

    return 1e11 * 0.5 * (1 - sin(10*t/env->T));
}

void analytic_potential_kick(struct env *env, struct mparticle *p, struct tparticle *t)
{
    int i;
    double dx, dy, dz;
    double r, r2, r3;
    if (p != NULL)
    {
    }

    if (t != NULL)
    {
        dx = t->r[0];
        dy = t->r[1];
        dz = t->r[2];

        r2 = dx*dx + dy*dy + dz*dz + env->eps2;
        r  = sqrt(r2);

        r3 = r*r2;

        t->v[0] -= 0;
        t->v[1] -= 0;
        t->v[2] -= 0;
        t->phi  -= 0;
    }
}

void accelT(struct env *env, double dt)
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

            t[i].v[0] -= p[j].m * dx / r3 * dt;
            t[i].v[1] -= p[j].m * dy / r3 * dt;
            t[i].v[2] -= p[j].m * dz / r3 * dt;
            t[i].phi  -= p[j].m / r;
        }

        analytic_potential_kick(env, NULL, &t[i]);

        t[i].v[0] *= env->dt;
        t[i].v[1] *= env->dt;
        t[i].v[2] *= env->dt;

#if 0
        dx = t[i].r[0];
        dy = t[i].r[1];
        dz = t[i].r[2];

        r2 = dx*dx + dy*dy + dz*dz;
        r  = sqrt(r2);
        r3 = r*r2;

        t[i].v[0] -= pot(env) * dx / r3 * dt;
        t[i].v[1] -= pot(env) * dy / r3 * dt;
        t[i].v[2] -= pot(env) * dz / r3 * dt;
#endif
    }
}

double time_step(struct env *env)
{
    //fprintf(stderr, "%f\n", exp(-0.5*env->t));
    //return (env->dt-0.001/env->units.T) * (1. / (1. + exp(-0.5*env->t))) + 0.001/env->units.T;
    return env->dt;
}

void process_step(struct env *env)
{
    tprofile(env, 0.1, 10, 100, 1);
    energy(env);
    //CoM(env);
    env->CoM[0] = env->pm[0].r[0];
    env->CoM[1] = env->pm[0].r[1];
    env->CoM[2] = env->pm[0].r[2];
    angmom(env);
    capture(env);
    save_snapshot(env);
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
    env.image.img = calloc(3 * env.image.nr * env.image.nc, sizeof(*env.image.img));

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

#if 0
    env.T = 4*fft;
    env.Nsteps = 200;
    env.dt = 0.01;
    env.step = 0;
#endif

#if 0
    droplet(&env, 1,1,1, 0.1, 0.001, 0.01);
    env.T = env.units.T * 4;
    env.Nsteps = 200;
    env.dt = (env.T / env.Nsteps) / 100;
    env.step = 0;
#endif

#if 0
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

#if 1
    env.Nm = 13;
    env.Nt = 0;
    env.eps2 = 0.000;
    fft = solsystem(&env);

    env.extent = 100 / env.units.L; //5.0 / env.units.L;
    //env.T = 10 / env.units.T; //40 * 365;
    env.T = 14611 / env.units.T; //40 * 365;
    env.Nsteps = 200;
    env.dt = pow(2,-7) / env.units.T;
    //env.dt = 0.0001 / env.units.T;
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

#if 0
    energy(&env);
    for (i=0; i < env.Nt; i++)
    {
        R0[i] = sqrt(pow(env.pt[i].r[0],2)
                   + pow(env.pt[i].r[1],2)
                   + pow(env.pt[i].r[2],2));
        E0[i] = env.pt[i].E;
    }
#endif

    //tprofile(&env, 0.01-1e-3, .1, 100, 1);

    double dt=0;
    env.t = 0;
    dt = time_step(&env);

    //energy(&env);
    //angmom(&env);

    process_step(&env);

    accelM(&env);
    kick(&env,  dt/2);

    for (env.t=0; env.t <= env.T; env.t += dt)
    {
        dt = time_step(&env);

        drift(&env, dt);
        accelM(&env);
        //accelT(&env);

        tt += dt;
        if (tt >= env.dt_step)
        {
            tt -= env.dt_step;
            env.step++;

            fprintf(stderr, "%ld %f %f\n", env.step, env.t, dt);

            kick(&env, dt/2.);
            process_step(&env);
            kick(&env, dt/2.);

        }
        else
        {
            kick(&env, dt);
        }
    }

    if (tt)
    {
        kick(&env, -dt/2.);
        fprintf(stderr, "%ld %f\n", env.step, env.t);
        process_step(&env);
    }
#endif

    energy(&env);
    for (i=0; i < env.Nt; i++)
        E1[i] = env.pt[i].E;

    int Nbins=100;
    double *Pu = calloc(Nbins, sizeof(*Pu));
    double *Pd = calloc(Nbins, sizeof(*Pd));

    if (env.Nt)
    {
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
    }

    double L = env.units.L;
    double V = env.units.L / env.units.T;
    for (i=0; i < env.Nm; i++)
    {
#if 0
        printf("PS %i %18.15e %18.15e %18.15e %18.15e %18.15e %18.15e %18.15e\n", 
                i, env.pm[i].r[0], env.pm[i].r[1], env.pm[i].r[2],
                   env.pm[i].v[0], env.pm[i].v[1], env.pm[i].v[2], env.pm[i].m);

#else
        printf("PS %i %18.15e %18.15e %18.15e %18.15e %18.15e %18.15e\n", 
                i, env.pm[i].r[0]-env.pm[0].r[0], env.pm[i].r[1]-env.pm[0].r[1], env.pm[i].r[2]-env.pm[0].r[2],
                   env.pm[i].v[0]-env.pm[0].v[0], env.pm[i].v[1]-env.pm[0].v[1], env.pm[i].v[2]-env.pm[0].v[2]);
#endif
    }

    return 0;
}

