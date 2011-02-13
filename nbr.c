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

void clear_velocity(struct env *env)
{
    int i;

    struct mparticle *p = env->pm;
    struct tparticle *t = env->pt;

    for (i=0; i < env->Nt; i++)
    {
        t[i].v[0] = 
        t[i].v[1] = 
        t[i].v[2] = 0;
    }

    for (i=0; i < env->Nm; i++)
    {
        p[i].v[0] = 
        p[i].v[1] = 
        p[i].v[2] = 0;
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

    #pragma omp parallel for private(j, dx,dy,dz,r2,r,r3) shared(env, t,p)
    for (i=0; i < env->Nt; i++)
    {
        for (j=0; j < env->Nm; j++)
        {
            dx = t[i].r[0] - p[j].r[0];
            dy = t[i].r[1] - p[j].r[1];
            dz = t[i].r[2] - p[j].r[2];

            r2 = dx*dx + dy*dy + dz*dz + p[j].eps2;
            r  = sqrt(r2);

            r3 = r*r2;

            t[i].v[0] -= p[j].m * dx / r3 * env->dt;
            t[i].v[1] -= p[j].m * dy / r3 * env->dt;
            t[i].v[2] -= p[j].m * dz / r3 * env->dt;
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
#else
    env.Nm = 2;
    env.Nt = 10000;
    float Rbin = 0.1;
    env.eps2 = 1.05; //2*Rbin / 1.5;
    env.extent = 1;
    fft = ball(&env, -0.8,0,0, 0.001, 1e11, Rbin, 0.2, 1, -2);
    //fft = tbr(&env, 1e11, Rbin, .2, 1.0, -2);
#endif
    env.T = 4 * fft;
    env.Nsteps = 200;
    env.dt = 0.00001;
    env.step = 0;

    printf("Maximum free-fall time is ~%.3e Gyr\n", fft);
    printf("Simulation time is ~%.3f Gyr\n", env.T);

#if 1

    env.dt_step = env.T / env.Nsteps;
    double tt=0;

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
        //energy(&env);
        capture(&env);
        save_snapshot(&env);
    }
#endif

    return 0;
}
