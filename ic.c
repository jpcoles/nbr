#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "nbr.h"

#define MODE_NORMAL              0
#define MODE_BINARY              1
#define MODE_STATIC_BINARY       2
#define MODE_RMAX_SHELL          4
#define MODE_DROPLET             8
#define MODE_2D                  16
#define MODE_XY                  32
#define MODE_YZ                  64
#define MODE_XZ                  128
#define ALL_MODES               255

void units(double *M, double *L, double *T)
{
    double Myr  = 1e6 * 365 * 24 * 60 * 60;   // [s]
    double Gyr  = 1e9 * 365 * 24 * 60 * 60;   // [s]
    double kpc  = 3.08568025e19;              // [m]
    double Msun = 1.98892e30;                 // [kg]

    double Gsi = 6.67300e-11*pow(kpc,-3)*Msun*pow(Gyr,2); // [kpc^3 msun^-1 gyr^-2]

    *M = 2.3262e5; // [msun]
    *L = 1.0     ; // [kpc]
    *T = 1.0 / sqrt(Gsi * pow(*L,-3) * *M);
}

double new_point(double Rmin, double Rmax, double slope)
{
    double c, r;

    // Normalization constant
    c = (slope == -3)
      ? log(Rmax/Rmin)
      : (pow(Rmax,3+slope) - pow(Rmin,3+slope)) / (3+slope);

    do
    {
        r = (slope == -3)
          ? Rmin * exp(drand48()*c)
          : pow((3+slope)*drand48()*c + pow(Rmin, 3+slope), 1/(3+slope));
    } while (r < Rmin);

    return r;
}

double enc_mass(double r, double M, double Rmin, double Rmax, double slope)
{
    // Normalization constant
    double c = (slope == -3)
             ? log(Rmax/Rmin)
             : (pow(Rmax,3+slope) - pow(Rmin,3+slope)) / (3+slope);

    double m  = (slope == -3)
              ? log(r/Rmin) * M/c
              : ((pow(r,3+slope) - pow(Rmin,3+slope)) / (3+slope)) * M/c; 
    return m;
}

double ball(struct env *env, 
            float x0, float y0, float z0, float R,
            float M, float Rbin, float Rmin, float Rmax, float slope)
{
    int i;

    double em,ff=0,FFmax=0;

    double x,y,z,t,w;

    double Mt, Lt, Tt;
    units(&Mt, &Lt, &Tt);

    //assert(env->Nm == 2);

    if (env->Nm) env->pm = malloc(env->Nm * sizeof(*env->pm));
    if (env->Nt) env->pt = malloc(env->Nt * sizeof(*env->pt));

    double v, mu;
    double m = M / env->Nm;
    double r,d;

    for (i=0; i < env->Nt; i++)
    {
#if 0
        if (prop->mode & MODE_RMAX_SHELL) 
            r = Rmax;
        else
#endif
        r = new_point(R/100, R, slope);

        z = 2.0 * drand48() - 1.0;

        t = 2.0 * M_PI * drand48();
        w = sqrt(1 - z*z);
        x = w * cos(t);
        y = w * sin(t);


        env->pt[i].r[0] = (x0 + x * r) / Lt;
        env->pt[i].r[1] = (y0 + y * r) / Lt;
        env->pt[i].r[2] = (z0 + z * r) / Lt;
        //printf("%f %f %f\n", env->pt[i].r[0], env->pt[i].r[1], env->pt[i].r[2]);

        env->pt[i].v[0] = 0 / (Lt/Tt);
        env->pt[i].v[1] = 0 / (Lt/Tt);
        env->pt[i].v[2] = 0 / (Lt/Tt);

        r = sqrt(pow(x0 + x * r, 2)
               + pow(y0 + y * r, 2)
               + pow(z0 + z * r, 2));

        em = enc_mass(r, M, Rmin, Rmax, slope);
        ff = M_PI/2 * sqrt((pow(r / Lt,3))/ (2*em/Mt));
        //ff = M_PI/4 * sqrt((2*pow(r / Lt,3))/ (em/Mt));
        if (ff > FFmax) 
        {
            //cerr << ff << " " << r << " " << (M) << endl;
            FFmax = ff;
        }
    }

    if (env->Nm)
    {
        x = Rbin/Lt;
        d = 2*x;
        v = sqrt(m/Mt * d*x/pow(d*d + env->eps2, 1.5));

        env->pm[0].r[0] = Rbin / Lt;
        env->pm[0].r[1] = 0 / Lt;
        env->pm[0].r[2] = 0 / Lt;
        env->pm[0].v[0] = 0 / (Lt/Tt);
        env->pm[0].v[1] = v;
        env->pm[0].v[2] = 0 / (Lt/Tt);
        env->pm[0].m      = m / Mt;

        env->pm[1].r[0] = -env->pm[0].r[0];
        env->pm[1].r[1] =  env->pm[0].r[1];
        env->pm[1].r[2] =  env->pm[0].r[2];
        env->pm[1].v[0] =  env->pm[0].v[0];
        env->pm[1].v[1] = -env->pm[0].v[1];
        env->pm[1].v[2] =  env->pm[0].v[2];
        env->pm[1].m    =  env->pm[0].m;
    }

    return FFmax;
}

double tbr(struct env *env, float M, float Rbin, float Rmin, float Rmax, float slope)
{
    int i;

    double em,ff=0,FFmax=0;

    double x,y,z,t,w;

    double Mt=1, Lt=1, Tt=1;
    //units(&Mt, &Lt, &Tt);

    //assert(env->Nm == 2);

    if (env->Nm) env->pm = malloc(env->Nm * sizeof(*env->pm));
    if (env->Nt) env->pt = malloc(env->Nt * sizeof(*env->pt));

    double v, mu;
    double m = M / env->Nm;
    double r,d;

    for (i=0; i < env->Nt; i++)
    {
#if 0
        if (prop->mode & MODE_RMAX_SHELL) 
            r = Rmax;
        else
#endif
        r = new_point(Rmin, Rmax, slope);

        z = 2.0 * drand48() - 1.0;

        t = 2.0 * M_PI * drand48();
        w = sqrt(1 - z*z);
        x = w * cos(t);
        y = w * sin(t);


        env->pt[i].r[0] = x * (r / Lt);
        env->pt[i].r[1] = y * (r / Lt);
        env->pt[i].r[2] = z * (r / Lt);
        //printf("%f %f %f\n", env->pt[i].r[0], env->pt[i].r[1], env->pt[i].r[2]);

        env->pt[i].v[0] = 0 / (Lt/Tt);
        env->pt[i].v[1] = 0 / (Lt/Tt);
        env->pt[i].v[2] = 0 / (Lt/Tt);

        em = enc_mass(r, M, Rmin/100, Rmax, slope);
        //em = enc_mass(r, M, Rmin, Rmax, slope);
        ff = M_PI/4 * sqrt((2*pow(r / Lt,3))/ (em/Mt));
        //printf("%f, %f, %f\n", ff, r, em);
        if (ff > FFmax) 
        {
            //cerr << ff << " " << r << " " << (M) << endl;
            FFmax = ff;
        }
    }

    if (env->Nm)
    {
        x = Rbin/Lt;
        d = 2*x;
        v = sqrt(m/Mt * d*x/pow(d*d + env->eps2, 1.5));

        env->pm[0].r[0] = Rbin / Lt;
        env->pm[0].r[1] = 0 / Lt;
        env->pm[0].r[2] = 0 / Lt;
        env->pm[0].v[0] = 0 / (Lt/Tt);
        env->pm[0].v[1] = v;
        env->pm[0].v[2] = 0 / (Lt/Tt);
        env->pm[0].m      = m / Mt;

        env->pm[1].r[0] = -env->pm[0].r[0];
        env->pm[1].r[1] =  env->pm[0].r[1];
        env->pm[1].r[2] =  env->pm[0].r[2];
        env->pm[1].v[0] =  env->pm[0].v[0];
        env->pm[1].v[1] = -env->pm[0].v[1];
        env->pm[1].v[2] =  env->pm[0].v[2];
        env->pm[1].m    =  env->pm[0].m;
    }

    return FFmax;
}

double single_mass(struct env *env, float M, float Rmin, float Rmax, float slope)
{
    int i;

    double em,ff=0,FFmax=0;

    double x,y,z,t,w;

    double Mt=1, Lt=1, Tt=1;
    //units(&Mt, &Lt, &Tt);

    assert(env->Nm == 1);

    if (env->Nm) env->pm = malloc(env->Nm * sizeof(*env->pm));
    if (env->Nt) env->pt = malloc(env->Nt * sizeof(*env->pt));

    double v, mu;
    double m = M / env->Nm;
    double r,d;

    for (i=0; i < env->Nt; i++)
    {
#if 0
        if (prop->mode & MODE_RMAX_SHELL) 
            r = Rmax;
        else
#endif
        r = new_point(Rmin, Rmax, slope);

        z = 2.0 * drand48() - 1.0;

        t = 2.0 * M_PI * drand48();
        w = sqrt(1 - z*z);
        x = w * cos(t);
        y = w * sin(t);


        env->pt[i].r[0] = x * (r / Lt);
        env->pt[i].r[1] = y * (r / Lt);
        env->pt[i].r[2] = z * (r / Lt);
        //printf("%f %f %f\n", env->pt[i].r[0], env->pt[i].r[1], env->pt[i].r[2]);

        env->pt[i].v[0] = 0 / (Lt/Tt);
        env->pt[i].v[1] = 0 / (Lt/Tt);
        env->pt[i].v[2] = 0 / (Lt/Tt);

        em = enc_mass(r, M, Rmin/100, Rmax, slope);
        //em = enc_mass(r, M, Rmin, Rmax, slope);
        ff = M_PI/4 * sqrt((2*pow(r / Lt,3))/ (em/Mt));
        //printf("%f, %f, %f\n", ff, r, em);
        if (ff > FFmax) 
        {
            //cerr << ff << " " << r << " " << (M) << endl;
            FFmax = ff;
        }
    }

    if (env->Nm)
    {
        env->pm[0].r[0] = 0;
        env->pm[0].r[1] = 0 / Lt;
        env->pm[0].r[2] = 0 / Lt;
        env->pm[0].v[0] = 0 / (Lt/Tt);
        env->pm[0].v[1] = 0 / (Lt/Tt);
        env->pm[0].v[2] = 0 / (Lt/Tt);
        env->pm[0].m      = m / Mt;
    }

    return FFmax;
}

double wall(struct env *env, float M)
{
    int i;

    double em,ff=0,FFmax=0;

    double x,y,z,t,w;

    double Mt, Lt, Tt;
    units(&Mt, &Lt, &Tt);

    assert(env->Nm == 1);

    if (env->Nm) env->pm = malloc(env->Nm * sizeof(*env->pm));
    if (env->Nt) env->pt = malloc(env->Nt * sizeof(*env->pt));

    double v, mu;
    double m = M / env->Nm;
    double r,d;

    for (i=0; i < env->Nt; i++)
    {
        x = (-0.5 - 0.01*(i/10)) / Lt;
        y = (0.2 + 0.01*(i%10)) / Lt;
        z = 0;

        env->pt[i].r[0] = x;
        env->pt[i].r[1] = y;
        env->pt[i].r[2] = z;
        //printf("%f %f %f\n", env->pt[i].r[0], env->pt[i].r[1], env->pt[i].r[2]);

        env->pt[i].v[0] = 0 / (Lt/Tt);
        env->pt[i].v[1] = 0 / (Lt/Tt);
        env->pt[i].v[2] = 0 / (Lt/Tt);

#if 0
        //em = enc_mass(r, M, Rmin/100, Rmax, slope);
        //em = enc_mass(r, M, Rmin, Rmax, slope);
        ff = M_PI/4 * sqrt((2*pow(r / Lt,3))/ (em/Mt));
        printf("%f, %f, %f\n", ff, r, em);
        if (ff > FFmax) 
        {
            //cerr << ff << " " << r << " " << (M) << endl;
            FFmax = ff;
        }
#endif
    }

    if (env->Nm)
    {
        env->pm[0].r[0] = 0;
        env->pm[0].r[1] = 0 / Lt;
        env->pm[0].r[2] = 0 / Lt;
        env->pm[0].v[0] = 0 / (Lt/Tt);
        env->pm[0].v[1] = 0 / (Lt/Tt);
        env->pm[0].v[2] = 0 / (Lt/Tt);
        env->pm[0].m      = m / Mt;
    }

    return FFmax;
}

double hernquist_v(double r)
{
    return 0;
}

double hernquist_r(double Rmin, double Rmax)
{
    double c, r, rho;
    int accept;

    // Normalization constant
    c = Rmin * pow(1 + Rmin,3);

    do
    {
        r = drand48() * (Rmax-Rmin) + Rmin;
        rho = c / r * pow(1 + r, 3);
        accept = drand48() <= rho;
    } while (!accept);

    printf("%f\n", r);

    return r;
}

float sgrstream(struct env *env)
{
    int i;

    double em,ff=0,FFmax=0;

    double x,y,z,t,w;

    double Mt=1, Lt=1, Tt=1;
    //units(&Mt, &Lt, &Tt);

    //assert(env->Nm == 2);

    if (env->Nm) env->pm = malloc(env->Nm * sizeof(*env->pm));
    if (env->Nt) env->pt = malloc(env->Nt * sizeof(*env->pt));

    double v, mu;
    double M = 1;
    double m = M / env->Nm;
    double r,d;
    double Rmin = 0.01;
    double Rmax = 0.1;

    for (i=0; i < env->Nt; i++)
    {
        r = hernquist_r(Rmin, Rmax);
        v = hernquist_v(r);

        z = 2.0 * drand48() - 1.0;

        t = 2.0 * M_PI * drand48();
        w = sqrt(1 - z*z);
        x = w * cos(t);
        y = w * sin(t);


        env->pt[i].r[0] = x * (r / Lt);
        env->pt[i].r[1] = y * (r / Lt);
        env->pt[i].r[2] = z * (r / Lt);
        //printf("%f %f %f\n", env->pt[i].r[0], env->pt[i].r[1], env->pt[i].r[2]);

        env->pt[i].v[0] = 0 / (Lt/Tt);
        env->pt[i].v[1] = 0 / (Lt/Tt);
        env->pt[i].v[2] = 0 / (Lt/Tt);

        //em = enc_mass(r, M, Rmin/100, Rmax, slope);
        //em = enc_mass(r, M, Rmin, Rmax, slope);
        //ff = M_PI/4 * sqrt((2*pow(r / Lt,3))/ (em/Mt));
        //printf("%f, %f, %f\n", ff, r, em);
        //if (ff > FFmax) 
       // {
        //    //cerr << ff << " " << r << " " << (M) << endl;
        //    FFmax = ff;
       // }
    }

    if (env->Nm)
    {
        //x = Rbin/Lt;
        //d = 2*x;
        //v = sqrt(m/Mt * d*x/pow(d*d + env->eps2, 1.5));

        env->pm[0].r[0] = 0 / Lt;
        env->pm[0].r[1] = 0 / Lt;
        env->pm[0].r[2] = 0 / Lt;
        env->pm[0].v[0] = 0 / (Lt/Tt);
        env->pm[0].v[1] = 0;
        env->pm[0].v[2] = 0 / (Lt/Tt);
        env->pm[0].m      = m / Mt;

//      env->pm[1].r[0] = -env->pm[0].r[0];
//      env->pm[1].r[1] =  env->pm[0].r[1];
//      env->pm[1].r[2] =  env->pm[0].r[2];
//      env->pm[1].v[0] =  env->pm[0].v[0];
//      env->pm[1].v[1] = -env->pm[0].v[1];
//      env->pm[1].v[2] =  env->pm[0].v[2];
//      env->pm[1].m    =  env->pm[0].m;
    }

    return FFmax;


}
