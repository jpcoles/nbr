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

void units_solsystem(double *M, double *L, double *T)
{
    double days = 24 * 60 * 60;               // [s]
    double AU   = 1.49597870691e11;           // [m]
    double Msun = 1.98892e30;                 // [kg]

    double Gsi = 6.67300e-11*pow(AU,-3)*Msun*pow(days,2); // [kpc^3 msun^-1 gyr^-2]

    //*L = 1.0;   // [AU]
    *L = 1e-1;   // [AU]
    *T = 1e3;   // [days] / sqrt(Gsi * pow(*L,-3) * *M);
    *M = pow(*T,-2) * pow(Gsi,-1) * pow(*L,3); // [Msun]

    fprintf(stderr, "G=%e\n", *M * pow(*T,2) * pow(Gsi,1) * pow(*L,-3));
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

float solsystem(struct env *env)
{
    int i;
    double Mt=1, Lt=1, Tt=1;
    units_solsystem(&Mt, &Lt, &Tt);
    env->units.M = Mt;
    env->units.L = Lt;
    env->units.T = Tt;
    env->units.G = 1;

    double Msun = 1.98892e30;                 // [kg]
    double Mjup = 1.8987e27;

    struct {
        double x,y,z;
        double vx,vy,vz;
        double M;
    } SS[] = 
    {
     // Sun
     { 0.,0.,0., 0., 0., 0., 1},

     // Mercury
     { 1.370754491185286E-01, -4.264661186200457E-01, -4.741366654969732E-02, 
       2.114746190733265E-02,  1.004581022861040E-02, -1.122428604918364E-03,   3.302e23 / Msun},
     // Venus
     { 7.229565616872776E-01,  5.608001294512701E-02, -4.098550983391910E-02, 
      -1.646858780384053E-03,  2.007293064458840E-02,  3.677906364857320E-04,  48.685e23 / Msun},
     // Earth
     {-9.472951247336898E-01,  2.926531282918977E-01,  1.678530616122745E-05, 
      -5.365105883437778E-03, -1.650359033144609E-02, -1.504100835864763E-06, 5.9736e24 / Msun},
     // Moon
     {-9.467458988729617E-01,  2.902726107565438E-01, -1.779573682223667E-04, 
      -4.769970905362414E-03, -1.634099625429580E-02,  2.968006365406434E-05, 734.9e20 / Msun},
     // Mars
     { 8.323613394797202E-01,  1.231253451781860E+00,  5.268023383023045E-03, 
      -1.105805890906632E-02,  9.028381206953562E-03,  4.615746913514457E-04, 6.4185e23 / Msun},
     // Jupiter
     { //-5.287494142904107E+00, -1.243507235129249E+00, -3.932468568369884E-01, 
       //1.687413721789376E-03, -6.977384452214245E-03,  5.658301989168687E-04, 1.89813e27 / Msun},
       -4.817300396248593E+00, -2.537352081929040E+00,  1.183594107357213E-01, 
       3.430544476066865E-03, -6.330892933515803E-03, -5.081277966460052E-05, 1.89813e27 / Msun},
     // Saturn
     { 7.002859392650318E+00,  5.953410760786836E+00, -3.825476620350814E-01, 
      -3.912528910160785E-03,  4.233888993422562E-03,  8.166441689583057E-05, 5.68319e26 / Msun},
     // Uranus
     {-1.818744549523943E+01, -2.186144292636072E+00,  2.279359559634153E-01, 
       4.420435097205414E-04, -4.093741565202423E-03, -2.092784146624804E-05, 86.8103e23 / Msun},
     // Neptune
     {-1.539752139377567E+01, -2.610796272562613E+01,  8.922320120651117E-01, 
       2.686374686169427E-03, -1.581855900310080E-03, -2.968604130174336E-05, 1.0241e26 / Msun},
     // Pluto
     {-3.039874622789153E+01,  1.918475637076636E+00,  8.585829777861537E+00, 
       3.999928092384586E-04, -3.325479950192453E-03,  2.292679556968999E-04, 1.314e22 / Msun},

     // X
     { 0.0, 200.0,0.0,
       0.0,0.0,0.0,     25*Mjup / Msun},
     // Y
     { 0.0, -200.0,0.0,
       0.0,0.0,0.0,     25*Mjup / Msun},
    };

    //assert(env->Nm == 10);

    if (env->Nm) env->pm = malloc(env->Nm * sizeof(*env->pm));
    if (env->Nt) env->pt = malloc(env->Nt * sizeof(*env->pt));

    double M = 0;

    env->pm[0].r[0] = SS[0].x / Lt;
    env->pm[0].r[1] = SS[0].y / Lt;
    env->pm[0].r[2] = SS[0].z / Lt;
    env->pm[0].v[0] = SS[0].vx / (Lt/Tt);
    env->pm[0].v[1] = SS[0].vy / (Lt/Tt);
    env->pm[0].v[2] = SS[0].vz / (Lt/Tt);
    env->pm[0].m    = SS[0].M / Mt;

    fprintf(stderr, "M0 is %e\n", env->pm[0].m);

    M = env->pm[0].m;

    double cx=0, cy=0, cz=0;
    cx = env->pm[0].m*env->pm[0].r[0];
    cy = env->pm[0].m*env->pm[0].r[1];
    cz = env->pm[0].m*env->pm[0].r[2];

    for (i=1; i < env->Nm; i++)
    {
        env->pm[i].r[0] = SS[i].x / Lt;
        env->pm[i].r[1] = SS[i].y / Lt;
        env->pm[i].r[2] = SS[i].z / Lt;
        env->pm[i].v[0] = SS[i].vx / (Lt/Tt);
        env->pm[i].v[1] = SS[i].vy / (Lt/Tt);
        env->pm[i].v[2] = SS[i].vz / (Lt/Tt);

#if 1
        if (i==(env->Nm-1) || i==(env->Nm-2))
        {
            double vx, vy, vz;
            double r = sqrt(pow(env->pm[i].r[0],2)
                          + pow(env->pm[i].r[1],2)
                          + pow(env->pm[i].r[2],2));
            double theta = atan2(env->pm[i].r[1], env->pm[i].r[0]);

            vx = 0;
            vy = +sqrt(M/fabs(r));
            vz = 0;

            env->pm[i].v[0] = vx*cos(theta) - vy*sin(theta);
            env->pm[i].v[1] = vx*sin(theta) + vy*cos(theta);
            env->pm[i].v[2] = 0;
        }
#endif


        double r = sqrt(pow(env->pm[i].r[0],2)
                      + pow(env->pm[i].r[1],2)
                      + pow(env->pm[i].r[2],2));

        env->pm[i].m    = SS[i].M / Mt;
        fprintf(stderr, "M%i is %e at R=%e\n", i, env->pm[i].m, r);
        if (i<10)
        {
            M += env->pm[i].m;
        }

        cx += env->pm[i].m*env->pm[i].r[0];
        cy += env->pm[i].m*env->pm[i].r[1];
        cz += env->pm[i].m*env->pm[i].r[2];
    }

    cx /= M;
    cy /= M;
    cz /= M;

    //env->pm[10].v[1] = +sqrt(M/fabs(env->pm[10].r[0]));

    fprintf(stderr, "COM %e %e %e\n", cx, cy, cz);

    //env->pm[10].v[1] = +sqrt(M/fabs(env->pm[10].r[0]));
    //env->pm[11].v[1] = -sqrt(M/fabs(env->pm[11].r[0]));
        
    
    return 0;
}
