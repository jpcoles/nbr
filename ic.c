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

#if 0
double ball(prop_t *prop, TipsyDarkParticle *p, int p_len, int stride)
{
    int    N       = prop->N;
    double slope   = prop->slope;
    double M       = prop->M;
    double Rmax    = prop->Rmax;
    double Rmin    = prop->Rmin;
    //double Rbinary = prop->Rbinary;
    double soft    = prop->soft;

    int32_t i=0;
    double x=0,y=0,z=0,r,t,w;
    double x0=0,y0=0,z0=0;
    double FFmax = 0; /* Max free fall time */

    double Mt, Lt, Tt;
    units(&Mt, &Lt, &Tt);

    double m = M / N;


    cerr << "Tt=" << Tt << endl;


    //double Tmin = HUGE_VAL;


    if (prop->mode & MODE_DROPLET)
    {
        z0 = 2.0 * drand48() - 1.0;
        t  = 2.0 * M_PI * drand48();
        w  = sqrt(1 - z*z);
        x0 = w * cos(t);
        y0 = w * sin(t);
    }

    //for (i=0; i < N-2; i++)
    for (i=0; i < N; i+=stride)
    {
        if (prop->mode & MODE_RMAX_SHELL) 
            r = Rmax;
        else
            r = new_point(Rmin, Rmax, slope);

        if (prop->mode & MODE_DROPLET)
        {
            x = x0 + prop->Rdroplet*(2*drand48()-1) * Lt;
            y = y0 + prop->Rdroplet*(2*drand48()-1) * Lt;
            z = z0 + prop->Rdroplet*(2*drand48()-1) * Lt;
        }
        else
        {
            if (prop->mode & MODE_2D)
                z = 0;
            else
                z = 2.0 * drand48() - 1.0;

            t = 2.0 * M_PI * drand48();
            w = sqrt(1 - z*z);
            x = w * cos(t);
            y = w * sin(t);

            if (prop->mode & MODE_YZ)
            {
                z = x;
                x = 0;
            }

            if (prop->mode & MODE_XZ)
            {
                z = y;
                y = 0;
            }

        }

        //double m = enc_mass(r, M, Rmin, Rmax, slope);
        //m /= 2; // Each particle gets half the total mass.

        //m = M / 2;

        p[i].pos[0] = x * (r / Lt);
        p[i].pos[1] = y * (r / Lt);
        p[i].pos[2] = z * (r / Lt);

        p[i].vel[0] = 0 / (Lt/Tt);
        p[i].vel[1] = 0 / (Lt/Tt);
        p[i].vel[2] = 0 / (Lt/Tt);
        p[i].eps    = soft / Lt;
        p[i].mass   = m / Mt;
        //cerr << p[i+0].mass << endl;
        p[i].phi    = 0;

#if 0
        x = Rbinary;
        double mu = (m/Mt)/2; //pow(M/Mt,2) / (2 * M/Mt);
        double v  = sqrt(mu / fabs(2*x/Lt));

        if (prop->mode & MODE_STATIC_BINARY) v = 0;

        cerr << v << endl;

        p[i+1].pos[0] = x / Lt;
        p[i+1].pos[1] = 0 / Lt;
        p[i+1].pos[2] = 0 / Lt;
        p[i+1].vel[0] = 0 / (Lt/Tt);
        p[i+1].vel[1] = v;
        p[i+1].vel[2] = 0 / (Lt/Tt);
        p[i+1].eps    = soft / Lt;
        p[i+1].mass   = m / Mt;
        p[i+1].phi    = 0;

        p[i+2].pos[0] = -x / Lt;
        p[i+2].pos[1] = 0 / Lt;
        p[i+2].pos[2] = 0 / Lt;
        p[i+2].vel[0] = 0 / (Lt/Tt);
        p[i+2].vel[1] = -v;
        p[i+2].vel[2] = 0 / (Lt/Tt);
        p[i+2].eps    = soft / Lt;
        p[i+2].mass   = m / Mt;
        p[i+2].phi    = 0;
#endif

        double em = enc_mass(r, M, Rmin, Rmax, slope);
        double ff = M_PI/4 * sqrt((2*pow(r / Lt,3))/ (em/Mt));
        if (ff > FFmax) 
        {
            cerr << ff << " " << r << " " << (M) << endl;
            FFmax = ff;
        }

#if 0
        double tt = 2*M_PI*(Rbinary/Lt) / v;
        if (tt < Tmin) 
        {
            Tmin = tt;
            cerr << tt << " " << v << endl;
        }
#endif
    }

    return FFmax;


    //cerr << "Max crossing time is ~" << (2*sqrt(3/M_PI /* / G=1 */ * maxR/kpc)) << " Gyr" << endl;
    //cerr << "twiddle crossing time is ~" << (4./3.*sqrt(N/(M/Msun)) * pow(maxR/kpc,1.5)) << " Gyr" << endl;

    //return M_PI/4 * sqrt((2*pow(r * R / Lt,3))/ (3*M/Mt));
}
#endif

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
        printf("%f, %f, %f\n", ff, r, em);
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
        printf("%f, %f, %f\n", ff, r, em);
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

#if 0

void help(char *track)
{
    if (track != NULL)
        cerr << track << endl;

    cerr << "Usage: tbr OPTIONS <N> <mass> <Rmax> <Rmin> <slope> <softening> <out file>" << endl
         << endl
         << "where OPTIONS are:" << endl
         << "   --seed=N            Initialize random seed to integer N." << endl
         << "                       Default is to base the seed on the time(NULL)." << endl
         << "   --mode=[0,1]        Select a configuration:" << endl
         << "                           0: Normal, binary orbiting pair and tracer particle." << endl
         << "                           1: Binary pair is not orbiting." << endl
         << "                           2: All particles are at Rmax." << endl
         << "                           3: All particles are in a small box." << endl
         << endl;

    exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
    ofTipsy out;
    TipsyHeader       h;

    prop_t prop;

    prop.mode = MODE_NORMAL;
    prop.Rdroplet = 0;
    prop.Rbinary = 0;
    prop.seed = time(NULL);

    static struct option long_options[] = {
        {"seed", required_argument, 0, 0},
        {"mode", required_argument, 0, 0},
        {"drad", required_argument, 0, 0},
        {"2d", no_argument, 0, 0},
        {"xy", no_argument, 0, 0},
        {"yz", no_argument, 0, 0},
        {"xz", no_argument, 0, 0},
        {"Rbin", required_argument, 0, 0},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    /*--------------------------------------------------------------------------
     * Process the command line flags
     *------------------------------------------------------------------------*/
    while (1)
    {
        int option_index = 0;
        int c = getopt_long(argc, argv, "h",
                            long_options, &option_index);

        if (c == -1) break;

        switch (c)
        {
            case 0:
                if (!strcmp("seed", long_options[option_index].name))
                {
                    prop.seed = atoi(optarg);
                }
                if (!strcmp("mode", long_options[option_index].name))
                {
                    prop.mode |= atoi(optarg);
                    if (prop.mode & ~ALL_MODES) help("1");
                }
                if (!strcmp("drad", long_options[option_index].name))
                {
                    prop.Rdroplet = atof(optarg);
                }
                if (!strcmp("2d", long_options[option_index].name))
                {
                    prop.mode |= MODE_2D;
                }
                if (!strcmp("yz", long_options[option_index].name))
                    prop.mode |= MODE_YZ;
                if (!strcmp("xz", long_options[option_index].name))
                    prop.mode |= MODE_XZ;
                if (!strcmp("xy", long_options[option_index].name))
                    prop.mode |= MODE_XY;

                if (!strcmp("Rbin", long_options[option_index].name))
                    prop.Rbinary = atof(optarg);

                break;

            case 'h': /* Fall through */
            case '?': /* Fall through */
            default:
                help("2");
        }
    }

    srand48(prop.seed);

    if (argc-optind < 7) help("3");

    // -- Change these values --
    //double R    = 20;
    //double M    = 1e11;
    //double soft = 1;
    // -------------------------

    prop.N       = atoi(argv[optind+0]);   // N
    prop.M       = atof(argv[optind+1]);   // M
    prop.Rmax    = atof(argv[optind+2]);   // Rmax
    prop.Rmin    = atof(argv[optind+3]);   // Rmin
    prop.slope   = atof(argv[optind+4]);   // slope
    prop.soft    = atof(argv[optind+5]);   // soft

    string mark(argv[optind+6]);


    assert(prop.slope >= 0);
    prop.slope *= -1;

    TipsyDarkParticle *d;

    double ff;

    if (prop.mode & (MODE_BINARY | MODE_STATIC_BINARY))
    {
        if (prop.Rbinary == 0) help("4");
        ff = tbr(&prop, &d);
    }
    else
    {
        d = new TipsyDarkParticle[prop.N];
        ff = ball(&prop, d,prop.N,1);
    }

    h.h_time = 0;
    h.h_nBodies = prop.N;
    h.h_nDims = 3;
    h.h_nDark = prop.N;
    h.h_nStar = 0;
    h.h_nSph = 0;

    out.open(mark.c_str(), "standard");
    mark += ".mark";

    cerr << mark << endl;
    ofstream markout(mark.c_str(), ios::out);

    out << h;
    markout << " " << h.h_nDark << " " << h.h_nSph<< " " << h.h_nStar << endl;
    for (int32_t i=0; i < prop.N; i++) 
    {
        out << d[i];
        if (d[i].mass < 1) markout << (i+1) << endl;
    }
    markout.close();
    out.close();

    delete d;

    cout << setprecision(10) << ff << endl;
    cerr << setprecision(10) << ff << endl;

    cerr << "N="        << prop.N << endl
         << "M="        << prop.M << endl
         << "Rmax="     << prop.Rmax << endl
         << "Rmin="     << prop.Rmin << endl
         << "Rbinary="  << prop.Rbinary << endl
         << "slope="    << prop.slope << endl
         << "soft="     << prop.soft << endl
         << "seed="     << prop.seed << endl
         << "mode="     << prop.mode << endl
         << "mark="     << mark << endl
         << "Tff="      << ff << endl;

    return EXIT_SUCCESS;
}


void help(char *track)
{
    if (track != NULL)
        cerr << track << endl;

    cerr << "Usage: tbr OPTIONS <N> <mass> <Rmax> <Rmin> <slope> <softening> <out file>" << endl
         << endl
         << "where OPTIONS are:" << endl
         << "   --seed=N            Initialize random seed to integer N." << endl
         << "                       Default is to base the seed on the time(NULL)." << endl
         << "   --mode=[0,1]        Select a configuration:" << endl
         << "                           0: Normal, binary orbiting pair and tracer particle." << endl
         << "                           1: Binary pair is not orbiting." << endl
         << "                           2: All particles are at Rmax." << endl
         << "                           3: All particles are in a small box." << endl
         << endl;

    exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
    ofTipsy out;
    TipsyHeader       h;

    prop_t prop;

    prop.mode = MODE_NORMAL;
    prop.Rdroplet = 0;
    prop.Rbinary = 0;
    prop.seed = time(NULL);

    static struct option long_options[] = {
        {"seed", required_argument, 0, 0},
        {"mode", required_argument, 0, 0},
        {"drad", required_argument, 0, 0},
        {"2d", no_argument, 0, 0},
        {"xy", no_argument, 0, 0},
        {"yz", no_argument, 0, 0},
        {"xz", no_argument, 0, 0},
        {"Rbin", required_argument, 0, 0},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    /*--------------------------------------------------------------------------
     * Process the command line flags
     *------------------------------------------------------------------------*/
    while (1)
    {
        int option_index = 0;
        int c = getopt_long(argc, argv, "h",
                            long_options, &option_index);

        if (c == -1) break;

        switch (c)
        {
            case 0:
                if (!strcmp("seed", long_options[option_index].name))
                {
                    prop.seed = atoi(optarg);
                }
                if (!strcmp("mode", long_options[option_index].name))
                {
                    prop.mode |= atoi(optarg);
                    if (prop.mode & ~ALL_MODES) help("1");
                }
                if (!strcmp("drad", long_options[option_index].name))
                {
                    prop.Rdroplet = atof(optarg);
                }
                if (!strcmp("2d", long_options[option_index].name))
                {
                    prop.mode |= MODE_2D;
                }
                if (!strcmp("yz", long_options[option_index].name))
                    prop.mode |= MODE_YZ;
                if (!strcmp("xz", long_options[option_index].name))
                    prop.mode |= MODE_XZ;
                if (!strcmp("xy", long_options[option_index].name))
                    prop.mode |= MODE_XY;

                if (!strcmp("Rbin", long_options[option_index].name))
                    prop.Rbinary = atof(optarg);

                break;

            case 'h': /* Fall through */
            case '?': /* Fall through */
            default:
                help("2");
        }
    }

    srand48(prop.seed);

    if (argc-optind < 7) help("3");

    // -- Change these values --
    //double R    = 20;
    //double M    = 1e11;
    //double soft = 1;
    // -------------------------

    prop.N       = atoi(argv[optind+0]);   // N
    prop.M       = atof(argv[optind+1]);   // M
    prop.Rmax    = atof(argv[optind+2]);   // Rmax
    prop.Rmin    = atof(argv[optind+3]);   // Rmin
    prop.slope   = atof(argv[optind+4]);   // slope
    prop.soft    = atof(argv[optind+5]);   // soft

    string mark(argv[optind+6]);


    assert(prop.slope >= 0);
    prop.slope *= -1;

    TipsyDarkParticle *d;

    double ff;

    if (prop.mode & (MODE_BINARY | MODE_STATIC_BINARY))
    {
        if (prop.Rbinary == 0) help("4");
        ff = tbr(&prop, &d);
    }
    else
    {
        d = new TipsyDarkParticle[prop.N];
        ff = ball(&prop, d,prop.N,1);
    }

    h.h_time = 0;
    h.h_nBodies = prop.N;
    h.h_nDims = 3;
    h.h_nDark = prop.N;
    h.h_nStar = 0;
    h.h_nSph = 0;

    out.open(mark.c_str(), "standard");
    mark += ".mark";

    cerr << mark << endl;
    ofstream markout(mark.c_str(), ios::out);

    out << h;
    markout << " " << h.h_nDark << " " << h.h_nSph<< " " << h.h_nStar << endl;
    for (int32_t i=0; i < prop.N; i++) 
    {
        out << d[i];
        if (d[i].mass < 1) markout << (i+1) << endl;
    }
    markout.close();
    out.close();

    delete d;

    cout << setprecision(10) << ff << endl;
    cerr << setprecision(10) << ff << endl;

    cerr << "N="        << prop.N << endl
         << "M="        << prop.M << endl
         << "Rmax="     << prop.Rmax << endl
         << "Rmin="     << prop.Rmin << endl
         << "Rbinary="  << prop.Rbinary << endl
         << "slope="    << prop.slope << endl
         << "soft="     << prop.soft << endl
         << "seed="     << prop.seed << endl
         << "mode="     << prop.mode << endl
         << "mark="     << mark << endl
         << "Tff="      << ff << endl;

    return EXIT_SUCCESS;
}

#endif
