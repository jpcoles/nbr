#ifndef NBR_H
#define NBR_H

#include <stdint.h>

struct mparticle
{
    double r[3];       // Position
    double v[3];       // Momentum
    float m;            // Mass
    float eps2;
    float phi;
};

struct tparticle
{
    double r[3];       // Position
    double v[3];       // Momentum
};

struct image
{
    uint32_t nr;
    uint32_t nc;
    unsigned char *img;
};

struct env
{
    int Nt;            // Number of test particles
    int Nm;            // Number of massive particles

    struct tparticle *pt;
    struct mparticle *pm;

    float dt;          // time step
    float t, T;        // current time, total sim time.
    float dt_step;
    long step;        // current step
    float Nsteps;        // current step

    float eps2;

    float extent;

    struct image image;
    char tag[256];

};


#endif