#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <png.h>
#include "io.h"

char *make_fname(struct env *env, char *ext, char *buf, size_t buflen)
{
    if (buf == NULL)
    {
	buflen = strlen(env->tag)+1+10+1+strlen(ext)+1;
    	buf = malloc(buflen);
    }
    snprintf(buf, buflen, "%s.%05ld.%s", env->tag, env->step, ext);
    return buf;
}

//==============================================================================
//                               save_snapshot
//==============================================================================
void save_snapshot(struct env *env)
{
    char *fname = make_fname(env, "png", NULL, 0);
    save_image_png(fname, env->image.img, env->image.nr, env->image.nc);
    free(fname);
}

//==============================================================================
//                               save_image_png
//==============================================================================
int save_image_png(char *fname, unsigned char *image, uint32_t nrows, uint32_t ncols)
{
    int i;
    FILE *fp = fopen(fname, "wb");

    if (fp == NULL)
    {
        fprintf(stderr, "Can't open %s\n", fname);
        return 1;
    }

    png_structp png_ptr = png_create_write_struct
       (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) return 1;

    png_init_io(png_ptr, fp);

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
    {
       png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
       return 1;
    }

    png_set_IHDR(png_ptr, info_ptr, ncols, nrows,
           8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
           PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    png_write_info(png_ptr, info_ptr);

    int row_stride = ncols * 3;

    for (i=0; i < nrows; i++)
    {
        png_bytep row_pointer = & (image[i * row_stride]);
        png_write_row(png_ptr, row_pointer);
    }

    png_write_end(png_ptr, info_ptr);

    png_destroy_write_struct(&png_ptr, &info_ptr);

    fclose(fp);

    return 0;
}

//==============================================================================
//                                  capture
//==============================================================================
void capture(struct env *env)
{
    struct mparticle *p = env->pm;
    struct tparticle *t = env->pt;

    //memset(env->image.img, 0, 3*env->image.nc*env->image.nr*sizeof(*(env->image.img)));

    uint32_t nr = env->image.nr;
    uint32_t nc = env->image.nc;
    size_t i;
    for (i=0; i < env->Nm; i++)
    {
        int32_t c = ( (p[i].r[0]-env->CoM[0]) + env->extent) / (2.0*env->extent) * nc;
        int32_t r = (-(p[i].r[1]-env->CoM[1]) + env->extent) / (2.0*env->extent) * nr;
        if (!(0 <= r && r < nr)) continue;
        if (!(0 <= c && c < nc)) continue;
        env->image.img[3*(r*nc + c) + 0] = (255);
        env->image.img[3*(r*nc + c) + 1] = (0);
        env->image.img[3*(r*nc + c) + 2] = (0);
    }

    for (i=0; i < env->Nt; i++)
    {
        int32_t c = ( t[i].r[0] + env->extent) / (2.0*env->extent) * nc;
        int32_t r = (-t[i].r[1] + env->extent) / (2.0*env->extent) * nr;
        if (!(0 <= r && r < nr)) continue;
        if (!(0 <= c && c < nc)) continue;
        env->image.img[3*(r*nc + c) + 0] = (255);
        env->image.img[3*(r*nc + c) + 1] = (255);
        env->image.img[3*(r*nc + c) + 2] = (255);
    }
}
