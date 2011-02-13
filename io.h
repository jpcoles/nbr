#ifndef IO_H
#define IO_H

#include <stdint.h>
#include "nbr.h"

void save_snapshot(struct env *env);
int save_image_png(char *fname, unsigned char *image, uint32_t nrows, uint32_t ncols);
void capture(struct env *env);
char *make_fname(struct env *env, char *ext, char *buf, size_t buflen);

#endif
