#ifndef INPUT_H
#ifdef __cplusplus
extern "C" {
#endif
#define INPUT_H

#include<stdio.h>
#include"types.h"
/* generic file- or pathname buffer length */
#define BLEN 200
int get_a_line(FILE *fp, char *buf);

void read_from_file(mdsys_t * sys, int * nprint,char * restfile, char* trajfile, char * ergfile,char * line);

void read_restfile(char *restfile, mdsys_t *sys);

#ifdef __cplusplus
}
#endif

#endif
