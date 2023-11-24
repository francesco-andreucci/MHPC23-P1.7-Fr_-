#ifndef MEMORY_H
#ifdef __cplusplus
extern "C" {
#endif
#define MEMORY_H

#include"types.h"

void cleanup(FILE *erg,FILE *traj, mdsys_t sys);
void allocate(mdsys_t sys);

#ifdef __cplusplus
}
#endif

#endif