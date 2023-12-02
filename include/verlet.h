#ifndef VERLET_H
#ifdef __cplusplus
extern "C" {
#endif
#define VERLET_H

#include"types.h"

void velverlet1(mdsys_t *sys);

void velverlet2(mdsys_t *sys);

#ifdef __cplusplus
}
#endif

#endif 
