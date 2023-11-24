#ifndef UTILITIES_H
#ifdef __cplusplus
extern "C" {
#endif
#define UTILITIES_H

double wallclock();

double pbc(double x, const double boxby2);

void ekin(mdsys_t *sys);

void azzero(double *d, const int n);

#ifdef __cplusplus
}
#endif


#endif