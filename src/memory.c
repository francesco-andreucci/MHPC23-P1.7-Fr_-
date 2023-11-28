#include"types.h"
#include <stdlib.h>
//#include <sys/time.h>
#include <stdio.h>

 /* allocate memory*/ 
void memalloc(mdsys_t *sys){
    sys->rx=(double *)malloc(sys->natoms*sizeof(double));
    sys->ry=(double *)malloc(sys->natoms*sizeof(double));
    sys->rz=(double *)malloc(sys->natoms*sizeof(double));
    sys->vx=(double *)malloc(sys->natoms*sizeof(double));
    sys->vy=(double *)malloc(sys->natoms*sizeof(double));
    sys->vz=(double *)malloc(sys->natoms*sizeof(double));
    sys->fx=(double *)malloc(sys->natoms*sizeof(double));
    sys->fy=(double *)malloc(sys->natoms*sizeof(double));
    sys->fz=(double *)malloc(sys->natoms*sizeof(double));
#ifdef LJMD_MPI
    sys->cx=(double *)malloc(sys->natoms*sizeof(double));
    sys->cy=(double *)malloc(sys->natoms*sizeof(double));
    sys->cz=(double *)malloc(sys->natoms*sizeof(double));
#endif
}

 /* clean up: close files, free memory */
void cleanup(FILE *erg,FILE *traj, mdsys_t sys){
if(sys.mpirank==0){
    fclose(erg);
    fclose(traj);
}
    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);
#ifdef LJMD_MPI
    free(sys.cx);
    free(sys.cy);
    free(sys.cz);
#endif
}
