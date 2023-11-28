#include"types.h"
#include"utilities.h"
#include"input.h"
#include"output.h"
#include"verlet.h"
#include"memory.h"
#include"force_comp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef LJMD_MPI
#include "mpi.h"
#endif

#define LJMD_VERSION 0.1



int main(int argc, char **argv){
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;
    mdsys_t sys;
    double t_start;
    sys.nsize=1;
    sys.mpirank=0;
#ifdef LJMD_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &sys.nsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &sys.mpirank);
#endif

    if(sys.mpirank==0){
    printf("LJMD version %3.1f\n", LJMD_VERSION);

    t_start = wallclock();

    /* read input file */
    read_from_file( &sys, &nprint,restfile,trajfile,ergfile,line);
    sys.nfi=0;
    } //endif myrank=0
#ifdef LJMD_MPI
    MPI_Bcast(&sys.natoms, 1, MPI_INT, 0,MPI_COMM_WORLD);
    MPI_Bcast(&sys.mass, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(&sys.epsilon, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(&sys.sigma, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(&sys.rcut, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(&sys.box, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(&sys.nsteps, 1, MPI_INT, 0,MPI_COMM_WORLD);
    MPI_Bcast(&sys.dt, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(&sys.nfi, 1, MPI_INT, 0,MPI_COMM_WORLD);
#endif

    /* allocate memory */
    memalloc(&sys);
    /*
    sys.cx=(double *)malloc(sys.natoms*sizeof(double));
    sys.cy=(double *)malloc(sys.natoms*sizeof(double));
    sys.cz=(double *)malloc(sys.natoms*sizeof(double));
*/
    /*only rank 0 reads input*/
if(sys.mpirank==0){
    /* read restart */
    read_restfile(restfile, &sys);
}// endif myrank=0
#ifdef LJMD_MPI
  //  MPI_Bcast(sys.rx, sys.natoms, MPI_DOUBLE, 0,MPI_COMM_WORLD);
  //  MPI_Bcast(sys.ry, sys.natoms, MPI_DOUBLE, 0,MPI_COMM_WORLD);
  //  MPI_Bcast(sys.rz, sys.natoms, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(sys.vx, sys.natoms, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(sys.vy, sys.natoms, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(sys.vz, sys.natoms, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(sys.fx, sys.natoms, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(sys.fy, sys.natoms, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(sys.fz, sys.natoms, MPI_DOUBLE, 0,MPI_COMM_WORLD);
#endif
    /* initialize forces and energies.*/

    force(&sys);

    ekin(&sys);

if(sys.mpirank==0){
    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");

    printf("Startup time: %10.3fs\n", wallclock()-t_start);
    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys, erg, traj);

    /* reset timer */
    t_start = wallclock();
}


    /**************************************************/
    /* main MD loop */

    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {
        if(sys.mpirank==0){
        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
            output(&sys, erg, traj);
        }
        /* propagate system and recompute energies */
        velverlet1(&sys);
        /*compute forces and potential energy*/
        force(&sys);
        velverlet2(&sys);
        ekin(&sys);
    }
    /**************************************************/
if(sys.mpirank==0){
    /* clean up: close files, free memory */
    printf("Simulation Done. Run time: %10.3fs\n", wallclock()-t_start);
}
if(sys.mpirank==0){
    cleanup(erg,traj,sys);
}else{
    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);
    free(sys.cx);
    free(sys.cy);
    free(sys.cz);

}
#ifdef LJMD_MPI
    MPI_Finalize();
#endif
    return 0;
}
