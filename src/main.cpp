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


#ifdef _MPI
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
#ifdef _MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &sys.nsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &sys.mpirank);
#endif

    if(sys.mpirank==0){
    printf("LJMD version %3.1f\n", LJMD_VERSION);

    t_start = wallclock();

    /* read input file */
    read_from_file( &sys, &nprint,restfile,trajfile,ergfile,line);
    } //endif myrank=0
    sys.nfi=0;
    //brodcasts
    //reset of rank
    /* allocate memory */
    memalloc(&sys);
if(sys.myrank==0){
    /* read restart */
    read_restfile(restfile, &sys);
}// endif myrank=0

    /* initialize forces and energies.*/

    force(&sys);

    ekin(&sys);

if(sys.myrank==0){
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

        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
            output(&sys, erg, traj);

        /* propagate system and recompute energies */
        velverlet1(&sys);
        /*compute forces and potential energy*/
        force(&sys);
        velverlet2(&sys);
        ekin(&sys);
    }
    /**************************************************/

    /* clean up: close files, free memory */
    printf("Simulation Done. Run time: %10.3fs\n", wallclock()-t_start);


    cleanup(erg,traj,sys);

#ifdef _MPI
    MPI_Finalize();
#endif
    return 0;
}
