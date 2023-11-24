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


#define LJMD_VERSION 0.1



int main(int argc, char **argv){
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;
    mdsys_t sys;
    double t_start;

    printf("LJMD version %3.1f\n", LJMD_VERSION);

    t_start = wallclock();

    /* read input file */
    read_from_file( &sys, &nprint,restfile,trajfile,ergfile,line);


    /* allocate memory */
    allocate(sys);

    /* read restart */

    read_restfile(restfile, &sys);
    /* initialize forces and energies.*/
    sys.nfi=0;
    force(&sys);
    ekin(&sys);

    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");

    printf("Startup time: %10.3fs\n", wallclock()-t_start);
    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys, erg, traj);

    /* reset timer */
    t_start = wallclock();

    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
            output(&sys, erg, traj);

        /* propagate system and recompute energies */
        velverlet1(&sys);
        velverlet2(&sys);
        ekin(&sys);
    }
    /**************************************************/

    /* clean up: close files, free memory */
    printf("Simulation Done. Run time: %10.3fs\n", wallclock()-t_start);

    cleanup(erg,traj,sys);

    return 0;
}
