#include <math.h>

#include "types.h"
#include "utilities.h"

#ifdef LJMD_MPI
#include"mpi.h"
#endif

void force(mdsys_t *sys) {
    int i, j;
    double rx, ry, rz,rsq;
    double r6, rinv,ffac;

    double c12 = 4.0 * sys->epsilon * pow(sys->sigma, 12.0);
    double c6 = 4.0 * sys->epsilon * pow(sys->sigma, 6.0);
    double rcsq = sys->rcut * sys->rcut;
#ifdef LJMD_MPI
    double epot=0.0;
    azzero(sys->cx,sys->natoms);
    azzero(sys->cy,sys->natoms);
    azzero(sys->cz,sys->natoms);
    MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0,MPI_COMM_WORLD);
#else
    /* zero energy and forces */
    sys->epot = 0.0;
    azzero(sys->fx, sys->natoms);
    azzero(sys->fy, sys->natoms);
    azzero(sys->fz, sys->natoms);
#endif
    for (i = 0; i < (sys->natoms) - 1; i+=sys->nsize) {
        int ii=i+sys->mpirank;
        if (ii >= (sys->natoms - 1)) break;
        for (j = ii + 1; j < (sys->natoms); ++j) {
            rx = pbc(sys->rx[ii] - sys->rx[j], 0.5 * sys->box);
            ry = pbc(sys->ry[ii] - sys->ry[j], 0.5 * sys->box);
            rz = pbc(sys->rz[ii] - sys->rz[j], 0.5 * sys->box);
            rsq = (rx * rx) + (ry * ry) + (rz * rz);
            if (rsq < rcsq) {
                rinv = 1.0 / rsq;
                r6 = rinv * rinv * rinv;
                ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rinv;


#ifdef LJMD_MPI
                epot += r6 * (c12 * r6 - c6);
                sys->cx[ii] += rx * ffac;
                sys->cx[j] -= rx * ffac;
                sys->cy[ii] += ry * ffac;
                sys->cy[j] -= ry * ffac;
                sys->cz[ii] += rz * ffac;
                sys->cz[j] -= rz * ffac;
#else
                sys->epot += r6 * (c12 * r6 - c6);
                sys->fx[i] += rx * ffac;
                sys->fx[j] -= rx * ffac;
                sys->fy[i] += ry * ffac;
                sys->fy[j] -= ry * ffac;
                sys->fz[i] += rz * ffac;
                sys->fz[j] -= rz * ffac;
#endif
            }
        }
    }
#ifdef LJMD_MPI
    MPI_Reduce(sys->cx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(sys->cy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(sys->cz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
}
