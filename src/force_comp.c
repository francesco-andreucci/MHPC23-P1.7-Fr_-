#include <math.h>
#include "types.h"
#include "utilities.h"

#ifdef LJMD_MPI
#include"mpi.h"
#endif

#ifdef _OPENMP
    #include "omp.h"
#endif

void force(mdsys_t *sys) {
    
    int tid, start, end;
    int i, j;
    double rx, ry, rz, rsq;
    double r6, rinv, ffac;

    double c12 = 4.0 * sys->epsilon * pow(sys->sigma, 12.0);
    double c6 = 4.0 * sys->epsilon * pow(sys->sigma, 6.0);
    double rcsq = sys->rcut * sys->rcut;

    double epot = 0.0;

#ifdef LJMD_MPI
    //double epot=0.0;
    MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0,MPI_COMM_WORLD);
#endif

#ifdef _OPENMP
    #pragma omp parallel reduction(+:epot) firstprivate(c12, c6, rcsq, i, j, tid, start, end, rx, ry, rz, rsq, r6, rinv, ffac)
#endif
{ 
    double *fx, *fy, *fz;
    // double *cx, *cy, *cz;
  
    /* reading thread id and defining default tid when omp is absent */
    #ifdef _OPENMP
        tid = omp_get_thread_num();
    #else
        tid = 0;
    #endif

    fx = sys->fx + (tid * sys->natoms);
    fy = sys->fy + (tid * sys->natoms);
    fz = sys->fz + (tid * sys->natoms);

    #ifdef LJMD_MPI
    azzero(sys->cx,sys->natoms);
    azzero(sys->cy,sys->natoms);
    azzero(sys->cz,sys->natoms);
    #else
    azzero(fx, sys->natoms);
    azzero(fy, sys->natoms);
    azzero(fz, sys->natoms);
    #endif

        for (i = 0; i < (sys->natoms - 1); i += sys->tmax * sys->nsize)
        {
            int ii = i + tid + sys->tmax * sys->mpirank;
            if(ii >= (sys->natoms - 1)) break;

            for (j = ii + 1; j < sys->natoms; ++j)
            {
                rx = pbc(sys->rx[ii] - sys->rx[j], 0.5 * sys->box);
                ry = pbc(sys->ry[ii] - sys->ry[j], 0.5 * sys->box);
                rz = pbc(sys->rz[ii] - sys->rz[j], 0.5 * sys->box);
                rsq = (rx * rx) + (ry * ry) + (rz * rz);
                    
                if (rsq < rcsq)
                {
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
                    epot += r6 * (c12 * r6 - c6);
                    fx[ii] += rx * ffac;
                    fx[j] -= rx * ffac;
                    fy[ii] += ry * ffac;
                    fy[j] -= ry * ffac;
                    fz[ii] += rz * ffac;
                    fz[j] -= rz * ffac;
                #endif
                }
            }
        }
        
        #ifdef _OPENMP
            #pragma omp barrier
        #endif
        
        int i = 1 + sys->natoms / sys->tmax;
        start = tid * i;
        end = start + i;

        if(end > sys->natoms) 
        {
            end = sys->natoms; 
        }

        for(i = 1; i < sys->tmax; ++i)
        {   
            int offset = i * sys->natoms;
            for(int j = start; j < end; ++j)
            {   
                #ifdef LJMD_MPI
                sys->cx[j] += sys->cx[offset + j];
                sys->cy[j] += sys->cy[offset + j];
                sys->cz[j] += sys->cz[offset + j];
                #else
                sys->fx[j] += sys->fx[offset + j];
                sys->fy[j] += sys->fy[offset + j];
                sys->fz[j] += sys->fz[offset + j];
                #endif
            }
        }
    }
    sys->epot = epot;

#ifdef LJMD_MPI
    MPI_Reduce(sys->cx, sys->fx, sys->natoms*sys->tmax, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(sys->cy, sys->fy, sys->natoms*sys->tmax, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(sys->cz, sys->fz, sys->natoms*sys->tmax, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
}