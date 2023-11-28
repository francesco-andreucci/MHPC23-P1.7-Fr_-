#include <math.h>
#include "types.h"
#include "utilities.h"
#include <stdio.h>
#ifdef _OPENMP
    #include "omp.h"
#endif

void force(mdsys_t *sys) {
    
    
    double epot = 0.0;

    #ifdef _OPENMP
        sys->tmax = omp_get_max_threads();
    #endif

    #ifdef _OPENMP
        #pragma omp parallel reduction(+:epot) //private(i, j, tid, rsq, ffac, r6, rinv) 
    #endif
    { 
    int i, j, tid, start, end, offset;
    double *fx, *fy, *fz;
    double c12,c6,rcsq;


    /* zero energy and forces */
    epot = 0.0;
    c12 = 4.0 * sys->epsilon * pow(sys->sigma, 12.0);
    c6 = 4.0 * sys->epsilon * pow(sys->sigma, 6.0);
    rcsq = sys->rcut * sys->rcut;

    // double r, rsq, ffac;
    // double rx, ry, rz;
    // int start, end;
  
    #ifdef _OPENMP
        sys->tmax = omp_get_num_threads();
        tid = omp_get_thread_num();
    #else
        tid = 0;
    #endif

            // double *fx, *fy, *fz;
            // int start, end;

            fx = sys->fx + (tid * sys->natoms);
            fy = sys->fy + (tid * sys->natoms);
            fz = sys->fz + (tid * sys->natoms);

            azzero(fx, sys->tmax*sys->natoms);
            azzero(fy, sys->tmax*sys->natoms);
            azzero(fz, sys->tmax*sys->natoms);

            for (int i = 0; i < (sys->natoms - 1); i += sys->tmax)
            {
                int ii, j;
                double rx1, ry1, rz1;
                ii = i + tid;

                //printf("sys->tmax: %d\t | i: %d\t | tid: %d | ii: %d", sys->tmax, i, tid, ii);
                // printf("sys->tmax: %d\n", sys->tmax);
                if (ii >= (sys->natoms - 1))
                    break;

                rx1 = sys->rx[ii];
                ry1 = sys->ry[ii];
                rz1 = sys->rz[ii];

                for (int j = ii + 1; j < sys->natoms; ++j)
                {
                    double rx, ry, rz, rsq;
                    rx = pbc(rx1 - sys->rx[j], 0.5 * sys->box);
                    ry = pbc(ry1 - sys->ry[j], 0.5 * sys->box);
                    rz = pbc(rz1 - sys->rz[j], 0.5 * sys->box);
                    rsq = (rx * rx) + (ry * ry) + (rz * rz);
                    
                    //printf("tid: %d, j: %d, ii: %d, rsq: %lf, rcsq: %lf, rsq < rcsq: %d\n", tid, j, ii, rsq, rcsq, rsq < rcsq);

                    if (rsq < rcsq)
                    {
                        double r6, rinv, ffac;
                        rinv = 1.0 / rsq;
                        r6 = rinv * rinv * rinv;
                        ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rinv;
                        epot += r6 * (c12 * r6 - c6);
                        fx[ii] += rx * ffac;
                        fx[j] -= rx * ffac;
                        fy[ii] += ry * ffac;
                        fy[j] -= ry * ffac;
                        fz[ii] += rz * ffac;
                        fz[j] -= rz * ffac;
                    }
                    //printf("tid: %d, fx[%d]: %lf, fy[%d]: %lf, fz[%d]: %lf\n", tid, ii,  fx[ii], ii, fy[ii], ii, fz[ii]);

                }
            }
        
        #ifdef _OPENMP
            #pragma omp barrier
        #endif

        // i = (sys->natoms - 1) / sys->tmax;
        i = 1 + (sys->natoms / sys->tmax);
        start = tid * i;
        end = start + i;

        // printf("tid: %d, i: %d, start: %d, end: %d, sys->natoms: %d, sys->tmax: %d\n", tid, i, start, end, sys->natoms, sys->tmax);

        // printf("tid: %d, end: %d, end > sys->natoms: %d\n", tid, end, end > sys->natoms);

        if(end > sys->natoms) end = sys->natoms; 
        
        // printf("tid: %d, end: %d\n", tid, end);


        // #ifdef _OPENMP
        //     #pragma omp for
        // #endif

        for(i = 1; i < sys->tmax; ++i)
        {   
            offset = i * sys->natoms;
            for(j = start; j < end; ++j)
            {   
                sys->fx[j] += sys->fx[offset + j];
                sys->fy[j] += sys->fy[offset + j];
                sys->fz[j] += sys->fz[offset + j];
            }
        }
    }
    sys->epot = epot;
}