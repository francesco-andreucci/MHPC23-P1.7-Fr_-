#include <math.h>
#include "types.h"
#include "utilities.h"

#ifdef LJMD_OMP
    #include "omp.h"
#endif

void force(mdsys_t *sys) {
    double r, rsq, ffac;
    double rx, ry, rz;
    double rx1, ry1, rz1;

    int i, j;

    double epot = 0.0;
    #ifdef LJMD_OMP
        tmax = omp_get_max_threads();
        tid = omp_get_thread_num();
    #else
        tid = 0;
        tmax = 1;
    #endif

        /* zero energy and forces */
        double r6, rinv;
        double c12 = 4.0 * sys->epsilon * pow(sys->sigma, 12.0);
        double c6 = 4.0 * sys->epsilon * pow(sys->sigma, 6.0);
        double rcsq = sys->rcut * sys->rcut;

    #ifdef LJMD_OMP
        #pragma omp parallel reduction(+:epot) //private(rx,ry,rz,rsq)
    #endif
    { 
        double *fx, *fy, *fz;
        int start, end;
        epot = 0.0;

        fx = sys->fx + tid * sys->natoms;
        fy = sys->fy + tid * sys->natoms;
        fz = sys->fz + tid * sys->natoms;

        azzero(fx, sys->natoms);
        azzero(fy, sys->natoms);
        azzero(fz, sys->natoms);

        for(int i = 0; i < (sys->natoms - 1); i += tmax)
        {   
            int ii = i + tid;
            if(ii >= (sys->natoms - 1)) break;

            rx1=sys->rx[ii];
            ry1=sys->ry[ii];
            rz1=sys->rz[ii];
            
            for(j = ii + 1; j < sys->natoms; ++j) 
            {
                rx = pbc(rx1 - sys->rx[j], 0.5 * sys->box);
                ry = pbc(ry1 - sys->ry[j], 0.5 * sys->box);
                rz = pbc(rz1 - sys->rz[j], 0.5 * sys->box);
                rsq = (rx * rx) + (ry * ry) + (rz * rz);
                
                if(rsq < rcsq) 
                {
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
            }
        }
            

    }
    sys->epot = epot;
}