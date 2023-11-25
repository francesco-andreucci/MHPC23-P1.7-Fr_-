// unit test for ekin function
#include "gtest/gtest.h"
#include "memory.h"
#include "types.h"
#include "input.h"

class InputTest: public ::testing::Test {

protected:

    mdsys_t *sys;

    void SetUp()
    {
        sys = new mdsys_t;
    }

    void TearDown()
        {
            delete sys;
        }
};

TEST_F(InputTest, case1)
{
    ASSERT_NE(sys,nullptr);
    int nprint;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *traj,*erg;

    read_from_file(sys, &nprint,restfile,trajfile,ergfile,line);
    memalloc(sys);
    read_restfile(restfile, sys);
    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");

    ASSERT_DOUBLE_EQ(sys->natoms, 2);
    ASSERT_DOUBLE_EQ(sys->mass, 39.948);
    ASSERT_DOUBLE_EQ(sys->epsilon, 0.2379);        
    ASSERT_DOUBLE_EQ(sys->sigma, 3.405);
    ASSERT_DOUBLE_EQ(sys->box, 17.1580);
    ASSERT_DOUBLE_EQ(sys->nsteps, 1);
    ASSERT_DOUBLE_EQ(sys->rcut, 8.5);
    ASSERT_DOUBLE_EQ(sys->dt, 5.0);

    /* atom 0*/
    ASSERT_DOUBLE_EQ(sys->rx[0], 3.00000000);
    ASSERT_DOUBLE_EQ(sys->ry[0], -1.00000000);
    ASSERT_DOUBLE_EQ(sys->rz[0], 4.00000000);

    ASSERT_DOUBLE_EQ(sys->vx[0], 1.00000);
    ASSERT_DOUBLE_EQ(sys->vy[0], 0.1300);
    ASSERT_DOUBLE_EQ(sys->vz[0], 2.65430);

    /* atom 1*/
    ASSERT_DOUBLE_EQ(sys->rx[1], 1.70000);
    ASSERT_DOUBLE_EQ(sys->ry[1], 0.3000);
    ASSERT_DOUBLE_EQ(sys->rz[1], 4.10000);

    ASSERT_DOUBLE_EQ(sys->vx[1], 0.00000);
    ASSERT_DOUBLE_EQ(sys->vy[1], 1.1300);
    ASSERT_DOUBLE_EQ(sys->vz[1], 1.65430);

    cleanup(erg,traj,*sys);
}