// unit test for ekin function
#include "gtest/gtest.h"
#include "types.h"
#include "utilities.h"
#include "constants.h"
#include <sys/time.h>

class EkinTest: public ::testing::Test {

protected:

    mdsys_t *sys;

    void SetUp()
    {
        sys = new mdsys_t;
        sys->natoms = 2;
        sys->mass = 1.0;
        sys->ekin = 0.0;
        sys->vx = new double[2];
        sys->vy = new double[2];
        sys->vz = new double[2];

        /* atom 1 velocities */   
        sys->vx[0] = 0.000000;
        sys->vy[0] = 0.000000;
        sys->vz[0] = 0.000000;

         /* atom 2 velocities */   
        sys->vx[1] = 0.000000;
        sys->vy[1] = 0.000000;
        sys->vz[1] = 0.000000;
    }

    void TearDown()
        {
            delete[] sys->vx;
            delete[] sys->vy;
            delete[] sys->vz;

            delete sys;
        }
};

TEST_F(EkinTest, case1)
{
    ASSERT_NE(sys,nullptr);
    ASSERT_DOUBLE_EQ(sys->ekin, 0.000000);
    ASSERT_DOUBLE_EQ(sys->vx[0], 0.000000);
    ASSERT_DOUBLE_EQ(sys->vy[0], 0.000000);
    ASSERT_DOUBLE_EQ(sys->vz[0], 0.000000);
    ekin(sys);
    ASSERT_DOUBLE_EQ(sys->ekin, 0.000000);
    ASSERT_DOUBLE_EQ(sys->temp, 0.000000);
}