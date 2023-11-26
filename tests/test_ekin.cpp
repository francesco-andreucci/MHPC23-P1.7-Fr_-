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
        sys->mass = 39.948000000000000;
        sys->ekin = 0.0;
        sys->vx = new double[2];
        sys->vy = new double[2];
        sys->vz = new double[2];
    }

    void TearDown()
        {
            delete[] sys->vx;
            delete[] sys->vy;
            delete[] sys->vz;

            delete sys;
        }
};

TEST_F(EkinTest, scenario1)
{   
    /* Velocities for atom 0 */
    sys->vx[0] = 1.181201103804874;
    sys->vy[0] = -0.051201103804874;
    sys->vz[0] = 2.640361453553472;

    /* Velocities for atom 1 */
    sys->vx[1] = -0.181201103804874;
    sys->vy[1] = 1.311201103804873;
    sys->vz[1] = 1.668238546446529;

    /* Checking the initialization of the structure from SetUp */
    ASSERT_NE(sys,nullptr);
    ASSERT_DOUBLE_EQ(sys->ekin, 0.000000);

    /* Double checking init of velocities for atomes 0 and 1 */
    // atom 0
    ASSERT_DOUBLE_EQ(sys->vx[0], 1.181201103804874);
    ASSERT_DOUBLE_EQ(sys->vy[0], -0.051201103804874);
    ASSERT_DOUBLE_EQ(sys->vz[0], 2.640361453553472);

    // atom 1
    ASSERT_DOUBLE_EQ(sys->vx[1],  -0.181201103804874);
    ASSERT_DOUBLE_EQ(sys->vy[1], 1.311201103804873);
    ASSERT_DOUBLE_EQ(sys->vz[1], 1.668238546446529);

    // calling ekin function
    ekin(sys);

    // checking our expectation from ekin function
    ASSERT_DOUBLE_EQ(sys->ekin, 616046.482825497281738);
    ASSERT_DOUBLE_EQ(sys->temp, 206670828.564033240079880);
}


TEST_F(EkinTest, scenario2)
{   
    /* Velocities for atom 0 */
    sys->vx[0] = 1.099999173846238;
    sys->vy[0] = 1.129999938803425;
    sys->vz[0] = -1.654301407521224;

    /* Velocities for atom 1 */
    sys->vx[1] = 0.300000826153762;
    sys->vy[1] = 3.130000061196575;
    sys->vz[1] = 2.654301407521224;

    /* Checking the initialization of the structure from SetUp */
    ASSERT_NE(sys, nullptr);
    ASSERT_DOUBLE_EQ(sys->ekin, 0.000000);

    /* Double checking init of velocities for atomes 0 and 1 */
    // atom 0
    ASSERT_DOUBLE_EQ(sys->vx[0], 1.099999173846238);
    ASSERT_DOUBLE_EQ(sys->vy[0], 1.129999938803425);
    ASSERT_DOUBLE_EQ(sys->vz[0], -1.654301407521224);

    // atom 1
    ASSERT_DOUBLE_EQ(sys->vx[1], 0.300000826153762);
    ASSERT_DOUBLE_EQ(sys->vy[1], 3.130000061196575);
    ASSERT_DOUBLE_EQ(sys->vz[1], 2.654301407521224);

    // calling ekin function
    ekin(sys);

    // checking our expectation from ekin function
    ASSERT_DOUBLE_EQ(sys->ekin, 1057697.201570168137550);
    ASSERT_DOUBLE_EQ(sys->temp, 354835492.308598577976227);
}