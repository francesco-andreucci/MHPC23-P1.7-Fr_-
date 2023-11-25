#include "gtest/gtest.h"
#include "constants.h"
#include "verlet.h"
#include "types.h"


class VerletTest: public ::testing::Test {

    protected:
       mdsys_t *sys;

       void SetUp() {
           sys = new mdsys_t;
           sys -> natoms=2;
           sys -> dt=5.0;
           sys -> mass = sys-> dt/mvsq2e;
           sys -> rx = new double [2];
           sys -> ry = new double [2];
           sys -> rz = new double [2];
           sys -> vx = new double [2];
           sys -> vy = new double [2];
           sys -> vz = new double [2];
           sys -> fx = new double [2];
           sys -> fy = new double [2];
           sys -> fz = new double [2];

           //initialization
           sys -> rx [0] = -1.0;
           sys -> ry [0] = 0.0;
           sys -> rz [0] = 0.0;
           sys -> vx [0] = 0.0;
           sys -> vy [0] = 0.0;
           sys -> vz [0] = 0.0;
           sys -> fx [0] = 1.0;
           sys -> fy [0] = 0.0;
           sys -> fz [0] = 0.0;

           sys -> rx [1] = 0.0;
           sys -> ry [1] = 0.0;
           sys -> rz [1] = 0.0;
           sys -> vx [1] = 0.5;
           sys -> vy [1] = 1.0;
           sys -> vz [1] = 1.5;
           sys -> fx [1] = 2.0;
           sys -> fy [1] = 1.0;
           sys -> fz [1] = 3.0;
       }

       void TearDown(){
           delete[] sys -> rx;
           delete[] sys -> ry;
           delete[] sys -> rz;
           delete[] sys -> vx;
           delete[] sys -> vy;
           delete[] sys -> vz;
           delete[] sys -> fx;
           delete[] sys -> fy;
           delete[] sys -> fz;
           delete[] sys;
       }


};


TEST_F(VerletTest, Verlet1)
       {
        ASSERT_NE(sys, nullptr);
        ASSERT_DOUBLE_EQ(sys -> rx[0], -1.0);
        ASSERT_DOUBLE_EQ(sys -> vx[0], 0.0);
        ASSERT_DOUBLE_EQ(sys -> fx[0], 1.0);

        ASSERT_DOUBLE_EQ(sys -> rx[1], 0.0);
        ASSERT_DOUBLE_EQ(sys -> vx[1], 0.5);
        ASSERT_DOUBLE_EQ(sys -> fx[1], 2.0);
        ASSERT_DOUBLE_EQ(sys -> ry[1], 0.0);
        ASSERT_DOUBLE_EQ(sys -> vy[1], 1.0);
        ASSERT_DOUBLE_EQ(sys -> fy[1], 1.0);
        ASSERT_DOUBLE_EQ(sys -> rz[1], 0.0);
        ASSERT_DOUBLE_EQ(sys -> vz[1], 1.5);
        ASSERT_DOUBLE_EQ(sys -> fz[1], 3.0);

        velverlet1(sys);

        ASSERT_DOUBLE_EQ(sys -> rx[0], 1.5 );
        ASSERT_DOUBLE_EQ(sys -> vx[0], 0.5 );
        ASSERT_DOUBLE_EQ(sys -> fx[0], 1.0 );

        ASSERT_DOUBLE_EQ(sys -> rx[1], 7.5);
        ASSERT_DOUBLE_EQ(sys -> vx[1], 1.5);
        ASSERT_DOUBLE_EQ(sys -> fx[1], 2.0);
        ASSERT_DOUBLE_EQ(sys -> ry[1], 7.5);
        ASSERT_DOUBLE_EQ(sys -> vy[1], 1.5);
        ASSERT_DOUBLE_EQ(sys -> fy[1], 1.0);
        ASSERT_DOUBLE_EQ(sys -> rz[1], 15.0);
        ASSERT_DOUBLE_EQ(sys -> vz[1], 3.0);
        ASSERT_DOUBLE_EQ(sys -> fz[1], 3.0);
       }


TEST_F(VerletTest, Verlet2)
       {
        ASSERT_NE(sys, nullptr);
        ASSERT_DOUBLE_EQ(sys -> rx[0], -1.0);
        ASSERT_DOUBLE_EQ(sys -> vx[0], 0.0);
        ASSERT_DOUBLE_EQ(sys -> fx[0], 1.0);

        ASSERT_DOUBLE_EQ(sys -> rx[1], 0.0);
        ASSERT_DOUBLE_EQ(sys -> vx[1], 0.5);
        ASSERT_DOUBLE_EQ(sys -> fx[1], 2.0);
        ASSERT_DOUBLE_EQ(sys -> ry[1], 0.0);
        ASSERT_DOUBLE_EQ(sys -> vy[1], 1.0);
        ASSERT_DOUBLE_EQ(sys -> fy[1], 1.0);
        ASSERT_DOUBLE_EQ(sys -> rz[1], 0.0);
        ASSERT_DOUBLE_EQ(sys -> vz[1], 1.5);
        ASSERT_DOUBLE_EQ(sys -> fz[1], 3.0);


        velverlet2(sys);


        ASSERT_DOUBLE_EQ(sys -> rx[0], -1.0);
        ASSERT_DOUBLE_EQ(sys -> vx[0], 0.5);
        ASSERT_DOUBLE_EQ(sys -> fx[0], 1.0);

        ASSERT_DOUBLE_EQ(sys -> rx[1], 0.0);
        ASSERT_DOUBLE_EQ(sys -> vx[1], 1.5);
        ASSERT_DOUBLE_EQ(sys -> fx[1], 2.0);
        ASSERT_DOUBLE_EQ(sys -> ry[1], 0.0);
        ASSERT_DOUBLE_EQ(sys -> vy[1], 1.5);
        ASSERT_DOUBLE_EQ(sys -> fy[1], 1.0);
        ASSERT_DOUBLE_EQ(sys -> rz[1], 0.0);
        ASSERT_DOUBLE_EQ(sys -> vz[1], 3.0);
        ASSERT_DOUBLE_EQ(sys -> fz[1], 3.0);
       }

