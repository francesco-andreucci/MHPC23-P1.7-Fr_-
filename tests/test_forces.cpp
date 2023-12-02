#include "gtest/gtest.h"
#include "types.h"
#include "force_comp.h"
#include "utilities.h"


#ifdef LJMD_MPI
#include "mpi.h"
#endif

#ifdef _OPENMP
    #include "omp.h"
#endif

class ForceTest: public ::testing::Test {

protected:

    mdsys_t *sys;

    void SetUp()
    {
        sys = new mdsys_t;
        #ifdef _OPENMP
        #pragma omp parallel
            sys->tmax = omp_get_max_threads();
            int tid = omp_get_thread_num();
        #else
            sys->tmax = 1;
            int tid = 0;
        #endif
        sys->epsilon=0.237900;
        sys->sigma=3.405000;
        sys->rcut=8.5;
        sys->natoms = 2;
        sys->box=17.158000;
        sys->dt = 5.0;
        sys->mass = 1.0 ;
        sys->mpirank=0;
        sys->nsize=1;
        sys->rx = new double[2];
        sys->vx = new double[2];
        sys->fx = new double[2*sys->tmax];
        sys->ry = new double[2];
        sys->vy = new double[2];
        sys->fy = new double[2*sys->tmax];
        sys->rz = new double[2];
        sys->vz = new double[2];
        sys->fz = new double[2*sys->tmax];
#ifdef LJMD_MPI
        sys->cx = new double[2*sys->tmax];
        sys->cy = new double[2*sys->tmax];
        sys->cz = new double[2*sys->tmax];
#endif
        sys->rx[0] = 0.0;
        sys->rx[1] = 0.0;
        sys->ry[0] = 0.0;
        sys->ry[1] = 0.0;
        sys->rz[0] = 0.0;
        sys->rz[1] = 0.0;
        sys->vx[0] = 0.0;
        sys->vx[1] = 0.0;
        sys->vy[0] = 0.0;
        sys->vy[1] = 0.0;
        sys->vz[0] = 0.0;
        sys->vz[1] = 0.0;
        sys->fx[0] = 0.0;
        sys->fx[1] = 0.0;
        sys->fy[0] = 0.0;
        sys->fy[1] = 0.0;
        sys->fz[0] = 0.0;
        sys->fz[1] = 0.0;

#ifdef LJMD_MPI
        sys->cx[0] = 0.0;
        sys->cx[1] = 0.0;
        sys->cy[0] = 0.0;
        sys->cy[1] = 0.0;
        sys->cz[0] = 0.0;
        sys->cz[1] = 0.0;
#endif
}

    void TearDown()
        {
            delete[] sys->rx;
            delete[] sys->vx;
            delete[] sys->fx;
            delete[] sys->ry;
            delete[] sys->vy;
            delete[] sys->fy;
            delete[] sys->rz;
            delete[] sys->vz;
            delete[] sys->fz;
        #ifdef LJMD_MPI
            delete[] sys->cx;
            delete[] sys->cz;
            delete[] sys->cy;
        #endif
            delete sys;
        }
};

TEST_F(ForceTest, force0)
{
//The particles are further apart than rcut, but inside the box: the force should be zero
#ifdef LJMD_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &sys->mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &sys->nsize);
#endif

    // define tollerance
    double toler = 10e-12;
    //The particles are further apart than rcut, but inside the box: the force should be zero
    sys->rx[0]=1.0;
    sys->rx[1]=1.1+sys->rcut;
    sys->ry[0]=1.0;
    sys->ry[1]=1.1+sys->rcut;
    sys->rz[0]=1.0;
    sys->rz[1]=1.1+sys->rcut;
if(sys->mpirank==0){
    ASSERT_NE(sys,nullptr);
    ASSERT_DOUBLE_EQ(sys->rx[0], 1.0);
    ASSERT_DOUBLE_EQ(sys->rx[1],1.1+sys->rcut);
    ASSERT_DOUBLE_EQ(sys->ry[0], 1.0);
    ASSERT_DOUBLE_EQ(sys->ry[1],1.1+sys->rcut);
    ASSERT_DOUBLE_EQ(sys->rz[0], 1.0);
    ASSERT_DOUBLE_EQ(sys->rz[1],1.1+sys->rcut);
}
    force(sys);
if(sys->mpirank==0){
    ASSERT_DOUBLE_EQ(sys->fx[0],0.0);
    ASSERT_DOUBLE_EQ(sys->fx[1],0.0);
    ASSERT_DOUBLE_EQ(sys->fy[0],0.0);
    ASSERT_DOUBLE_EQ(sys->fy[1],0.0);
    ASSERT_DOUBLE_EQ(sys->fz[0],0.0);
    ASSERT_DOUBLE_EQ(sys->fz[1],0.0);
}
}


TEST_F(ForceTest, forceno0)
{
    // define tolleranc
    double toler = 10e-12;

/*In this case the particles are close enough to interact: we compare the
 * output of the force command with reference values.
 */
    
    #ifdef LJMD_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &sys->mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &sys->nsize);
    #endif

    /*In this case the particles are close enough to interact: we compare the
    * output of the force command with reference values. */
    sys->rx[0]=1.0;
    sys->rx[1]=2.1;
    sys->ry[0]=1.0;
    sys->ry[1]=2.2;
    sys->rz[0]=1.0;
    sys->rz[1]=2.3;
if(sys->mpirank==0){
    ASSERT_NE(sys,nullptr);
    ASSERT_DOUBLE_EQ(sys->rx[0], 1.0);
    ASSERT_DOUBLE_EQ(sys->rx[1],2.1);
    ASSERT_DOUBLE_EQ(sys->ry[0], 1.0);
    ASSERT_DOUBLE_EQ(sys->ry[1],2.2);
    ASSERT_DOUBLE_EQ(sys->rz[0], 1.0);
    ASSERT_DOUBLE_EQ(sys->rz[1],2.3);
}
    force(sys);
if(sys->mpirank==0){
    EXPECT_NEAR(sys->fx[0],-1024.386173537659715, toler);
    EXPECT_NEAR(sys->fx[1],1024.386173537659715, toler);
    EXPECT_NEAR(sys->fy[0],-1117.512189313810723, toler);
    EXPECT_NEAR(sys->fy[1],1117.512189313810723, toler);
    EXPECT_NEAR(sys->fz[0],-1210.638205089961275, toler);
    EXPECT_NEAR(sys->fz[1],1210.638205089961275, toler);
}
}


TEST_F(ForceTest, forceno0pbc)
{
    // define tollerance
    double toler = 10e-12;

/*In this case the particles lie on the same yz plane: the interaction is only along the x direction, but the interparticle distance is nominally larger than box edge.
 */ #ifdef LJMD_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &sys->mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &sys->nsize);
   #endif

    sys->rx[0]=3.7;
    sys->rx[1]=23.558;
    sys->ry[0]=4.6;
    sys->ry[1]=4.6;
    sys->rz[0]=1.9;
    sys->rz[1]=1.9;

if(sys->mpirank==0){
    ASSERT_NE(sys,nullptr);
    ASSERT_DOUBLE_EQ(sys->rx[0], 3.7);
    ASSERT_DOUBLE_EQ(sys->rx[1],23.558);
    ASSERT_DOUBLE_EQ(sys->ry[0], 4.6);
    ASSERT_DOUBLE_EQ(sys->ry[1],4.6);
    ASSERT_DOUBLE_EQ(sys->rz[0], 1.9);
    ASSERT_DOUBLE_EQ(sys->rz[1],1.9);
}
    force(sys);
if(sys->mpirank==0){
    EXPECT_NEAR(sys->fx[0],-59.933619233522784, toler);
    EXPECT_NEAR(sys->fx[1],59.933619233522784, toler);
    EXPECT_NEAR(sys->fy[0],0.0, toler);
    EXPECT_NEAR(sys->fy[1],0.0, toler);
    EXPECT_NEAR(sys->fz[0],0.0, toler);
    EXPECT_NEAR(sys->fz[1],0.0, toler);
}
}

int main(int argc, char** argv) {

    ::testing::InitGoogleTest(&argc, argv);
#ifdef LJMD_MPI
    MPI_Init(&argc, &argv);
#endif
    int rank, size;
    #ifdef LJMD_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
   #endif


    int result = RUN_ALL_TESTS();
    printf("%d\n", result);




    #ifdef LJMD_MPI
        MPI_Finalize();
    #endif

    return 0;
}
