# LJMD python main
import ctypes as ct
# from mpi4py import MPI

# mimic c _mdsys structure
class mdlib_sys(ct.Structure):
    _fields_ = [
                ('natoms', ct.c_int),
                ('nfi', ct.c_int),
                ('nsteps', ct.c_int),
                ('dt', ct.c_double),
                ('mass', ct.c_double),
                ('epsilon', ct.c_double),
                ('sigma', ct.c_double),
                ('box', ct.c_double),
                ('rcut', ct.c_double),
                ('ekin', ct.c_double),
                ('epot', ct.c_double),
                ('temp', ct.c_double),
                ("rx", ct.POINTER(ct.c_double)),
                ("ry", ct.POINTER(ct.c_double)),
                ("rz", ct.POINTER(ct.c_double)),
                ("vx", ct.POINTER(ct.c_double)),
                ("vy", ct.POINTER(ct.c_double)),
                ("vz", ct.POINTER(ct.c_double)),
                ("fx", ct.POINTER(ct.c_double)),
                ("fy", ct.POINTER(ct.c_double)),
                ("fz", ct.POINTER(ct.c_double)),
                ("nsize", ct.c_int),
                ("mpirank", ct.c_int),
               ]
    
# Load the mdlib shared library
mdlib = ct.CDLL('../build/libmdlib.dylib') 

# Call the ekin
sys = mdlib_sys()

# Declare the ekin function signature
ekin = mdlib.ekin
ekin.argtypes = [ct.POINTER(mdlib_sys)]
ekin.restype = None

# Allocate memory
sys.rx = (ct.c_double * sys.nsize)()
sys.ry = (ct.c_double * sys.nsize)()
sys.rz = (ct.c_double * sys.nsize)()
sys.vx = (ct.c_double * sys.nsize)()
sys.vy = (ct.c_double * sys.nsize)()
sys.vz = (ct.c_double * sys.nsize)()
sys.fx = (ct.c_double * sys.nsize)()
sys.fy = (ct.c_double * sys.nsize)()
sys.fz = (ct.c_double * sys.nsize)()

# set up
sys.natoms = 2
sys.mass = 39.948000000000000

#/* Velocities for atom 0 */
sys.vx[0] = 1.181201103804874
sys.vy[0] = -0.051201103804874
sys.vz[0] = 2.640361453553472

#/* Velocities for atom 1 */
sys.vx[1] = -0.181201103804874
sys.vy[1] = 1.311201103804873
sys.vz[1] = 1.668238546446529

# Call the ekin function
ekin(ct.byref(sys))

# Access the updated values in sys
print("Ekin:", sys.ekin)
print("Temp:", sys.temp)