import ctypes as ct
import sys
import os 

#define some staff
BLEN = 256

# mimic c _mdsys structure
class mdlib_sys(ct.Structure):
   _fields_ = [
      ('nprint', ct.c_int),
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
      ("restfile", ct.c_char_p),
      ("trajfile", ct.c_char_p),
      ("ergfile", ct.c_char_p),
   ]

# Load the mdlib shared library
mdlib = ct.CDLL('../build/libmdlib.dylib') 

# Call the ekin
tsys = mdlib_sys()

# Declare the function signatures

# allocate memory
memalloc = mdlib.memalloc
memalloc.argtypes = [ct.POINTER(mdlib_sys)]
memalloc.restype = None

# kenetic energy
ekin = mdlib.ekin
ekin.argtypes = [ct.POINTER(mdlib_sys)]
ekin.restype = None

# force function
force = mdlib.force
force.argtypes = [ct.POINTER(mdlib_sys)]
force.restype = None

# azzero function
azzero = mdlib.azzero
azzero.argtypes = [ct.POINTER(ct.c_double), ct.c_int]
azzero.restype = None


# velverlet1 function
velverlet1 = mdlib.velverlet1
velverlet1.argtypes = [ct.POINTER(mdlib_sys)]
velverlet1.restype = None

# ververlet2 function
velverlet2 = mdlib.velverlet2
velverlet2.argtypes = [ct.POINTER(mdlib_sys)]
velverlet2.restype = None

# output function
output = mdlib.output
output.argtypes = (
      ct.POINTER(mdlib_sys),
      ct.POINTER(ct.c_void_p),
      ct.POINTER(ct.c_void_p), 
    )
output.restype = None

# clean up function
cleanup = mdlib.cleanup
cleanup.argtypes = [ct.POINTER(ct.c_void_p), ct.POINTER(ct.c_void_p), ct.POINTER(mdlib_sys)]
cleanup.restype = None

# Allocate memory
memalloc(ct.byref(tsys))

# set nfi to zero
tsys.nfi=0

# read from file 
def read_file(inputfile):
    file = open(inputfile,"r")
    row = [line.rstrip('\n').split(" ")[0] for line in file]

    tsys.natoms = ct.c_int(int(row[0]))
    tsys.mass = ct.c_double(float(row[1]))
    tsys.epsilon = ct.c_double(float(row[2]))
    tsys.sigma = ct.c_double(float(row[3]))
    tsys.rcut = ct.c_double(float(row[4]))
    tsys.box = ct.c_double(float(row[5]))
    tsys.restfile = "".join(row[6].split()).encode()
    tsys.trajfile = "".join(row[7].split()).encode()
    tsys.ergfile = "".join(row[8].split()).encode()
    tsys.nsteps = ct.c_int(int(float(row[9])))
    tsys.dt = ct.c_double(float(row[10]))
    tsys.nprint = ct.c_int(int(row[11]))



if __name__ == '__main__':
	
	if len(sys.argv) != 1:
		inputfile = str(sys.argv[1])
		read_file(inputfile)
	else:
            print("Usage: python3 main.py input_file.inp")

# Test reading file
print(f"natoms: {tsys.natoms}")
print(f"mass: {tsys.mass}")
print(f"epsilon: {tsys.epsilon}")
print(f"sigma: {tsys.sigma}")
print(f"rcut: {tsys.rcut}")
print(f"box: {tsys.box}")
print(f"restfiles: {tsys.restfile}")
print(f"trajfile: {tsys.trajfile}")
print(f"ergfile: {tsys.ergfile}")
print(f"nsteps: {tsys.nsteps}")
print(f"dt: {tsys.dt}")
print(f"nprint: {tsys.nprint}")