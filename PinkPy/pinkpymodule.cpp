/*  Author: Joseph Glavan    j.glavan4@gmail.com
    
    Module: PinkPy
    Generates 3-dimensional 1/f pink noise, returning it as an array that can then be passed to Python
    
    Adapted from an unidentified source - thank you mystery person!
*/

#include <Python.h> // import Python API
#include <arrayobject.h> // need to import numpy array API --> I think this gets covered in the init function?????
#include "pinkGen.h"
#include "pinkGen.cpp"

int NX, NY, NZ, SeedVal, NSIZE;
float StdSExp, StdTExp;
long NN[4];
float * FFT_Data;


static PyObject * pink_noise_wrapper (PyObject *self, PyObject *args) { // from some witchcraft in the symbol table this will get called from Python (I think)
    // declare I/O variables
    PyObject *ret;
    
    // Parse Python tuple or dictionary of arguments into a C pointer
    if (!PyArg_ParseTuple(args, "iiiffi", &NX, &NY, &NZ, &StdSExp, &StdTExp, &SeedVal))
        return NULL;
    
    // Do C stuff
    npy_intp ret_shape[] = {NX, NY, NZ};  // should be an array of the length of each dimension e.g. [window_width, window_height, nFrames]
    NSIZE = (NX * NY * NZ);
    NN[0]=0; NN[1]=NX; NN[2]=NY; NN[3]=NZ;
    InitData();
    PrepareStimulus(FFT_Data,StdSExp,StdTExp);
    complex_to_float(FFT_Data, NSIZE);
    
    // Convert return value from C to Python object (more specifically, a numpy ndarray)
    float * reduced = new float[NSIZE];
    for (int i=0;i<NSIZE;i++)
        reduced[i] = FFT_Data[i*2];
    ret = PyArray_SimpleNewFromData(3, ret_shape, NPY_FLOAT, reduced);
    
    // Clean up
    free(FFT_Data);
    delete[] reduced;
    return ret;
}

static PyMethodDef PinkMethods[] = {
    {"pinkNoise", pink_noise_wrapper, METH_VARARGS, "int nx,ny,nz; float s,t; int seed; \n Generate 3D 1/f pink noise, returned as a numpy ndarray."},
        {NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initpinkpy (void) {
    (void) Py_InitModule("pinkpy", PinkMethods);
    import_array();
}
