// This file provides the interface between Python and the C++
// implementation of the fast marching method in fast_marching.cpp

#include "stdio.h"

#include "Python.h"
#include "numpy/noprefix.h"
#include "fast_marching.h"

static PyObject *distance_method(PyObject *self, PyObject *args);

static PyMethodDef fmm_methods[] =
{
    {"cFastMarcher", (PyCFunction)distance_method, METH_VARARGS,
     "Use the python wrapper to this function"
     "Returns the signed distance or travel time from "
     "the zero level set of phi. "
    },
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initcfmm(void)
{
    PyObject* m;
    m = Py_InitModule3("cfmm", fmm_methods,
        "a c extension for calculating the signed distance and travel "
        "time from the zero level set of a function");
    if (m == NULL)
        return;
    import_array();
}

static PyObject *distance_method(PyObject *self, PyObject *args)
{
  // when we get here we should have:
  // -- phi, dx, flag, and speed
  // -- and the input error checking should be done

  PyObject *pphi, *pdx, *pflag, *pspeed;
  int       self_test;
  PyArrayObject *phi, *dx, *flag, *speed, *distance;

  if (!PyArg_ParseTuple(args, "OOOOi", &pphi, &pdx, &pflag,
                        &pspeed, &self_test))
  {
    return NULL;
  }

  if (! (self_test==0 || self_test==1))
  {
    PyErr_SetString(PyExc_ValueError, "self_test must be 0 or 1");
    return NULL;
  }

  phi = (PyArrayObject *)PyArray_FROMANY(pphi, PyArray_DOUBLE, 1,
                                         10, NPY_IN_ARRAY);
  if (!phi)
  {
    PyErr_SetString(PyExc_ValueError,
                    "phi must be a 1 to 12-D array of doubles");
    return NULL;
  }

  dx = (PyArrayObject *)PyArray_FROMANY(pdx, PyArray_DOUBLE, 1,
                                        1, NPY_IN_ARRAY);
  if (!dx)
  {
    PyErr_SetString(PyExc_ValueError, "dx must be a 1D array of doubles");
    Py_XDECREF(phi);
    return NULL;
  }

  flag = (PyArrayObject *)PyArray_FROMANY(pflag, PyArray_INT, 1,
                                          10, NPY_IN_ARRAY);
  if (!flag)
  {
    PyErr_SetString(PyExc_ValueError,
                    "flag must be a 1D to 12-D array of integers");
    Py_XDECREF(phi);
    Py_XDECREF(dx);
    return NULL;
  }

  if (pspeed != Py_None)
  {
    speed = (PyArrayObject *)PyArray_FROMANY(pspeed, PyArray_DOUBLE, 1,
                                             10, NPY_IN_ARRAY);
    if (!speed)
    {
      PyErr_SetString(PyExc_ValueError,
                      "speed must be a 1D to 12-D array of doubles");
      Py_XDECREF(phi);
      Py_XDECREF(dx);
      Py_XDECREF(flag);
      return NULL;
    }

    if (! PyArray_SAMESHAPE(phi,speed))
    {
      PyErr_SetString(PyExc_ValueError,
                      "phi and speed must have the same shape");
      Py_XDECREF(phi);
      Py_XDECREF(dx);
      Py_XDECREF(flag);
      Py_XDECREF(speed);
      return NULL;
    }
  }
  else speed=0;

  if (! (PyArray_NDIM(phi)==(npy_intp)PyArray_DIM(dx,0))) // ?!
  {
    PyErr_SetString(PyExc_ValueError, "dx must be of length len(phi.shape)");
    Py_XDECREF(phi);
    Py_XDECREF(dx);
    Py_XDECREF(flag);
    Py_XDECREF(speed);
    return NULL;
  }

  for (int i=0; i<PyArray_DIM(dx,0); i++)
  {
    double d = *(double *)PyArray_GETPTR1(dx,i);
    if (d<=0.0)
    {
      PyErr_SetString(PyExc_ValueError, "dx must be greater than zero");
      Py_XDECREF(phi);
      Py_XDECREF(dx);
      Py_XDECREF(flag);
      Py_XDECREF(speed);
      return NULL;
    }
  }

  if (! PyArray_SAMESHAPE(phi,flag))
  {
    PyErr_SetString(PyExc_ValueError, "phi and flag must have the same shape");
    Py_XDECREF(phi);
    Py_XDECREF(dx);
    Py_XDECREF(flag);
    Py_XDECREF(speed);
    return NULL;
  }

  int shape[MaximumDimension];
  npy_intp shape2[MaximumDimension];
  for (int i=0; i<PyArray_NDIM(phi); i++)
  {
    shape[i] = PyArray_DIM(phi,i);
    shape2[i] = PyArray_DIM(phi,i);
  }

  // make a new array for the return value
  distance = (PyArrayObject *)PyArray_ZEROS(PyArray_NDIM(phi),
                                            shape2, PyArray_DOUBLE, 0);
  if (! distance) return NULL;

  // create a level set object to do the calculation
  double * local_phi        = (double *) PyArray_DATA(phi);
  double * local_dx         = (double *) PyArray_DATA(dx);
  int    * local_flag       = (int *)    PyArray_DATA(flag);
  double * local_speed      = 0;
  if (speed) local_speed    = (double *) PyArray_DATA(speed);
  double * local_distance   = (double *) PyArray_DATA(distance);

  fastMarcher *fm = new fastMarcher(
    local_phi,
    local_dx,
    local_flag,
    local_speed,
    local_distance,
    PyArray_NDIM(phi),
    shape,
    self_test);

  int error = fm->getError();
  delete fm;

  Py_DECREF(phi);
  Py_DECREF(flag);
  Py_DECREF(dx);
  Py_XDECREF(speed);

  switch (error)
  {
  case 0:      // no-error
    break;
  case 1:      // unknown error
    // we should never get here
    PyErr_SetString(PyExc_ValueError,
                    "an unknown error has occurred in the scikit-fmm c "
                    "extension module");
    Py_XDECREF(distance);
    return NULL;
  case 2:
    PyErr_SetString(PyExc_ValueError,
                    "the array phi contains no zero contour (no zero level set)");
    Py_XDECREF(distance);
    return NULL;
  }

  // python wrapper adds mask back
  return (PyObject *)distance;
}

