// This file provides the interface between Python and the C++
// implementation of the fast marching method in fast_marching.cpp

#include "Python.h"
#include "numpy/noprefix.h"

#include "distance_marcher.h"
#include "travel_time_marcher.h"
#include "extension_velocity_marcher.h"

#include <stdexcept>

#define DISTANCE              0
#define TRAVEL_TIME           1
#define EXTENSION_VELOCITY    2

static PyObject *distance_method(PyObject *self, PyObject *args);

static PyMethodDef fmm_methods[] =
{
    {"cFastMarcher", (PyCFunction)distance_method, METH_VARARGS,
     "Entry point for scikit-fmm c extension"
     "Use the python wrapper to this function"
    },
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "cfmm",                               /* m_name */
    "c extension module for scikit-fmm",  /* m_doc */
    -1,                                   /* m_size */
    fmm_methods,                          /* m_methods */
    NULL,                                 /* m_reload */
    NULL,                                 /* m_traverse */
    NULL,                                 /* m_clear */
    NULL,                                 /* m_free */
};
#endif

PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
PyInit_cfmm(void)
#else
initcfmm(void)
#endif
{
    PyObject* m;

#if PY_MAJOR_VERSION >= 3
    m = PyModule_Create(&moduledef);

    if (m == NULL)
        return NULL;

    import_array();
    return m;
#else
    m = Py_InitModule3("cfmm", fmm_methods,
        "c extension module for scikit-fmm");
    if (m == NULL)
        return;

    import_array();
#endif
}

static PyObject *distance_method(PyObject *self, PyObject *args)
{
  // when we get here we should have:
  // -- phi, dx, flag, and speed
  // -- and the input error checking should be done

  PyObject *pphi, *pdx, *pflag, *pspeed, *pext_mask;
  int       self_test, mode, order, periodic;
  PyArrayObject *phi, *dx, *flag, *speed, *distance, *f_ext, *ext_mask;
  double narrow=0;
  distance = 0;
  f_ext    = 0;
  speed    = 0;
  ext_mask = 0;



  if (!PyArg_ParseTuple(args, "OOOOOiiidi", &pphi, &pdx, &pflag,
                        &pspeed, &pext_mask, &self_test, &mode,
                        &order, &narrow, &periodic))
  {
    return NULL;
  }

  if (! (self_test==0 || self_test==1))
  {
    PyErr_SetString(PyExc_ValueError, "self_test must be 0 or 1");
    return NULL;
  }


  if (! (order==1 || order==2))
  {
    PyErr_SetString(PyExc_ValueError, "order must be 1 or 2");
    return NULL;
  }

  if (! (mode==DISTANCE ||
         mode==TRAVEL_TIME ||
         mode==EXTENSION_VELOCITY))
  {
    PyErr_SetString(PyExc_ValueError, "invalid mode flag");
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

  flag = (PyArrayObject *)PyArray_FROMANY(pflag, PyArray_LONG, 1,
                                          10, NPY_IN_ARRAY);
  if (!flag)
  {
    PyErr_SetString(PyExc_ValueError,
                    "flag must be a 1D to 12-D array of integers");
    Py_XDECREF(phi);
    Py_XDECREF(dx);
    return NULL;
  }

  if (mode == TRAVEL_TIME || mode == EXTENSION_VELOCITY)
  {
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
  }

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

  if (mode == EXTENSION_VELOCITY)
  {
    f_ext = (PyArrayObject *)PyArray_ZEROS(PyArray_NDIM(phi),
                                           shape2, PyArray_DOUBLE, 0);
    if (! f_ext) return NULL;

    ext_mask = (PyArrayObject *)PyArray_FROMANY(pext_mask, PyArray_LONG, 1,
                                                10, NPY_IN_ARRAY);
    if (! ext_mask)
      {
        PyErr_SetString(PyExc_ValueError,
                        "ext_mask must be a 1D, 12-D array of integers");
        Py_XDECREF(phi);
        Py_XDECREF(dx);
        Py_XDECREF(flag);
        Py_XDECREF(speed);
        return NULL;
      }

  }

  // create a level set object to do the calculation
  double * local_phi        = (double *) PyArray_DATA(phi);
  double * local_dx         = (double *) PyArray_DATA(dx);
  long   * local_flag       = (long *)   PyArray_DATA(flag);
  long    * local_ext_mask   = 0;
  if (ext_mask) local_ext_mask = (long *) PyArray_DATA(ext_mask);
  double * local_speed      = 0;
  if (speed) local_speed    = (double *) PyArray_DATA(speed);
  double * local_distance   = (double *) PyArray_DATA(distance);
  int error;

  baseMarcher *marcher = 0;
  switch (mode)
  {
    case DISTANCE:
    {
      marcher = new distanceMarcher(
        local_phi,
        local_dx,
        local_flag,
        local_distance,
        PyArray_NDIM(phi),
        shape,
        self_test,
        order,
        narrow,
        periodic);
    }
    break;
    case TRAVEL_TIME:
    {
      marcher = new travelTimeMarcher(
        local_phi,
        local_dx,
        local_flag,
        local_distance,
        PyArray_NDIM(phi),
        shape,
        self_test,
        order,
        local_speed,
        narrow,
        periodic);
    }
    break;
    case EXTENSION_VELOCITY:
    {
      double * local_fext = (double *) PyArray_DATA(f_ext);
      marcher = new extensionVelocityMarcher(
        local_phi,
        local_dx,
        local_flag,
        local_distance,
        PyArray_NDIM(phi),
        shape,
        self_test,
        order,
        local_ext_mask,
        local_speed,
        local_fext,
        narrow,
        periodic);
    }
    break;
  default: error=1;
  }

  try {
      marcher->march();
      error = marcher->getError();
      delete marcher;
    } catch (const std::exception& exn) {
      // propagate error
      PyErr_SetString(PyExc_RuntimeError, exn.what());
      Py_XDECREF(phi);
      Py_XDECREF(dx);
      Py_XDECREF(flag);
      Py_XDECREF(speed);
      Py_XDECREF(ext_mask);
      return NULL;
    }

  Py_DECREF(phi);
  Py_DECREF(flag);
  Py_DECREF(dx);
  Py_XDECREF(speed);
  Py_XDECREF(ext_mask);

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

  if (mode == EXTENSION_VELOCITY)
  {
    return Py_BuildValue("NN", distance, f_ext);
  }
  return (PyObject *)distance;
}
