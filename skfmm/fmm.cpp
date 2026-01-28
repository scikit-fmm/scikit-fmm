// This file provides the interface between Python and the C++
// implementation of the fast marching method in fast_marching.cpp

#include "Python.h"
#include "numpy/ndarrayobject.h"

#ifndef NPY_IN_ARRAY
#define NPY_IN_ARRAY (NPY_ARRAY_ALIGNED | NPY_ARRAY_ENSUREARRAY | NPY_ARRAY_NOTSWAPPED | NPY_ARRAY_BEHAVED)
#endif

#include "distance_marcher.h"
#include "travel_time_marcher.h"
#include "travel_time_marcher_genes.h"
#include "extension_velocity_marcher.h"

#include <stdexcept>
#include <cstdio>

#define DISTANCE              0
#define TRAVEL_TIME           1
#define EXTENSION_VELOCITY    2
#define TRAVEL_TIME_GENES     3

static PyObject* distance_method(PyObject* self, PyObject* args);

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
	  // TODO are these NULL members pointers or not?
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

    if (m == nullptr)
        return nullptr;

    import_array();
    return m;
#else
    m = Py_InitModule3("cfmm", fmm_methods,
        "c extension module for scikit-fmm");
    if (m == nullptr)
        return;

    import_array();
#endif
}

static PyObject* distance_method(PyObject* self, PyObject* args)
{
  // when we get here we should have:
  // -- phi, dx, flag, and speed
  // -- and the input error checking should be done

  PyObject* pphi;
  PyObject* pdx;
  PyObject* pflag;
  PyObject* pspeed;
  PyObject* pext_mask;
  PyObject* pspeeds = nullptr;
  PyObject* pdrivers = nullptr;

  int       self_test, mode, order, periodic;

  PyArrayObject* phi;
  PyArrayObject* dx;
  PyArrayObject* flag;

  PyArrayObject* speed = nullptr;
  PyArrayObject* distance = nullptr;
  PyArrayObject* f_ext = nullptr;
  PyArrayObject* ext_mask = nullptr;
  PyArrayObject* speeds = nullptr; // for genetics extension
  PyArrayObject* drivers = nullptr; // for genetics extension

  double narrow = 0;

  // Read the arguments given to the Python interpreter into the PyObject
  // pointers above:
  if (!PyArg_ParseTuple(args, "OOOOOiiidi|OO", &pphi, &pdx, &pflag,
                        &pspeed, &pext_mask, &self_test, &mode,
                        &order, &narrow, &periodic, &pspeeds, &pdrivers))
  {
    return nullptr;
  }

  if (! (self_test==0 || self_test==1))
  {
    PyErr_SetString(PyExc_ValueError, "self_test must be 0 or 1");
    return nullptr;
  }


  if (! (order==1 || order==2))
  {
    PyErr_SetString(PyExc_ValueError, "order must be 1 or 2");
    return nullptr;
  }

  if (! (mode==DISTANCE ||
         mode==TRAVEL_TIME ||
         mode==EXTENSION_VELOCITY ||
         mode==TRAVEL_TIME_GENES))
  {
    PyErr_SetString(PyExc_ValueError, "invalid mode flag");
    return nullptr;
  }

  phi = (PyArrayObject *)PyArray_FROMANY(pphi, NPY_DOUBLE, 1,
                                         12, NPY_IN_ARRAY);
  if (!phi)
  {
    PyErr_SetString(PyExc_ValueError,
                    "phi must be a 1 to 12-D array of doubles");
    return nullptr;
  }

  dx = (PyArrayObject *)PyArray_FROMANY(pdx, NPY_DOUBLE, 1,
                                        1, NPY_IN_ARRAY);
  if (!dx)
  {
    PyErr_SetString(PyExc_ValueError, "dx must be a 1D array of doubles");
    Py_XDECREF(phi);
    return nullptr;
  }

  flag = (PyArrayObject *)PyArray_FROMANY(pflag, NPY_LONGLONG, 1,
                                          12, NPY_IN_ARRAY);
  if (!flag)
  {
    PyErr_SetString(PyExc_ValueError,
                    "flag must be a 1D to 12-D array of integers");
    Py_XDECREF(phi);
    Py_XDECREF(dx);
    return nullptr;
  }

  if (mode == TRAVEL_TIME || mode == EXTENSION_VELOCITY)
	{
		speed = (PyArrayObject *)PyArray_FROMANY(pspeed, NPY_DOUBLE, 1, 12, NPY_IN_ARRAY);
		if (!speed)
		{
			PyErr_SetString(PyExc_ValueError,
											"speed must be a 1D to 12-D array of doubles");
			Py_XDECREF(phi);
			Py_XDECREF(dx);
			Py_XDECREF(flag);
			return nullptr;
		}

		if (! PyArray_SAMESHAPE(phi,speed))
		{
			PyErr_SetString(PyExc_ValueError,
											"phi and speed must have the same shape");
			Py_XDECREF(phi);
			Py_XDECREF(dx);
			Py_XDECREF(flag);
			Py_XDECREF(speed);
			return nullptr;
		}
	}

  if (mode == TRAVEL_TIME_GENES) {
    int speeds_dim = PyArray_NDIM(phi) + 1;
    speeds = (PyArrayObject *)PyArray_FROMANY(pspeeds, NPY_DOUBLE, 0, 0, NPY_IN_ARRAY);
    std::printf("pspeeds = %x\n", pspeeds);                                          
    std::printf("speeds = %x\n", speeds);                                          
    if (!speeds)
    {
      PyErr_SetString(PyExc_ValueError,
                      "speeds not initialised");
      Py_XDECREF(phi);
      Py_XDECREF(dx);
      Py_XDECREF(flag);
      return nullptr;
    }

    PyArrayObject *tmp = (PyArrayObject *)pdrivers;
    printf("drivers dtype = %d\n", PyArray_TYPE(tmp));
    drivers = (PyArrayObject *)PyArray_FROMANY(pdrivers, NPY_DOUBLE, 1, 12, NPY_IN_ARRAY);

    std::printf("pdrivers = %x\n", pdrivers);                                          
    std::printf("drivers = %x\n", drivers);                                          

    if (!drivers) {
      PyErr_SetString(PyExc_ValueError,
                      "drivers not initialised");
      Py_XDECREF(phi);
      Py_XDECREF(dx);
      Py_XDECREF(flag);
      Py_XDECREF(speeds);
      return nullptr;
    }

    if (! PyArray_SAMESHAPE(phi,drivers))
    {
      // print shapes of phi and drivers
      std::printf("drivers size %d\n", drivers->nd);
      std::printf("drivers array ndim %d\n", PyArray_SHAPE(drivers));
      std::printf("phi array ndim %d\n", PyArray_SHAPE(phi));
      std::printf("%d\n", Py_TYPE(pdrivers));
      std::printf("%d\n", Py_TYPE(drivers));
      std::printf("%d\n", Py_TYPE(phi));
      PyErr_SetString(PyExc_ValueError,
                      "phi and drivers must have the same shape");
      Py_XDECREF(phi);
      Py_XDECREF(dx);
      Py_XDECREF(flag);
      Py_XDECREF(speeds);
      Py_XDECREF(drivers);
      return nullptr;
    }
  }

  if (! (PyArray_NDIM(phi)==(npy_intp)PyArray_DIM(dx,0))) // ?!
  {
    PyErr_SetString(PyExc_ValueError, "dx must be of length len(phi.shape)");
    Py_XDECREF(phi);
    Py_XDECREF(dx);
    Py_XDECREF(flag);
    Py_XDECREF(speed);
    return nullptr;
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
      return nullptr;
    }
  }

  if (! PyArray_SAMESHAPE(phi,flag))
  {
    PyErr_SetString(PyExc_ValueError, "phi and flag must have the same shape");
    Py_XDECREF(phi);
    Py_XDECREF(dx);
    Py_XDECREF(flag);
    Py_XDECREF(speed);
    return nullptr;
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
                                            shape2, NPY_DOUBLE, 0);
  if (!distance) return nullptr;

  if (mode == EXTENSION_VELOCITY)
  {
    f_ext = (PyArrayObject *)PyArray_ZEROS(PyArray_NDIM(phi),
                                           shape2, NPY_DOUBLE, 0);
    if (!f_ext) return nullptr;

    ext_mask = (PyArrayObject *)PyArray_FROMANY(pext_mask, NPY_LONGLONG, 1,
                                                12, NPY_IN_ARRAY);
    if (!ext_mask)
    {
      PyErr_SetString(PyExc_ValueError,
                      "ext_mask must be a 1D, 12-D array of integers");
      Py_XDECREF(phi);
      Py_XDECREF(dx);
      Py_XDECREF(flag);
      Py_XDECREF(speed);
      return nullptr;
    }
  }

  // create a level set object to do the calculation
  double* local_phi        = (double *) PyArray_DATA(phi);
  double* local_dx         = (double *) PyArray_DATA(dx);
  long long   *local_flag       = (long long *)   PyArray_DATA(flag);
  long long   *local_ext_mask   = nullptr;
  if (ext_mask) local_ext_mask = (long long *) PyArray_DATA(ext_mask);
  double* local_speed      = nullptr;
  if (speed) local_speed    = (double* ) PyArray_DATA(speed);
  //speeds and drivers used for genetics extension
  unsigned* local_drivers      = nullptr;
  if (drivers) local_drivers    = (unsigned *)PyArray_DATA(drivers);
  double* local_speeds    = nullptr;
  if (speeds) local_speeds    = (double *)PyArray_DATA(speeds);
  double* local_distance   = (double *)PyArray_DATA(distance);
  int error;

  baseMarcher* marcher = nullptr;
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
      double*  local_fext = (double *) PyArray_DATA(f_ext);
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
    case TRAVEL_TIME_GENES:
    {
      marcher = new travelTimeMarcherGenes(
        local_phi,
        local_dx,
        local_flag,
        local_distance,
        PyArray_NDIM(phi),
        shape,
        self_test,
        order,
        local_drivers,
        local_speeds,
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
      return nullptr;
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
    return nullptr;
  case 2:
    PyErr_SetString(PyExc_ValueError,
                    "the array phi contains no zero contour (no zero level set)");
    Py_XDECREF(distance);
    return nullptr;
  }

  if (mode == EXTENSION_VELOCITY)
  {
    return Py_BuildValue("NN", distance, f_ext);
  }
  return (PyObject *)distance;
}
