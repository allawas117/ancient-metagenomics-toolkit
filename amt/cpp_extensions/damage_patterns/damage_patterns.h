#ifndef DAMAGE_PATTERNS_H
#define DAMAGE_PATTERNS_H

#include <Python.h>

static PyObject* analyze_damage_wrapper(PyObject* self, PyObject* args);

PyMODINIT_FUNC PyInit_damage_patterns(void);

#endif // DAMAGE_PATTERNS_H