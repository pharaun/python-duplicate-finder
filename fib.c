#include <Python.h>

static PyObject * fib(PyObject *self, PyObject *args) {
    int a = 0;
    int b = 1;
    int c; 
    int n;

    if (!PyArg_ParseTuple(args, "i", &n))
	return NULL;

    PyObject *list = PyList_New(0);
    PyObject *number;

    while(b < n){
	number = PyInt_FromLong(b);
	PyList_Append(list, number);
	Py_DECREF(number);

	c = a+b;
	a = b;
	b = c;
    }

    return list;
}

PyMethodDef methods[] = {
    {"fib", fib, METH_VARARGS, "Returns a fibonacci sequence as a list"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initfib() {
    (void) Py_InitModule("fib", methods);   
}
