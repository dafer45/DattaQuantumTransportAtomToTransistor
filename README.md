# DattaQuantumTransportAtomToTransistor

This project aims to provide implementations of the exercises in the book "Quantum Transport: Atom to Transistor" by S. Datta (2005), using the C++ library [TBTK](https://github.com/dafer45/TBTK/) for second quantized models.
The project is meant to demonstrate the ease with which a wide range of quantum mechanical calculations can be performed using TBTK.
It also aims to be an entry point for students interested in quantum transport.
Aiding such students by providing an implementation that puts focus on higher level concepts rather than irrelevant numerical details.

Each exercise is implemented in a separate folder called E_X_X, where X_X is the exercise number.
To compile and the exercises, the [TBTK](https://github.com/dafer45/TBTK/) library must first be installed.
Once TBTK is installed, the exercises can be compiled and run as follows
```bash
cd E_X_X
cmake .
make
./build/Application
```

## Currently completed exercises
Exercise 2.1
