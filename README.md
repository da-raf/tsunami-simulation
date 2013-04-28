tsunami-simulation
==================

This project will solve riemann problems to simulate a tsunami.
It was developed in an university course at the TUM (Technische Universität München)

If you want to use this solver, just include 'FWave.hpp' from the 'src' directory.
The only function you need is computeNetUpdates (see doxygen docs).


----------------------
DOCUMENTATION

The documentation can be built as follows:

doxygen doxygen.conf

This will give you html documentation inside the 'html' directory.
If you want a pdf:

cd latex
make

Now open refman.pdf


----------------------
UNIT TESTING

The unit tests can be built and executed as follows:

cd src
cxxtestgen --error-printer FWaveSolver.hpp -o runner.cpp
g++ runner.cpp -o runner
./runner
