Place for some 4fun gr sims i have and might do.

geodesic-gr directory has code for simulating particle movements in Schwarzschild spacetime around a black hole.
The particles are randomly generated. Simulation settings are in setting.txt

How it works:
Sagemath does symbolic calculation to get the geodesic equation for desired metric (schwarzschild metric is the only case that has been implemented currently) and then generated c++ header file with the functions based on the generated geodesic equations. The main c++ then uses these equations in the header files to simulate the particles.
Python code is for data generation and visualization.

Original idea was to be able to input any metric in sagemath then generate geodesic equation and be able to simulate based on the metric. First step was making the current version. Not interested in continuing this project right now, since just making this was quite satisfactory and I achieved my goal of learning bit of sagemath and doing some code generation using it.