import sympy as sp
from sympy.printing.cxx import cxxcode

n = 4

L = Manifold(n, 'L', latexname=r'\mathcal{L}', start_index=0, structure='Lorentzian')


#Different coodrinate charts but having multiple with transit between them can make sagemath
#die cuz it starts to calculate stuff for every chart lmao

#X.<t,x,y,z>=L.chart()
Spherical.<t,r,theta,phi>=L.chart(r't r:(0,+oo) theta:(0,pi):\theta phi:(0,2*pi):\phi')

#transit_Spherical_to_X = Spherical.transition_map(X, [t, r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta)])
#transit_Spherical_to_X.set_inverse(t, sqrt(x^2+y^2+z^2), atan2(sqrt(x^2+y^2),z), atan2(y, x))



#L.set_default_chart(Spherical)
#L.set_default_frame(Spherical.frame())


#Metric
g = L.metric('g')




#Metric functions
#Schwarzschild metric in natural units:
#f = (1 - r_s/r) = 1- 2M/r
#We assume M=1.0 for most calculations currently

var('M')
f = Spherical.function(1-2*M/r)


g[0,0], g[1,1], g[2,2], g[3,3] = -f, 1/f, r^2, r^2 * sin(theta)^2



christoffel_symbols = g.christoffel_symbols()
v = sp.symbols('v0:4')

acceleration = [sp.Integer(0)]*n

iterator = christoffel_symbols.non_redundant_index_generator()

for i,j,k in iterator  :
    term = - christoffel_symbols[i, j, k].expr()._sympy_() * v[j] * v[k]
    if j!= k:
        term *= 2
    acceleration[i] += term

repl, reduced = sp.cse(acceleration, optimizations='basic')


for k in repl:
    print(k)

for j in reduced:
    print(j)


lines = []

for temp_var, expr in repl:
    lines.append(f"    const double {temp_var} = {cxxcode(expr)};")

for i, expr in enumerate(reduced):
    lines.append(f"    a[{i}] = {cxxcode(expr)};")

#for line in lines:
#    print(line)

prologue = """#pragma once
#include <array>
#include <cmath>

inline std::array<double,4> geodesic_acc(const std::array<double, 4>& x, const std::array<double, 4>& v, [[maybe_unused]] double M = 1.0){

    [[maybe_unused]] const double t = x[0];
    [[maybe_unused]] const double r = x[1];
    [[maybe_unused]] const double theta = x[2];
    [[maybe_unused]] const double phi = x[3];

    const double v0 = v[0];
    const double v1 = v[1];
    const double v2 = v[2];
    const double v3 = v[3];

    std::array<double, 4> a{};

"""

epilogue = """
    return a;
}
"""




#Velocity field
var('u0 u1 u2 u3')
u = L.vector_field(u0, u1, u2, u3, name='u')
var('timelike', 'spacelike', 'null')
timelike = solve(g(u,u).expr() == -1, u0)[1]
spacelike = solve(g(u,u).expr() == 1, u0)[1]
null = solve(g(u,u).expr() == 0, u0)[1]

timelike_rhs = cxxcode(timelike.rhs()._sympy_())
spacelike_rhs = cxxcode(spacelike.rhs()._sympy_())
null_rhs = cxxcode(null.rhs()._sympy_())


def vel_function(velocityrhs, signature):
    writing = f"inline std::array<double, 4> {signature}(const std::array<double, 4>& x, const std::array<double, 3>& three_velocity, [[maybe_unused]] double M = 1.0)" + """{
    [[maybe_unused]] const double t = x[0];
    [[maybe_unused]] const double r = x[1];
    [[maybe_unused]] const double theta = x[2];
    [[maybe_unused]] const double phi = x[3];

    const double u1 = three_velocity[0];
    const double u2 = three_velocity[1];
    const double u3 = three_velocity[2];

    std::array<double, 4> u{};
    u[1] = u1;
    u[2] = u2;
    u[3] = u3;\n""" + f"    u[0] = {velocityrhs};" + """
    return u;
    }
    """
    return writing







with open("geodesic.hpp", "w") as f:
    f.write(prologue + "\n".join(lines) + "\n" + epilogue)
    f.write("\n" + vel_function(timelike_rhs, "timelike_velocity") + "\n")
    f.write("\n" + vel_function(spacelike_rhs, "spacelike_velocity") + "\n")
    f.write("\n" + vel_function(null_rhs, "null_velocity") + "\n")














