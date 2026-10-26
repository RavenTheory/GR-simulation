#pragma once
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

    const double x0 = 2*M;
    const double x1 = -r + x0;
    const double x2 = 1.0/x1;
    const double x3 = 1.0/r;
    const double x4 = v1*x3;
    const double x5 = std::pow(v3, 2);
    const double x6 = std::sin(theta);
    const double x7 = std::cos(theta);
    a[0] = v0*x0*x2*x4;
    a[1] = -M*std::pow(v1, 2)*x2*x3 + M*std::pow(v0, 2)*x1/std::pow(r, 3) - std::pow(v2, 2)*x1 - x1*x5*std::pow(x6, 2);
    a[2] = -2*v2*x4 + x5*x6*x7;
    a[3] = -2*v3*(v2*x7/x6 + x4);

    return a;
}

inline std::array<double, 4> timelike_velocity(const std::array<double, 4>& x, const std::array<double, 3>& three_velocity, [[maybe_unused]] double M = 1.0){
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
    u[3] = u3;
    u[0] = -std::sqrt(-2*M*std::pow(r, 3)*std::pow(u2, 2) - 2*M*std::pow(r, 3)*std::pow(u3, 2)*std::pow(std::sin(theta), 2) - 2*M*r + std::pow(r, 4)*std::pow(u2, 2) + std::pow(r, 4)*std::pow(u3, 2)*std::pow(std::sin(theta), 2) + std::pow(r, 2)*std::pow(u1, 2) + std::pow(r, 2))/(2*M - r);
    return u;
    }
    

inline std::array<double, 4> spacelike_velocity(const std::array<double, 4>& x, const std::array<double, 3>& three_velocity, [[maybe_unused]] double M = 1.0){
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
    u[3] = u3;
    u[0] = -std::sqrt(-2*M*std::pow(r, 3)*std::pow(u2, 2) - 2*M*std::pow(r, 3)*std::pow(u3, 2)*std::pow(std::sin(theta), 2) + 2*M*r + std::pow(r, 4)*std::pow(u2, 2) + std::pow(r, 4)*std::pow(u3, 2)*std::pow(std::sin(theta), 2) + std::pow(r, 2)*std::pow(u1, 2) - std::pow(r, 2))/(2*M - r);
    return u;
    }
    

inline std::array<double, 4> null_velocity(const std::array<double, 4>& x, const std::array<double, 3>& three_velocity, [[maybe_unused]] double M = 1.0){
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
    u[3] = u3;
    u[0] = -r*std::sqrt(-2*M*r*std::pow(u2, 2) - 2*M*r*std::pow(u3, 2)*std::pow(std::sin(theta), 2) + std::pow(r, 2)*std::pow(u2, 2) + std::pow(r, 2)*std::pow(u3, 2)*std::pow(std::sin(theta), 2) + std::pow(u1, 2))/(2*M - r);
    return u;
    }
    
