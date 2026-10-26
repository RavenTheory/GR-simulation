#include <array>
#include <cmath>
#include <cstdio>
#include "geodesic.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;


using Vec4 = std::array<double, 4>;

struct State {Vec4 x, v; };




State derivative(const State& s) {
    return {s.v, geodesic_acc(s.x, s.v)};
}


State slope_args(const State& s, const State& k, double h) {
    State arg;
    for (int i=0; i<4; ++i) {
        arg.x[i] = s.x[i] + h * k.x[i];
        arg.v[i] = s.v[i] + h * k.v[i];
    }
    return arg;
}


State RK4(const State& s, double h) {
    State k1, k2, k3, k4;
    k1 = derivative(s);
    k2 = derivative(slope_args(s, k1, 0.5 * h));
    k3 = derivative(slope_args(s, k2, 0.5 * h));
    k4 = derivative(slope_args(s, k3, h));
    State result;
    for (int i = 0; i < 4; ++i) {
        result.x[i] = s.x[i] + (h/6.0) * (k1.x[i] + 2*k2.x[i] + 2*k3.x[i] + k4.x[i]);
        result.v[i] = s.v[i] + (h/6.0) * (k1.v[i] + 2*k2.v[i] + 2*k3.v[i] + k4.v[i]);
    }
    return result;
}


void write_frame(std::ofstream& file, std::size_t frame, const std::vector<State>& particles) {
    for (const auto& p : particles) {
        file << frame << ';'
             << p.x[0] << ';'
             << p.x[1] << ';'
             << p.x[2] << ';'
             << p.x[3] << '\n';
    }
}

bool is_infinite_or_nan(const State& s) {
    for (const auto& val : s.x) {
        if (std::isinf(val) || std::isnan(val)) {
            return true;
        }
    }
    for (const auto& val : s.v) {
        if (std::isinf(val) || std::isnan(val)) {
            return true;
        }
    }
    return false;
}

bool has_escaped(const State& s, double escape_radius = 30.0) {
    return s.x[1] > escape_radius;
}

void simulate(vector<State> states, double h, int steps, std::ofstream& file, double escape_radius= 30.0) {

    file << "frame;t;r;theta;phi\n";

    for (auto i=0; i<steps; ++i) {

        for(std::size_t j=0; j< states.size();) {
            states[j] = RK4(states[j], h);
            if (states[j].x[1] <= 2.01) { // event horizon at r=2
                states[j] = std::move(states.back());
                states.pop_back();
                continue;
            }
            if (is_infinite_or_nan(states[j])){
                states[j] = std::move(states.back());
                states.pop_back();
                continue;
            }
            if (has_escaped(states[j], escape_radius)) {
                states[j] = std::move(states.back());
                states.pop_back();
                continue;
            }
            ++j;
        }
        write_frame(file, i + 1, states);
    }
    file.close();
}

void read_settings(map<string, string>& settings) {
    string line;
    ifstream ReadFile("settings.txt");
    while (getline (ReadFile, line)) {
        if (line[0] == '#') { //skip line if it's commented
            continue;
        }
        string variable = "";
        string value = "";
        bool variable_name_area = true;
        for (int i=0; i<line.length(); i++) {
            if (line[i] == ' ' || line[i] == '\n') {
                continue;
            }
            if (line[i] == '=') {
                variable_name_area = false;
                continue;
            }
            if (variable_name_area) {
                variable += line[i];
            } else {
                value += line[i];
            }
        }
        settings[variable] = value;
    }
}

std::vector<State> read_initial_states(const std::string& filename)
{
    std::ifstream file(filename);

    if (!file)
        throw std::runtime_error("Could not open " + filename);

    std::vector<State> states;

    std::string line;

    // Skip header
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;

        State state;

        /*// id
        std::getline(ss, value, ',');
        state.id = std::stoull(value);
        */

        // r
        std::getline(ss, value, ',');
        state.x[1] = std::stod(value);

        // theta
        std::getline(ss, value, ',');
        state.x[2] = std::stod(value);

        // phi
        std::getline(ss, value, ',');
        state.x[3] = std::stod(value);

        // u1
        std::getline(ss, value, ',');
        state.v[1] = std::stod(value);

        // u2
        std::getline(ss, value, ',');
        state.v[2] = std::stod(value);

        // u3
        std::getline(ss, value, ',');
        state.v[3] = std::stod(value);

        // Initial coordinate time
        state.x[0] = 0.0;

        // Calculate this properly from g(u,u) = -1
        state.v[0] = timelike_velocity(state.x, {state.v[1], state.v[2], state.v[3]})[0];

        states.push_back(state);
    }

    return states;
}


int main() {
    map<string, string> settings;
    read_settings(settings);

    double M = std::stod(settings["M"]); // Mass of the black hole
    double h = std::stod(settings["h"]); // Step size
    int steps = std::stoi(settings["steps"]); // Number of steps
    double escape_radius = std::stod(settings["escape_radius"]); // Escape radius
    string initial_data_file = settings["data_file_name"]; // Initial data file

    std::vector<State> initial_states = read_initial_states(initial_data_file);

    std::ofstream output("trajectory.csv");

    simulate(initial_states, h, steps, output);



    return 0;
}















