// nbody.cpp  — simple student-style N-body
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <random>
#include <iomanip>

using namespace std;

struct Vec3 {
    double x=0, y=0, z=0;
};

struct Particle {
    double m=1.0;
    Vec3 pos, vel, force;
};

struct Config {
    double G = 6.674e-11;   // gravitational constant
    double dt = 1.0;        // time step
    int steps = 1000;       // number of iterations
    int dump_every = 10;    // how often to print a line
    double softening = 1e-3; // small value to avoid div-by-0 (unit-dependent)
    unsigned seed = 42;     // RNG for random init
};

struct NBodyState {
    vector<Particle> p;

    size_t size() const { return p.size(); }
    void resize(size_t n) { p.assign(n, Particle{}); }

    // --- init: random ---
    void init_random(size_t N, double mass_min, double mass_max,
                     double pos_span, double vel_span, unsigned seed=42) {
        resize(N);
        mt19937 rng(seed);
        uniform_real_distribution<double> Um(mass_min, mass_max);
        uniform_real_distribution<double> U(-1.0, 1.0);
        for (auto &b : p) {
            b.m = Um(rng);
            b.pos = { U(rng)*pos_span, U(rng)*pos_span, U(rng)*pos_span };
            b.vel = { U(rng)*vel_span, U(rng)*vel_span, U(rng)*vel_span };
            b.force = {0,0,0};
        }
    }

    // --- init: Sun–Earth–Moon (rough numbers) ---
    void init_sem() {
        resize(3);
        // masses (kg)
        p[0].m = 1.98847e30; // Sun
        p[1].m = 5.9722e24;  // Earth
        p[2].m = 7.342e22;   // Moon
        // positions (m)
        p[0].pos = {0,0,0};
        p[1].pos = {1.496e11, 0, 0};
        p[2].pos = {1.496e11 + 3.844e8, 0, 0};
        // velocities (m/s) approx tangential
        p[0].vel = {0,0,0};
        p[1].vel = {0, 29780, 0};
        p[2].vel = {0, 29780 + 1022, 0};
        for (auto &b : p) b.force = {0,0,0};
    }


    bool load_from_file(const std::string& path) {
        std::ifstream fin(path);
        if (!fin) return false;

        // read the first non-empty, non-comment line
        std::string line;
        while (std::getline(fin, line)) {
            // trim leading spaces
            std::size_t i = 0;
            while (i < line.size() && std::isspace((unsigned char)line[i])) ++i;
            if (i >= line.size() || line[i] == '#') continue; // skip blank or comment
            line.erase(0, i);
            break;
        }
        if (line.empty()) return false;

        std::istringstream iss(line);

        int N;
        if (!(iss >> N) || N <= 0) return false;
        resize((std::size_t)N);

        for (int i = 0; i < N; ++i) {
            Particle b;
            if (!(iss >> b.m
                      >> b.pos.x >> b.pos.y >> b.pos.z
                      >> b.vel.x >> b.vel.y >> b.vel.z
                      >> b.force.x >> b.force.y >> b.force.z)) {
                return false; // not enough numbers on the line
                      }
            p[i] = b;
        }
        return true;
    }

};

// --- math helpers ---
static inline Vec3 add(const Vec3&a,const Vec3&b){ return {a.x+b.x,a.y+b.y,a.z+b.z}; }
static inline Vec3 sub(const Vec3&a,const Vec3&b){ return {a.x-b.x,a.y-b.y,a.z-b.z}; }
static inline Vec3 mul(const Vec3&a,double s){ return {a.x*s,a.y*s,a.z*s}; }

static inline double norm2(const Vec3& v){ return v.x*v.x + v.y*v.y + v.z*v.z; }

// --- compute forces (O(N^2)) ---
void compute_forces(NBodyState& S, const Config& cfg) {
    // reset
    for (auto &b : S.p) b.force = {0,0,0};

    const double eps2 = cfg.softening * cfg.softening;

    int N = (int)S.size();
    for (int i=0;i<N;i++) {
        for (int j=i+1;j<N;j++) {
            Vec3 rij = sub(S.p[j].pos, S.p[i].pos);
            double r2 = norm2(rij) + eps2;
            double r  = sqrt(r2);
            if (r == 0) continue; // avoid divide-by-zero (identical positions)
            // |F| = G m1 m2 / r^2
            double Fmag = cfg.G * S.p[i].m * S.p[j].m / r2;
            // unit vector r_hat = rij / r
            Vec3 rhat = mul(rij, 1.0/r);
            Vec3 Fij  = mul(rhat, Fmag);
            // apply equal and opposite
            S.p[i].force = add(S.p[i].force, Fij);
            S.p[j].force = sub(S.p[j].force, Fij);
        }
    }
}

// --- integrate one step (Euler: v += a*dt; x += v*dt) ---
void step_euler(NBodyState& S, const Config& cfg) {
    for (auto &b : S.p) {
        Vec3 a = { b.force.x / b.m, b.force.y / b.m, b.force.z / b.m };
        b.vel = add(b.vel, mul(a, cfg.dt));
        b.pos = add(b.pos, mul(b.vel, cfg.dt)); // uses v_new per spec
    }
}

// --- write one TSV line: N then each particle's (m, x y z, vx vy vz, fx fy fz) ---
// Print positions centered on the Sun and scaled to gigameters (Gm = 1e9 m).
// Velocities left in m/s (or switch to km/s by dividing by 1e3 below).
void write_state_tsv(const NBodyState& S, std::ostream& out) {
    const double POS_SCALE = 1e9; // 1e9 m => values ~ 0..150 for Earth-Sun (fits chart)
    // const double VEL_SCALE = 1e3; // uncomment to print km/s instead of m/s

    // Center on the Sun (particle 0) so scene doesn't drift off frame
    Vec3 c = S.p[0].pos;

    out << S.size();
    out << std::setprecision(12);
    for (const auto& b : S.p) {
        double px = (b.pos.x - c.x) / POS_SCALE;
        double py = (b.pos.y - c.y) / POS_SCALE;
        double pz = (b.pos.z - c.z) / POS_SCALE;

        // If you want km/s, use (b.vel.x / VEL_SCALE). Keeping m/s is fine for plot.py.
        out << '\t' << b.m
            << '\t' << px << '\t' << py << '\t' << pz
            << '\t' << b.vel.x << '\t' << b.vel.y << '\t' << b.vel.z
            << '\t' << b.force.x << '\t' << b.force.y << '\t' << b.force.z;
    }
    out << '\n';
}


// --- main: parse args and run ---
int main(int argc, char** argv) {
    ios::sync_with_stdio(false);

    string mode = argv[1];
    Config cfg;
    cfg.dt = stod(argv[2]);
    cfg.steps = stoi(argv[3]);
    cfg.dump_every = stoi(argv[4]);


    NBodyState S;

    // decide init mode
    auto is_number = [](const string& s)->bool{
        if (s.empty()) return false;
        char* end=nullptr;
        strtod(s.c_str(), &end);
        // is it integer-ish?
        bool all_digit = all_of(s.begin(), s.end(), [](char c){ return isdigit((unsigned char)c); });
        return all_digit && end && *end=='\0';
    };

    if (mode == "sem") {
        S.init_sem();
    } else if (is_number(mode)) {
        size_t N = (size_t)stoull(mode);
        // simple random ranges: pick something sane (units arbitrary)
        S.init_random(N,
              1e22,   // mass_min: ~small dwarf planet
              1e27,   // mass_max: ~gas giant like Jupiter
              1.0e11, // pos_span: ~1 AU (distance Earth–Sun)
              1.0e4,  // vel_span: ~10 km/s (planet orbital speeds)
              cfg.seed);

    } else {
        if (!S.load_from_file(mode)) {
            cerr << "failed to load file: " << mode << "\n";
            return 2;
        }
    }

    // initial forces + initial dump
    compute_forces(S, cfg);
    write_state_tsv(S, cout);

    for (int step=1; step<=cfg.steps; ++step) {
        // 1) forces at current positions
        compute_forces(S, cfg);
        // 2) integrate
        step_euler(S, cfg);
        // 3) output occasionally
        if (step % cfg.dump_every == 0) {
            // recompute forces for logging (forces correspond to printed state)
            compute_forces(S, cfg);
            write_state_tsv(S, cout);
        }
    }
    return 0;
}
