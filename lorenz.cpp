#include <unsupported/Eigen/MPRealSupport>
#include <unsupported/Eigen/AutoDiff>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <cassert>

using namespace mpfr;
using namespace Eigen;

template <typename t>
using vec3 = Matrix<t, 3, 1>;

template <typename t>
using vecN = Matrix<t, Dynamic, 1>;

template <typename t>
using mat3 = Matrix<t, 3, 3>;

static mp_prec_t precision = 512;
static mpreal max_sqr_err("1e-200");

// Examples of combinatorics:
//
// once (8,2)-renormalizable
// LRRRRRRRR RLL
// twice (8,2)-renormalizable
// LRRRRRRRRRLLRLLRLLRLLRLLRLLRLLRLL RLLLRRRRRRRRLRRRRRRRR
//
// once (6,2)-renormalizale
// LRRRRRR RLL
// twice (6,2)-renormalizale
// LRRRRRRRLLRLLRLLRLLRLLRLL RLLLRRRRRRLRRRRRR


void error(const char *msg)
{
    std::cerr << __FILE__ << ':' << __LINE__ << ": " << msg << std::endl;
    exit(EXIT_FAILURE);
}

// A Lorenz map f is defined on the unit interval, minus the critical point c.
// It is normalized so that the critical values are 0 and 1; i.e. f(c-) = 1 and
// f(c+) = 0.
// The boundary values v0 = f(0) and v1 = f(1) are parameters (and so is c).
// The critical exponent at c is denoted alpha (and is typically fixed).
//
// f is (w0,w1)-renormalizable iff:
//
//  *   the first-return map to [l, r] is (affinely conjugate) a Lorenz map,
//  *   the itinerary of the left/right branch of the first-return map is
//      w0/w1,
//  *   l = f^{n1 - 1}(0) and r = f^{n0 - 1}(1), where n0/n1 is the length of the
//      word w0/w1.
//

template <typename scalar>
struct lorenz_map {
    vec3<scalar> parameters;
    mpreal d;

    lorenz_map(const mpreal &alpha_) : d(alpha_) { }
    lorenz_map(const vec3<scalar> &parameters_, const mpreal &alpha_)
        : parameters(parameters_), d(alpha_) {}

    mpreal alpha() const { return d; }
    const scalar &c() const { return parameters[0]; }
    const scalar &v0() const { return parameters[1]; }
    const scalar &v1() const { return parameters[2]; }
    void set_parameters(const scalar &c, const scalar &v0, const scalar &v1)
    {
        parameters[0] = c; parameters[1] = v0; parameters[2] = v1;
    }

    // Compute f(x)
    void apply(scalar &y, const scalar &x) const
    {
        if (x < c()) {
            // Left branch
            y = x / c();
            y = 1 - pow(1 - y, alpha());
            y = v0() + y * (1 - v0());
        } else if (x > c()) {
            // Right branch
            y = (x - c()) / (1 - c());
            y = pow(y, alpha());
            y = y * v1();
        } else {
            error("attempting to evaluate f at c");
        }
    }

    // Compute f^{-1}(x)
    void pull_back(scalar &y, const scalar &x, int left_branch) const
    {
        if (left_branch) {
            y = (x - v0()) / (1 - v0());
            y = 1 - pow(1 - y, 1 / alpha());
            y = y * c();
        } else {
            y = x / v1();
            y = pow(y, 1 / alpha());
            y = c() + (1 - c()) * y;
        }
    }

    // Compute f^n(x)
    void iterate(scalar &y, const scalar &x, size_t n) const
    {
        y = x;
        for (size_t i = 0; i < n; ++i)
            apply(y, y);
    }

    // Compute x, f(x), ..., f^n(x)
    void orbit(vecN<scalar> &y, const scalar &x, size_t n) const
    {
        if (y.size() < n)
            y.resize(n);

        y[0] = x;
        for (size_t i = 1; i < n; ++i)
            apply(y[i], y[i - 1]);
    }
};

void debug_vector(const std::vector<mpreal> &x)
{
    for (auto t : x)
        std::cerr << t << ' ';
    std::cerr << std::endl;
}

bool is_admissible(const std::string &w0, const std::string &w1)
{
    if (!(w0.size() > 0 && 'L' == w0[0] && w1.size() > 0 && 'R' == w1[0]))
        return false;

    for (auto s : w0)
        if (!(s == 'L' || 'R' == s))
            return false;

    for (auto s : w1)
        if (!(s == 'L' || 'R' == s))
            return false;

    return true;
}

// Return orbit of the full affine Lorenz map (with critical point c) whose
// itinerary is given by w.  The returned orbit is a pre-orbit of c, so
//
//      f(x[i]) = x[i + 1]  and  f(x[n - 1]) = c,
//
// where n is the length of the word w.
void realize_orbit_with_itinerary(
        std::vector<mpreal> &x,
        const std::string &w, const mpreal &c)
{
    size_t n = w.size();
    x.resize(n);

    // Pull back n times with the full affine Lorenz map
    for (mpreal y(c); n-- > 0; y = x[n])
        x[n] = w[n] == 'L' ? c - (1 - y) * c : c + y * (1 - c);
}

// Realize a (w0,w1)-renormalizable Lorenz map using the Thurston algorithm.
// The resulting map has the property that its parameters (v0, v1) are fixed
// under renormalization.
template <typename scalar>
void realize_renormalizable_map(
        lorenz_map<scalar> &f,
        const std::string &w0, const std::string &w1, const scalar &c,
        const mpreal &sqr_eps, bool verbose)
{
    // Find admissible initial guess to Thurston algorithm
    // TODO: why is it correct to choose pre-orbit of full affine map?
    std::vector<scalar> x0, x1;
    realize_orbit_with_itinerary(x0, (w1 + w0).substr(1), c);
    realize_orbit_with_itinerary(x1, (w0 + w1).substr(1), c);
    x0.push_back(c);
    x1.push_back(c);

    size_t n0 = w0.size(), n1 = w1.size(), n = n0 + n1;
    scalar sqr_err(1), tmp, v0, v1;
    size_t count = 0;
    for (; sqr_err > sqr_eps && count < 1e6; ++count) {
        // This is the Thurston pull back step: pull back each point in such a
        // way that a fixed point of the pull back is an actual orbit of the
        // desired Lorenz family.  The magic of this algorithm is that it does
        // converge to a fixed point.

        // *If* x0/x1 was the orbit of 0/1, then the second element is v0/v1
        v0 = x0[1];
        v1 = x1[1];
        f.set_parameters(c, v0, v1);

        // This ensures a fixed point of the Thurston algorithm fixes v0 and v1
        const scalar &l = x0[n1 - 1], &r = x1[n0 - 1];
        x0[n - 1] = l + v0 * (r - l);
        x1[n - 1] = l + v1 * (r - l);

        // debug_vector(x0);
        // debug_vector(x1);

        sqr_err = 0;
        for (size_t i = 1; i < n; ++i) {
            tmp = x0[i - 1];
            f.pull_back(x0[i - 1], x0[i], x0[i - 1] < c);
            sqr_err += sqr(tmp - x0[i - 1]);

            tmp = x1[i - 1];
            f.pull_back(x1[i - 1], x1[i], x1[i - 1] < c);
            sqr_err += sqr(tmp - x1[i - 1]);
        }
    }

    f.set_parameters(c, x0[1], x1[1]);
    if (verbose) {
        std::cerr << "thurston: #iterations = " << count << ", error = " <<
            sqrt(sqr_err) << std::endl;
        // debug_vector(x0);
        // debug_vector(x1);
    }
}

// Check if the map which sends x[i] to y[i] is increasing.
bool is_increasing(const std::vector<mpreal> &x, const std::vector<mpreal> &y)
{
    assert(x.size() == y.size());

    // Determine permutation which sorts x
    std::vector<size_t> perm(x.size());
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(),
            [&](size_t i, size_t j){ return x[i] < x[j]; });

    // Check if y is ordered the same way as x
    for (size_t i = 1; i < x.size(); ++i)
        if (y[perm[i - 1]] > y[perm[i]])
            return false;

    return true;
}

// Compute the orbit x0, f(x0), .., f^{n-1}(x0)
template <typename scalar>
void orbit(
        // output
        std::vector<mpreal> &x,
        // input
        const mpreal &x0, size_t n, const lorenz_map<scalar> &f)
{
    x.resize(n);
    x[0] = x0;

    for (size_t i = 1; i < n; ++i) {
        if (x[i - 1] < f.c()) {
            // Left branch
            x[i] = x[i - 1] / f.c();
            x[i] = 1 - pow(1 - x[i], f.alpha());
            x[i] = f.v0() + x[i] * (1 - f.v0());
        } else if (x[i - 1] > f.c()) {
            // Right branch
            x[i] = (x[i - 1] - f.c()) / (1 - f.c());
            x[i] = pow(x[i], f.alpha());
            x[i] = x[i] * f.v1();
        } else {
            error("attempting to evaluate f at c");
        }
    }
}

// Check that the critical orbits are ordered in the correct way for f to be
// (w0,w1)-quasi-renormalizable.
template <typename scalar>
bool is_quasi_renormalizable(
        const lorenz_map<scalar> &f, const std::string &w0, const std::string &w1)
{
    size_t n0 = w0.size(), n1 = w1.size();
    if (!(n0 > 0 && n1 > 0))
        return true;
    if (!('L' == w0[0] && 'R' == w1[0]))
        return false;

    // Realize (w0,w1)-quasi-renormalizable itineraries
    size_t n = n0 + n1 - 1;
    std::vector<scalar> y0, y1;
    realize_orbit_with_itinerary(y0, (w1 + w0).substr(1, n), f.c());
    realize_orbit_with_itinerary(y1, (w0 + w1).substr(1, n), f.c());
    y0.insert(y0.end(), y1.begin(), y1.end());

    // Determine the orbits of the critical values 0 and 1 under f
    std::vector<scalar> x0(n), x1(n);
    orbit(x0, 0, n, f);
    orbit(x1, 1, n, f);
    x0.insert(x0.end(), x1.begin(), x1.end());

    // f is quasi-renormalizable iff x0 and y0 are ordered the same way
    return is_increasing(x0, y0);
}

template <typename scalar>
struct renormalization_operator {
    typedef vec3<scalar> InputType;
    typedef vec3<scalar> ValueType;

    std::string w0, w1;
    mpreal alpha;

    renormalization_operator(const std::string &w0_, const std::string &w1_,
            const mpreal &alpha_)
        : w0(w0_), w1(w1_), alpha(alpha_) { }

    template <typename t>
    void operator() (const vec3<t> &x, vec3<t> *py) const
    {
        lorenz_map<t> f(x, alpha);

        size_t n0 = w0.size(), n1 = w1.size();
        assert(n0 > 0 && n1 > 0);

        // l = f^{n1 - 1)(0), vl = f^{n0}(l)
        t l(0), vl;
        f.iterate(l, l, n1 - 1);
        f.iterate(vl, l, n0);

        // r = f^{n0 - 1)(1), vr = f^{n1}(r)
        t r(1), vr;
        f.iterate(r, r, n0 - 1);
        f.iterate(vr, r, n1);

        // The return interval is [l, r] and the image of the first-return map
        // of f on this interval is [vl, vr]
        *py <<
            (f.c() - l) / (r - l),
            (vl - l) / (r - l),
            (vr - l) / (r - l);
    }
};

// Find a fixed point of the (w0,w1)-renormalization operator.
// The issue lies in determining c.  To improve chances of convergence we
// combine the Thurston algorithm (to determine v0, v1) and a Newton iteration
// (to determine c).
template <typename scalar>
void realize_fixed_point(lorenz_map<scalar> &f,
        const std::string &w0, const std::string &w1, const scalar &c0,
        bool verbose)
{
    renormalization_operator<scalar> op(w0, w1, f.alpha());
    AutoDiffJacobian< renormalization_operator<scalar> > renormalize(op);
    vec3<scalar> y, h;
    mat3<scalar> jacobian;
    scalar sqr_err(1), last_sqr_err(2), c = c0;
    size_t count = 0;

    while (count++ < 1e6) {
        // Determine v0 and v1
        realize_renormalizable_map(f, w0, w1, c, max_sqr_err, false);

        // Calculate Newton step
        renormalize(f.parameters, &y, &jacobian);
        auto lu = (jacobian - mat3<scalar>::Identity()).fullPivLu();
        h = lu.solve(y - f.parameters);
        // std::cerr << "x = " << f.parameters.transpose() << "\ny = " <<
        //     y.transpose() << "\nh = " << h.transpose() << std::endl;

        last_sqr_err = sqr_err;
        sqr_err = h[0] * h[0];
        if (sqr_err > last_sqr_err && sqr_err < max_sqr_err) {
            sqr_err = last_sqr_err;
            break;  // error is small but increasing, so stop
        }

        c -= h[0];

        if (c < 0 || c > 1)
            error("newton step moved critical point outside (0,1)");
    }

    f.parameters[0] = c;

    if (verbose) {
        std::cerr << "newton+thurston: #iterations = " << count << ", error = " <<
            sqrt(sqr_err) << std::endl;
    }
}

// Find a fixed point of the (w0,w1)-renormalization operator using Newton
// iteration
template <typename scalar>
void realize_fixed_point_newton(lorenz_map<scalar> &f,
        const std::string &w0, const std::string &w1,
        const lorenz_map<scalar> f0, bool verbose)
{
    renormalization_operator<scalar> op(w0, w1, f.alpha());
    AutoDiffJacobian< renormalization_operator<scalar> > renormalize(op);
    vec3<scalar> y, h;
    mat3<scalar> jacobian;
    scalar sqr_err(1), last_sqr_err(2);
    size_t count = 0;

    f = f0;
    if (!is_quasi_renormalizable(f, w0, w1)) {
        if (verbose) std::cerr << "starting guess not quasi-renormalizable\n";
        realize_renormalizable_map(f, w0, w1, f.c(), max_sqr_err, false);
    }

    while (count++ < 1e6) {
        // Calculate Newton step
        renormalize(f.parameters, &y, &jacobian);
        auto lu = (jacobian - mat3<scalar>::Identity()).fullPivLu();
        h = lu.solve(y - f.parameters);

        last_sqr_err = sqr_err;
        sqr_err = h.squaredNorm();
        if (sqr_err > last_sqr_err && sqr_err < max_sqr_err) {
            sqr_err = last_sqr_err;
            break;  // error is small but increasing, so stop
        }

        f.parameters -= h;

        if (f.c() < 0 || f.c() > 1)
            error("newton step moved critical point outside (0,1)");
        if (f.v0() < 0 || f.v0() > 1)
            error("newton step moved v0 outside (0,1)");
        if (f.v1() < 0 || f.v1() > 1)
            error("newton step moved v1 outside (0,1)");
    }

    if (verbose) {
        std::cerr << "newton: #iterations = " << count << ", error = " <<
            sqrt(sqr_err) << std::endl;
    }
}

template <typename scalar>
void realize_renormalizable_map_fast(
        lorenz_map<scalar> &f,
        const std::string &w0, const std::string &w1, const scalar &c,
        bool verbose)
{
    // Find initial guess using Thurston iteration
    mpreal sqr_eps("1e-6");
    realize_renormalizable_map(f, w0, w1, c, sqr_eps, verbose);

    renormalization_operator<scalar> op(w0, w1, f.alpha());
    AutoDiffJacobian< renormalization_operator<scalar> > renormalize(op);
    vec3<scalar> y, h;
    mat3<scalar> jacobian;
    scalar sqr_err(1), last_sqr_err(2);
    size_t count = 0;

    while (count++ < 1e6) {
        // Calculate Newton step in (v0,v1)-direction
        renormalize(f.parameters, &y, &jacobian);
        jacobian(0, 0) = 1;
        jacobian(0, 1) = 0;
        jacobian(0, 2) = 0;
        jacobian(1, 0) = 0;
        jacobian(2, 0) = 0;
        auto lu = (jacobian - mat3<scalar>::Identity()).fullPivLu();
        h = lu.solve(y - f.parameters);

        last_sqr_err = sqr_err;
        sqr_err = h[1] * h[1] + h[2] * h[2];
        if (sqr_err > last_sqr_err && sqr_err < max_sqr_err) {
            sqr_err = last_sqr_err;
            break;  // error is small but increasing, so stop
        }

        f.parameters[1] -= h[1];
        f.parameters[2] -= h[2];

        if (f.v0() < 0 || f.v0() > 1)
            error("newton step moved v0 outside (0,1)");
        if (f.v1() < 0 || f.v1() > 1)
            error("newton step moved v1 outside (0,1)");
    }

    if (verbose) {
        std::cerr << "thurston->newton: #iterations = " << count << ", error = " <<
            sqrt(sqr_err) << std::endl;
    }
}

void print_fixed_point(
        const std::string &w0, const std::string &w1, const mpreal &c,
        const mpreal &alpha)
{
    // Realize a (w0,w1)-renormalizable map with critical point c
    lorenz_map<mpreal> f0(alpha);
    realize_renormalizable_map(f0, w0, w1, c, max_sqr_err, true);
    std::cerr << "f = " << f0.c() << ' ' << f0.v0() << ' ' << f0.v1() <<
        std::endl;
    bool b = is_quasi_renormalizable(f0, w0, w1);
    std::cerr << "quasi-renormalizable: " << b << std::endl;

    vec3<mpreal> y;
    mat3<mpreal> jac;
    renormalization_operator<mpreal> op(w0, w1, f0.alpha());
    AutoDiffJacobian< renormalization_operator<mpreal> > renormalize(op);
    renormalize(f0.parameters, &y, &jac);

    std::cerr << "R(f) =\n\t" << y.transpose() << std::endl;
    std::cerr << "spec DR(f) =\n\t" << jac.eigenvalues().transpose() <<
        std::endl;

    lorenz_map<mpreal> f(alpha);
    realize_fixed_point(f, w0, w1, c, true);
    std::cerr << "fixed point f =\n\t" << f.parameters.transpose() << std::endl;
    renormalize(f.parameters, &y, &jac);
    std::cerr << "R(f) =\n\t" << y.transpose() << std::endl;
    std::cerr << "spec DR(f) =\n\t" << jac.eigenvalues().transpose() <<
        std::endl;

    realize_fixed_point_newton(f, w0, w1, f0, true);
    std::cerr << "fixed point f =\n\t" << f.parameters.transpose() << std::endl;
    renormalize(f.parameters, &y, &jac);
    std::cerr << "R(f) =\n\t" << y.transpose() << std::endl;
    std::cerr << "spec DR(f) =\n\t" << jac.eigenvalues().transpose() <<
        std::endl;
}

void print_monotone_type_fixed_points()
{
    std::vector<mpreal> alphas = { 1.25, 1.5, 1.75, 2, 2.5, 3, 4, 5, 10, 20 };
    std::string w0 = "LR";

    std::cout << "#w0\tw1\talpha\tc\tv0\tv1\teval0\teval1\teval2" << std::endl;

    for (;;) {
        for (const mpreal &alpha : alphas) {
            lorenz_map<mpreal> f(alpha);
            f.set_parameters(0.5, 0, 1);
            std::string w1(w0.size(), 'L');
            w1[0] = 'R';

            while (w1.size() > 1) {
                realize_fixed_point(f, w0, w1, f.c(), true);

                renormalization_operator<mpreal> op(w0, w1, f.alpha());
                AutoDiffJacobian< renormalization_operator<mpreal> > renormalize(op);
                vec3<mpreal> y;
                mat3<mpreal> jac;
                renormalize(f.parameters, &y, &jac);
                auto evals = jac.eigenvalues();

                std::cout <<
                    w0 << '\t' <<
                    w1 << '\t' <<
                    f.alpha().toString() << '\t' <<
                    f.c().toString() << '\t' <<
                    f.v0().toString() << '\t' <<
                    f.v1().toString() << '\t' <<
                    evals[0].real().toString() << '\t' <<
                    evals[1].real().toString() << '\t' <<
                    evals[2].real().toString() <<
                    std::endl;

                w1.pop_back();
            }
        }

        w0.push_back('R');
    }
}

// Evaluate (w0,w1)-renormalization on words (rw0, rw1).
//
// Use case: if an actual Lorenz map f is twice (w0,w1)-renormalizable then its
// second (w0,w1)-renormalization equals its first (rw0,rw1)-renormalization,
// where (rw0,rw1) is the (w0,w1)-renormalization of (w0,w1).
void renormalize_words(
        std::string &rw0, std::string &rw1,
        const std::string &w0, const std::string &w1)
{
    std::string rw0_copy(rw0);
    rw0.clear();
    for (auto ch : rw0_copy)
        rw0.append('L' == ch ? w0 : w1);

    std::string rw1_copy(rw1);
    rw1.clear();
    for (auto ch : rw1_copy)
        rw1.append('L' == ch ? w0 : w1);
}

// For a range of critical points c, find the map which fixes (v0,v1) and print
// the critical point of its 'nrenorm' first renormalizations.
// The output is one line of critical points, followed by one line each for the
// critical points of the n-th renormalization.
void print_critical_point_movement(
        const std::string &w0, const std::string &w1, const mpreal &alpha)
{
    static int nrenorm = 2;
    static int ngrid = 100;

    std::string rw0("L"), rw1("R");
    lorenz_map<mpreal> f(alpha);
    vec3<mpreal> y;
    mpreal c;
    mpreal scale(ngrid + 1);
    scale = 1 / scale;

    // Print c coordinates
    for (int j = 0; j < ngrid; ++j) {
        c = (j + 1) * scale;
        std::cout << c.toString() << (j == ngrid - 1 ? '\n' : '\t');
    }

    // Print c(R^{i+1}(f)) coordinates
    for (int i = 0; i < nrenorm; ++i) {
        renormalize_words(rw0, rw1, w0, w1);

        for (int j = 0; j < ngrid; ++j) {
            c = (j + 1) * scale;
            realize_renormalizable_map(f, rw0, rw1, c, max_sqr_err, false);
            renormalization_operator<mpreal> renormalize(rw0, rw1, f.alpha());
            renormalize(f.parameters, &y);
            std::cout << y[0].toString() << (j == ngrid - 1 ? '\n' : '\t');
        }
    }
}

int main(int argc, char *argv[])
{
    mpreal::set_default_prec(precision);

#if 0
    if (argc != 4) {
        std::cerr << "usage: lorenz w0 w1 alpha\n\n";
        exit(EXIT_FAILURE);
    }

    std::string w0(argv[1]);
    std::string w1(argv[2]);
    mpreal alpha(argv[3]);
    print_critical_point_movement(w0, w1, alpha);
#elif 0
    if (argc != 1) {
        std::cerr << "usage: lorenz\n\n";
        exit(EXIT_FAILURE);
    }

    print_monotone_type_fixed_points();
#elif 0
    if (argc != 5) {
        std::cerr << "usage: lorenz w0 w1 c alpha\n\n";
        exit(EXIT_FAILURE);
    }

    std::string w0(argv[1]);
    std::string w1(argv[2]);
    mpreal c(argv[3]);
    mpreal alpha(argv[4]);

    if (!is_admissible(w0, w1)) {
        std::cerr << "inadmissible combinatorics" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (!(c > 0 && c < 1)) {
        std::cerr << "critical point must lie in (0,1)" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (!(alpha > 1)) {
        std::cerr << "critical exponent must be strictly larger than 1" << std::endl;
        exit(EXIT_FAILURE);
    }

    print_fixed_point(w0, w1, c, alpha);
#elif 1
    if (argc != 5) {
        std::cerr << "usage: lorenz w0 w1 c alpha\n\n";
        exit(EXIT_FAILURE);
    }

    std::string w0(argv[1]);
    std::string w1(argv[2]);
    mpreal c(argv[3]);
    mpreal alpha(argv[4]);

    lorenz_map<mpreal> f(alpha);
    // realize_renormalizable_map(f, w0, w1, c, max_sqr_err, true);
    realize_renormalizable_map_fast(f, w0, w1, c, true);

    std::cerr << "fixed point f =\n\t" << f.parameters.transpose() << std::endl;
#endif

    return EXIT_SUCCESS;
}
