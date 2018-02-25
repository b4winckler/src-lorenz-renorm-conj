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


void error(const char *msg)
{
    std::cerr << __FILE__ << ':' << __LINE__ << ' ' << msg << std::endl;
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
    int d;

    lorenz_map() : d(2) { }
    lorenz_map(const vec3<scalar> &parameters_)  : parameters(parameters_), d(2) {}

    int alpha() const { return d; }
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
            y = 1 - pow(1 - y, 2);
            y = v0() + y * (1 - v0());
        } else if (x > c()) {
            // Right branch
            y = (x - c()) / (1 - c());
            y = pow(y, 2);
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
            y = 1 - sqrt(1 - y);
            y = y * c();
        } else {
            y = x / v1();
            y = sqrt(y);
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

void pull_back_lorenz(
        // output
        mpreal &y,
        // input
        const mpreal &x, const mpreal &c, const mpreal &v0, const mpreal &v1,
        int left)
{
    if (left) {
        y = (x - v0) / (1 - v0);
        y = 1 - sqrt(1 - y);
        y = y * c;
    } else {
        y = x / v1;
        y = sqrt(y);
        y = c + (1 - c) * y;
    }
}

// Realize a (w0,w1)-renormalizable Lorenz map using the Thurston algorithm.
// The resulting map has the property that its parameters (v0, v1) are fixed
// under renormalization.
template <typename scalar>
void realize_renormalizable_map(
        lorenz_map<scalar> &f,
        const std::string &w0, const std::string &w1, const scalar &c,
        bool verbose)
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
    for (; sqr_err > max_sqr_err && count < 1e6; ++count) {
        // This is the Thurston pull back step: pull back each point in such a
        // way that a fixed point of the pull back is an actual orbit of the
        // desired Lorenz family.  The magic of this algorithm is that it does
        // converge to a fixed point.

        // *If* x0/x1 was the orbit of 0/1, then the second element is v0/v1
        v0 = x0[1];
        v1 = x1[1];

        // This ensures a fixed point of the Thurston algorithm fixes v0 and v1
        const scalar &l = x0[n1 - 1], &r = x1[n0 - 1];
        x0[n - 1] = l + v0 * (r - l);
        x1[n - 1] = l + v1 * (r - l);

        // debug_vector(x0);
        // debug_vector(x1);

        sqr_err = 0;
        for (size_t i = 1; i < n; ++i) {
            tmp = x0[i - 1];
            pull_back_lorenz(x0[i - 1], x0[i], c, v0, v1, x0[i - 1] < c);
            sqr_err += sqr(tmp - x0[i - 1]);

            tmp = x1[i - 1];
            pull_back_lorenz(x1[i - 1], x1[i], c, v0, v1, x1[i - 1] < c);
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

    renormalization_operator(const std::string &w0_, const std::string &w1_)
        : w0(w0_), w1(w1_) { }

    template <typename t>
    void operator() (const vec3<t> &x, vec3<t> *py) const
    {
        lorenz_map<t> f(x);

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
        const std::string &w0, const std::string &w1, bool verbose)
{
    renormalization_operator<scalar> op(w0, w1);
    AutoDiffJacobian< renormalization_operator<scalar> > renormalize(op);
    vec3<scalar> y, h;
    mat3<scalar> jacobian;
    scalar sqr_err(1);
    size_t count = 0;

    for (; sqr_err > max_sqr_err && count < 1e6; ++count) {
        // Determine v0 and v1
        realize_renormalizable_map(f, w0, w1, f.c(), false);

        // Take a Newton step
        renormalize(f.parameters, &y, &jacobian);
        auto lu = (jacobian - mat3<scalar>::Identity()).fullPivLu();
        h = lu.solve(y - f.parameters);
        f.parameters[0] -= h[0];

        sqr_err = h[0] * h[0];
    }

    if (verbose) {
        std::cerr << "newton: #iterations = " << count << ", error = " <<
            sqrt(sqr_err) << std::endl;
    }
}

void usage()
{
    std::cerr << "usage: lorenz w0 w1 c\n\n";
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
    if (argc != 4)
        usage();

    mpreal::set_default_prec(precision);

    std::string w0(argv[1]);
    std::string w1(argv[2]);
    mpreal c(atof(argv[3]));

    if (!is_admissible(w0, w1)) {
        std::cerr << "inadmissible combinatorics" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (!(c > 0 && c < 1)) {
        std::cerr << "critical point must lie in (0,1)" << std::endl;
        exit(EXIT_FAILURE);
    }

#if 0
    std::vector<mpreal> x;
    realize_orbit_with_itinerary(x, w0 + w1, c);
    debug_vector(x);
#endif

    // Realize a (w0,w1)-renormalizable map with critical point c
    lorenz_map<mpreal> f0;
    realize_renormalizable_map(f0, w0, w1, c, true);
    std::cerr << "f = " << f0.c() << ' ' << f0.v0() << ' ' << f0.v1() << std::endl;
    bool b = is_quasi_renormalizable(f0, w0, w1);
    std::cerr << "quasi-renormalizable: " << b << std::endl;

    vec3<mpreal> y;
    mat3<mpreal> jac;
    renormalization_operator<mpreal> op(w0, w1);
    AutoDiffJacobian< renormalization_operator<mpreal> > renormalize(op);
    renormalize(f0.parameters, &y, &jac);

    std::cerr << "R(f) =\n\t" << y.transpose() << std::endl;
    std::cerr << "spec DR(f) =\n\t" << jac.eigenvalues().transpose() <<
        std::endl;

    realize_fixed_point(f0, w0, w1, true);
    std::cerr << "fixed point f0 =\n\t" << f0.parameters.transpose() << std::endl;
    renormalize(f0.parameters, &y, &jac);
    std::cerr << "R(f0) =\n\t" << y.transpose() << std::endl;
    std::cerr << "spec DR(f0) =\n\t" << jac.eigenvalues().transpose() <<
        std::endl;


    return EXIT_SUCCESS;
}
