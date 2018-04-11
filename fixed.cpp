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
using vec2 = Matrix<t, 2, 1>;
template <typename t>
using vec3 = Matrix<t, 3, 1>;

template <typename t>
using vecN = Matrix<t, Dynamic, 1>;

template <typename t>
using mat2 = Matrix<t, 2, 2>;
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


#define error(msg) \
{ \
    std::cerr << __FILE__ << ':' << __LINE__ << ": " << (msg) << std::endl; \
    exit(EXIT_FAILURE); \
} \

// A Lorenz map f is defined on the unit interval, minus the critical point c.
// It is normalized so that the critical values are 0 and 1; i.e. f(c-) = 1 and
// f(c+) = 0.
// The boundary values v0 = f(0) and v1 = f(1) are parameters (and so is c).
// The critical exponent at c is denoted alpha (and is typically fixed).
//
// f is (w0,w1)-renormalizable iff:
//
//  *   the first-return map to [l, r] is (affinely conjugate to) a Lorenz map,
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

void thurston_guess(
        std::vector<mpreal> &x0, std::vector<mpreal> &x1,
        const std::string &w0, const std::string &w1, const mpreal &c)
{
    // Find admissible initial guess to Thurston algorithm
    // TODO: why is it correct to choose pre-orbit of full affine map?
    realize_orbit_with_itinerary(x0, (w1 + w0).substr(1), c);
    realize_orbit_with_itinerary(x1, (w0 + w1).substr(1), c);
    x0.push_back(c);
    x1.push_back(c);
}

void thurston_iterate(
        std::vector<mpreal> &x0, std::vector<mpreal> &x1,
        lorenz_map<mpreal> &f,
        const std::string &w0, const std::string &w1,
        const mpreal &sqr_eps, bool verbose)
{
    size_t n0 = w0.size(), n1 = w1.size(), n = n0 + n1;
    mpreal sqr_err(1), tmp, v0, v1;
    size_t count = 0;

    for (; sqr_err > sqr_eps && count < 1e6; ++count) {
        // *If* x0/x1 was the orbit of 0/1, then the second element is v0/v1
        v0 = x0[1];
        v1 = x1[1];
        f.parameters[1] = v0;
        f.parameters[2] = v1;

        // This ensures a fixed point of the Thurston algorithm fixes v0 and v1
        const mpreal &l = x0[n1 - 1], &r = x1[n0 - 1];
        x0[n - 1] = l + v0 * (r - l);
        x1[n - 1] = l + v1 * (r - l);

        // debug_vector(x0);
        // debug_vector(x1);

        sqr_err = 0;
        for (size_t i = 1; i < n; ++i) {
            tmp = x0[i - 1];
            f.pull_back(x0[i - 1], x0[i], x0[i - 1] < f.c());
            sqr_err += sqr(tmp - x0[i - 1]);

            tmp = x1[i - 1];
            f.pull_back(x1[i - 1], x1[i], x1[i - 1] < f.c());
            sqr_err += sqr(tmp - x1[i - 1]);
        }
    }

    f.parameters[1] = x0[1];
    f.parameters[2] = x1[1];

    if (verbose) {
        std::cerr << "thurston: #iterations = " << count << ", error = " <<
            sqrt(sqr_err) << std::endl;
        // debug_vector(x0);
        // debug_vector(x1);
    }
}

// Realize a (w0,w1)-renormalizable Lorenz map using the Thurston algorithm.
// The resulting map has the property that its parameters (v0, v1) are fixed
// under renormalization.
template <typename scalar>
void thurston(
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

template <typename scalar>
struct rop {
    typedef vec3<scalar> InputType;
    typedef vec3<scalar> ValueType;

    std::string w0, w1;
    mpreal alpha;

    rop(const std::string &w0_, const std::string &w1_, const mpreal &alpha_)
        : w0(w0_), w1(w1_), alpha(alpha_)
    {
        assert(w0.size() > 0 && w1.size() > 0);
    }

    template <typename t>
    void operator() (const vec3<t> &x, vec3<t> *py) const
    {
        lorenz_map<t> f(x, alpha);
        size_t n0 = w0.size(), n1 = w1.size();

        // l = f^{n1 - 1)(0), vl = f^{n0}(l)
        t l(0), vl;
        f.iterate(l, l, n1 - 1);
        f.iterate(vl, l, n0);

        // r = f^{n0 - 1)(1), vr = f^{n1}(r)
        t r(1), vr;
        f.iterate(r, r, n0 - 1);
        f.iterate(vr, r, n1);

        *py <<
            f.c() * (r - l - 1) + l,
            f.v0() * (r - l) - vl + l,
            f.v1() * (r - l) - vr + l;
    }
};

template <typename scalar>
void newton(lorenz_map<scalar> &f,
        const std::string &w0, const std::string &w1,
        const lorenz_map<scalar> f0, bool verbose)
{
    rop<scalar> op(w0, w1, f.alpha());
    AutoDiffJacobian< rop<scalar> > dop(op);
    vec3<scalar> y, h;
    mat3<scalar> jac;
    scalar sqr_err(1), last_sqr_err(2);
    size_t count = 0;

    f = f0;

    while (count++ < 1e3) {
        // Calculate Newton step
        dop(f.parameters, &y, &jac);
        h = jac.fullPivLu().solve(y);

        last_sqr_err = sqr_err;
        sqr_err = h.squaredNorm();
        if (sqr_err >= last_sqr_err && sqr_err < max_sqr_err) {
            sqr_err = last_sqr_err;
            break;  // error is small but increasing, so stop
        }

        std::cerr << f.parameters.transpose() << '\t' << h.transpose() <<
            '\t' << jac.eigenvalues().transpose() << std::endl;

        f.parameters -= h;

        if (f.c() < 0 || f.c() > 1)
            error("newton step moved critical point outside (0,1)");
        if (f.v0() < 0 || f.v0() > 1)
            error("newton step moved v0 outside (0,1)");
        if (f.v1() < 0 || f.v1() > 1)
            error("newton step moved v1 outside (0,1)");
    }

    if (count > 100)
        std::cerr << "newton took unusually long to coverge: " << count <<
            std::endl;

    if (verbose) {
        std::cerr << "newton: #iterations = " << count << ", error = " <<
            sqrt(sqr_err) << std::endl;
        std::cerr << "        eigenvalues = " << jac.eigenvalues().transpose()
            << std::endl;
    }
}

template <typename scalar>
void newton_thurston(lorenz_map<scalar> &f,
        const std::string &w0, const std::string &w1,
        const lorenz_map<scalar> f0, bool verbose)
{
    rop<scalar> op(w0, w1, f0.alpha());
    AutoDiffJacobian< rop<scalar> > dop(op);
    vec3<scalar> y, h;
    mat3<scalar> jac;
    scalar sqr_err(1), last_sqr_err(2);
    size_t count = 0;

    f = f0;

    int stage = 0;
    while (count++ < 1e3) {
        if (0 == stage)
            thurston(f, w0, w1, f.c(), mpreal("1e-50"), verbose);

        dop(f.parameters, &y, &jac);
        h = jac.fullPivLu().solve(y);

        last_sqr_err = sqr_err;
        sqr_err = 0 == stage ? h[0] * h[0] : h.squaredNorm();
        if (sqr_err >= last_sqr_err && sqr_err < max_sqr_err) {
            sqr_err = last_sqr_err;
            if (0 == stage) {
                if (verbose)
                    std::cerr << "newton+thurston: "
                        "entering second stage #iterations = " << count
                        << std::endl;
                stage = 1;
                continue;
            } else {
                break;
            }
        }

        if (verbose)
            std::cerr << "newton+thurston: step = " << h.transpose() <<
                std::endl;

        if (0 == stage)
            f.parameters[0] -= h[0];
        else
            f.parameters -= h;

        if (f.c() < 0 || f.c() > 1)
            error("newton step moved critical point outside (0,1)");
        if (f.v0() < 0 || f.v0() > 1)
            error("newton step moved v0 outside (0,1)");
        if (f.v1() < 0 || f.v1() > 1)
            error("newton step moved v1 outside (0,1)");
    }

    if (count > 100)
        std::cerr << "newton+thurston: took unusually long to coverge: " <<
            count << std::endl;

    if (verbose) {
        std::cerr << "newton+thurston: #iterations = " << count <<
            ", error = " << sqrt(sqr_err) << std::endl;
        std::cerr << "                 eigenvalues = "
            << jac.eigenvalues().transpose() << std::endl;
    }
}

int main(int argc, char *argv[])
{
    mpreal::set_default_prec(precision);

    if (argc != 5) {
        std::cerr << "usage: fixed w0 w1 c alpha\n\n";
        exit(EXIT_FAILURE);
    }

    std::string w0(argv[1]);
    std::string w1(argv[2]);
    mpreal c(argv[3]);
    mpreal alpha(argv[4]);

    // Realize a (w0,w1)-renormalizable map with critical point c
    lorenz_map<mpreal> f0(alpha);
    f0.parameters[0] = c;
    // thurston(f0, w0, w1, c, max_sqr_err, true);
    // std::cerr << "f0 = " << f0.c() << ' ' << f0.v0() << ' ' << f0.v1() <<
    //     std::endl;
    lorenz_map<mpreal> f = f0;
    // newton(f, w0, w1, f0, true);

    newton_thurston(f, w0, w1, f0, true);
    std::cerr << "R fixed point =\n" << f.c().toString() << '\n' <<
        f.v0().toString() << '\n' << f.v1().toString() <<
        std::endl;

    renormalization_operator<mpreal> op(w0, w1, f0.alpha());
    AutoDiffJacobian< renormalization_operator<mpreal> > renorm(op);
    vec3<mpreal> y;
    mat3<mpreal> jac;
    renorm(f.parameters, &y, &jac);
    std::cout << "DR eigenvalues =\n" << jac.eigenvalues();

    return EXIT_SUCCESS;
}
