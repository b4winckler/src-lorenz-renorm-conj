#define USE_MPFR 0

#include <iostream>

#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <unsupported/Eigen/AutoDiff>

#if USE_MPFR
#include <unsupported/Eigen/MPRealSupport>
#define real mpreal
using namespace mpfr;
#else
#define real double
#endif


#define error(msg) \
{ \
    std::cerr << __FILE__ << ':' << __LINE__ << ": " << (msg) << std::endl; \
    exit(EXIT_FAILURE); \
} \


using namespace Eigen;


template <typename t> using vec = Matrix<t, Dynamic, 1>;
template <typename t> using mat = Matrix<t, Dynamic, Dynamic>;


const size_t diffeo_size = 10;
const size_t lorenz_size = 3 + 2 * diffeo_size;


template <typename scalar>
void pure(scalar &y, const scalar &x, const scalar &s)
{
    scalar r = exp(s);
    y = x * (2 + (r - 1) * x) / (r + 1);
}

template <typename scalar>
void clear_lorenz(vec<scalar> &lorenz, scalar c = 0.5, scalar v0 = 0, scalar v1 = 1)
{
    lorenz.resize(lorenz_size);
    lorenz.setZero(lorenz_size);
    lorenz[0] = c;
    lorenz[1] = v0;
    lorenz[2] = v1;
}

template <typename scalar>
const scalar &crit(const vec<scalar> &lorenz)
{
    return lorenz[0];
}

template <typename scalar>
void apply(scalar &y, const scalar &x, const vec<scalar> &lorenz, const real &alpha)
{
    if (x < lorenz[0]) {
        // Left branch
        y = x / lorenz[0];
        y = 1 - pow(1 - y, alpha);

        for (size_t k = 0; k < diffeo_size; ++k)
            pure(y, y, lorenz[3 + k]);

        y = lorenz[1] + y * (1 - lorenz[1]);
    } else if (x > lorenz[0]) {
        // Right branch
        y = (x - lorenz[0]) / (1 - lorenz[0]);
        y = pow(y, alpha);

        for (size_t k = 0; k < diffeo_size; ++k)
            pure(y, y, lorenz[3 + diffeo_size + k]);

        y = y * lorenz[2];
    } else {
        error("cannot apply lorenz map");
    }
}

template <typename scalar>
void iterate(scalar &y, const scalar &x, size_t n, const vec<scalar> &lorenz, const real &alpha)
{
    y = x;
    for (size_t i = 0; i < n; ++i)
        apply(y, y, lorenz, alpha);
}

template <typename scalar>
struct thurston_op {
    typedef vec<scalar> InputType;
    typedef vec<scalar> ValueType;

    vec<real> &family2d;
    const real &alpha;
    const std::string &w0;
    const std::string &w1;

    thurston_op(vec<scalar> &family2d_, const real &alpha_,
            const std::string &w0_, const std::string &w1_)
        : family2d(family2d_), alpha(alpha_), w0(w0_), w1(w1_) { }

    void guess(vec<scalar> &x)
    {
    }

    template <typename t>
    void operator() (const vec<t> &x, vec<t> *py) const
    {
    }
};

template <typename scalar>
struct renorm_op {
    typedef vec<scalar> InputType;
    typedef vec<scalar> ValueType;

    const real &alpha;
    const std::string &w0;
    const std::string &w1;

    renorm_op(const real &alpha_, const std::string &w0_,
            const std::string &w1_)
        : alpha(alpha_), w0(w0_), w1(w1_) { }

    template <typename t>
    void operator() (const vec<t> &x, vec<t> *py) const
    {
        size_t n0 = w0.size(), n1 = w1.size();

        // l = f^{n1 - 1)(0), vl = f^{n0}(l)
        t l(0), vl;
        iterate(l, l, n1 - 1, x, alpha);
        iterate(vl, l, n0, x, alpha);

        // r = f^{n0 - 1)(1), vr = f^{n1}(r)
        t r(1), vr;
        iterate(r, r, n0 - 1, x, alpha);
        iterate(vr, r, n1, x, alpha);

        // The return interval is [l, r] and the image of the first-return map
        // of f on this interval is [vl, vr]
        vec<t> &y = *py;
        y[0] = (crit(x) - l) / (r - l);
        y[1] = (vl - l) / (r - l);
        y[2] = (vr - l) / (r - l);
    }
};

int main(int argc, char *argv[])
{
    if (argc != 5) {
        std::cerr << "usage: renorm w0 w1 c alpha\n\n";
        exit(EXIT_FAILURE);
    }

    std::string w0(argv[1]);
    std::string w1(argv[2]);
#if USE_MPFR
    real c(argv[3]);
    real alpha(argv[4]);
#else
    real c = atof(argv[3]);
    real alpha = atof(argv[4]);
#endif

    vec<real> f0, f1;
    renorm_op<real> renorm(alpha, w0, w1);
    AutoDiffJacobian< renorm_op<real> > drenorm(renorm);
    mat<real> jac(lorenz_size, lorenz_size);

    clear_lorenz(f0, c);
    clear_lorenz(f1);
    drenorm(f0, &f1, &jac);

    std::cerr << f0 << std::endl;
    std::cerr << f1 << std::endl;
    std::cerr << jac << std::endl;

    return EXIT_SUCCESS;
}
