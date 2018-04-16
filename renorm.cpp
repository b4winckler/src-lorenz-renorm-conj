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


const size_t diffeo_size = 5;

struct context {
    real alpha;
    vec<real> grid0;
    vec<real> grid1;
};


template <typename scalar>
void pure(scalar &y, const scalar &x, const scalar &s)
{
    scalar r = exp(s);
    y = x * (2 + (r - 1) * x) / (r + 1);
}

template <typename scalar>
void init_lorenz(vec<scalar> &lorenz, const context &ctx, scalar c = 0.5,
        scalar v0 = 0, scalar v1 = 1)
{
    size_t n = 3 + ctx.grid0.size() + ctx.grid1.size();
    lorenz.resize(n);
    lorenz[0] = c;
    lorenz[1] = v0;
    lorenz[2] = v1;
    lorenz.segment(3, ctx.grid0.size()) = ctx.grid0;
    lorenz.tail(ctx.grid1.size()) = ctx.grid1;
}

template <typename scalar>
const scalar &crit(const vec<scalar> &lorenz)
{
    return lorenz[0];
}

template <typename scalar>
void apply(scalar &y, const scalar &x, const vec<scalar> &lorenz,
        const context &ctx)
{
    // Adjust for critical point, then fold
    size_t i0, k0, k1;
    bool left = x < lorenz[0];
    if (left) {
        y = x / lorenz[0];
        y = 1 - pow(1 - y, ctx.alpha);
        i0 = k0 = 3;
        k1 = k0 + ctx.grid0.size() - 1;
    } else if (x > lorenz[0]) {
        y = (x - lorenz[0]) / (1 - lorenz[0]);
        y = pow(y, ctx.alpha);
        i0 = k0 = 3 + ctx.grid0.size();
        k1 = k0 + ctx.grid1.size() - 1;
    }

    // Interpolate diffeomorphism
    if (y > 0 && y < 1) {
        const vec<real> &grid = left ? ctx.grid0 : ctx.grid1;
        size_t k;
        // With autodiff library it is not possible (?) to take the floor of y
        // so we need to search to find the interval on which to interpolate.
        while (k1 > k0 + 1) {
            k = (k0 + k1) / 2;
            if (y < grid[k - i0]) {
                k1 = k;
            } else if (y > grid[k - i0]) {
                k0 = k;
            } else {
                k0 = k;
                k1 = k0 + 1;
            }
        }
        y = (grid[k1 - i0] - y) * lorenz[k0] + (y - grid[k0 - i0]) * lorenz[k1];
        y /= (grid[k1 - i0] - grid[k0 - i0]);
    }

    // Adjust for boundary value
    y = left ? lorenz[1] + y * (1 - lorenz[1]) : y = y * lorenz[2];
}

template <typename scalar>
void iterate(scalar &y, const scalar &x, size_t n, const vec<scalar> &lorenz, const context &ctx)
{
    y = x;
    for (size_t i = 0; i < n; ++i)
        apply(y, y, lorenz, ctx);
}

template <typename scalar>
struct thurston_op {
    typedef vec<scalar> InputType;
    typedef vec<scalar> ValueType;

    vec<real> &family2d;
    const context &ctx;
    const std::string &w0;
    const std::string &w1;

    thurston_op(vec<scalar> &family2d_, const context &ctx_,
            const std::string &w0_, const std::string &w1_)
        : family2d(family2d_), ctx(ctx_), w0(w0_), w1(w1_) { }

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

    const context &ctx;
    const std::string &w0;
    const std::string &w1;

    renorm_op(const context &ctx_, const std::string &w0_,
            const std::string &w1_)
        : ctx(ctx_), w0(w0_), w1(w1_) { }

    template <typename t>
    void operator() (const vec<t> &x, vec<t> *py) const
    {
        size_t n0 = w0.size(), n1 = w1.size();

        // l = f^{n1 - 1)(0), vl = f^{n0}(l)
        t l(0), vl;
        iterate(l, l, n1 - 1, x, ctx);
        iterate(vl, l, n0, x, ctx);

        // r = f^{n0 - 1)(1), vr = f^{n1}(r)
        t r(1), vr;
        iterate(r, r, n0 - 1, x, ctx);
        iterate(vr, r, n1, x, ctx);

        // The return interval is [l, r] and the image of the first-return map
        // of f on this interval is [vl, vr]
        vec<t> &y = *py;
        y[0] = (crit(x) - l) / (r - l);
        y[1] = (vl - l) / (r - l);
        y[2] = (vr - l) / (r - l);

        // Endpoints of diffeos
        y[3] = 0;
        y[3 + ctx.grid0.size()] = 0;
        y[3 + ctx.grid0.size() - 1] = 1;
        y[3 + ctx.grid0.size() + ctx.grid1.size() - 1] = 1;

        for (size_t k = 1; k < ctx.grid0.size() - 1; ++k) {
            t p(l + ctx.grid0[k] * (crit(x) - l));
            iterate(p, p, n0, x, ctx);
            y[3 + k] = (p - vl) / (r - vl);
        }

        for (size_t k = 1; k < ctx.grid1.size() - 1; ++k) {
            t p(crit(x) + ctx.grid0[k] * (r - crit(x)));
            iterate(p, p, n1, x, ctx);
            y[3 + ctx.grid0.size() + k] = (p - l) / (vr - l);
        }
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

    context ctx;
    ctx.alpha = alpha;
    ctx.grid0.setLinSpaced(diffeo_size, 0, 1);
    ctx.grid1.setLinSpaced(diffeo_size, 0, 1);

    vec<real> f0, f1;
    init_lorenz(f0, ctx, c, 0.1, 0.9);
    init_lorenz(f1, ctx);

    renorm_op<real> renorm(ctx, w0, w1);
    AutoDiffJacobian< renorm_op<real> > drenorm(renorm);
    mat<real> jac(f0.size(), f0.size());
    drenorm(f0, &f1, &jac);

    std::cerr << f0 << std::endl;
    std::cerr << f1 << std::endl;
    std::cerr << jac << std::endl;
    std::cerr << jac.eigenvalues() << std::endl;

    return EXIT_SUCCESS;
}
