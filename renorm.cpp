#define USE_MPFR 1

#include <iostream>
#include <limits>

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
    std::cerr << __FILE__ << ':' << __LINE__ << ": " << msg << std::endl; \
    exit(EXIT_FAILURE); \
}

#define clamp(x, x0, x1) ((x) < (x0) ? (x0) : ((x) > (x1) ? (x1) : (x)))


using namespace Eigen;


template <typename t> using vec = Matrix<t, Dynamic, 1>;
template <typename t> using mat = Matrix<t, Dynamic, Dynamic>;


struct context {
    real alpha;
    vec<real> grid0;
    vec<real> grid1;
};


template <typename scalar>
void init_lorenz(vec<scalar> &lorenz, const context &ctx, scalar c = 0.5,
        scalar v0 = 0, scalar v1 = 1)
{
    size_t n = 3 + ctx.grid0.size() + ctx.grid1.size();
    lorenz.resize(n);
    lorenz[0] = c;
    lorenz[1] = v0;
    lorenz[2] = v1;

    if (n > 3) {
        lorenz.segment(3, ctx.grid0.size()) = ctx.grid0;
        lorenz.tail(ctx.grid1.size()) = ctx.grid1;
    }
}

template <typename scalar>
void apply(scalar &y, const scalar &x, const vec<scalar> &lorenz,
        const context &ctx)
{
    // Adjust for critical point, then fold
    size_t i0, k0, k1;
    bool left_branch = x < lorenz[0];
    bool right_branch = x > lorenz[0];
    if (left_branch) {
        y = x / lorenz[0];
        y = 1 - pow(1 - y, ctx.alpha);
        i0 = k0 = 3;
        k1 = k0 + ctx.grid0.size() - 1;
    } else if (right_branch) {
        y = (x - lorenz[0]) / (1 - lorenz[0]);
        y = pow(y, ctx.alpha);
        i0 = k0 = 3 + ctx.grid0.size();
        k1 = k0 + ctx.grid1.size() - 1;
    } else {
        error("can't evaluate lorenz map at x = " << x <<
                " (lorenz = " << lorenz.head(3).transpose() << ")");
    }

    // Interpolate diffeomorphism
    if (lorenz.size() > 5 && y > 0 && y < 1) {
        const vec<real> &grid = left_branch ? ctx.grid0 : ctx.grid1;
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
    y = left_branch ? lorenz[1] + y * (1 - lorenz[1]) : y = y * lorenz[2];
}

template <typename scalar>
void pull_back(scalar &y, const scalar &x, bool left_branch,
        const vec<scalar> &lorenz, const context &ctx)
{
    // Adjust for boundary values
    size_t l, k0, k1;
    if (left_branch) {
        y = (x - lorenz[1]) / (1 - lorenz[1]);
        l = k0 = 3;
        k1 = k0 + ctx.grid0.size() - 1;
    } else {
        y = x / lorenz[2];
        l = k0 = 3 + ctx.grid0.size();
        k1 = k0 + ctx.grid1.size() - 1;
    }

    // HACK! Avoid rounding errors
    y = clamp(y, 0, 1);

    // Interpolate diffeomorphism
    if (lorenz.size() > 5 && y > 0 && y < 1) {
        const vec<real> &grid = left_branch ? ctx.grid0 : ctx.grid1;
        size_t k;
        // Search to find the interval on which to interpolate.
        while (k1 > k0 + 1) {
            k = (k0 + k1) / 2;
            if (y < lorenz[k]) {
                k1 = k;
            } else if (y > lorenz[k]) {
                k0 = k;
            } else {
                k0 = k;
                k1 = k0 + 1;
            }
        }
        y = (lorenz[k1] - y) * grid[k0 - l] + (y - lorenz[k0]) * grid[k1 - l];
        y /= (lorenz[k1] - lorenz[k0]);
    }

    // Fold, then adjust for critical point
    if (left_branch) {
        y = 1 - pow(1 - y, 1 / ctx.alpha);
        y = y * lorenz[0];
    } else {
        y = pow(y, 1 / ctx.alpha);
        y = lorenz[0] + (1 - lorenz[0]) * y;
    }
}

template <typename scalar>
void iterate(scalar &y, const scalar &x, size_t n, const vec<scalar> &lorenz, const context &ctx)
{
    y = x;
    for (size_t i = 0; i < n; ++i)
        apply(y, y, lorenz, ctx);
}

template <typename scalar>
void thurston_guess(vec<scalar> &shadow_cycle, const std::string &w0,
        const std::string &w1)
{
    size_t n = w0.size() + w1.size();
    scalar c(0.5);
    shadow_cycle.resize(2 * n);
#if 1   // Put critical values and their images in order
    shadow_cycle.head(n) = vec<scalar>::LinSpaced(n, 0, c);
    shadow_cycle.tail(n) = vec<scalar>::LinSpaced(n, 1, c);
#else   // Put whole orbits in kneading sequence order
    shadow_cycle[n - 1] = shadow_cycle[2 * n - 1] = c;

    for (size_t i = n - 1; i > 0; --i) {
        shadow_cycle[i - 1] = 'L' == knead0[i - 1]
            ?  c - (1 - shadow_cycle[i]) * c
            : c + shadow_cycle[i] * (1 - c);
        shadow_cycle[n + i - 1] = 'L' == knead1[i - 1]
            ?  c - (1 - shadow_cycle[n + i]) * c
            : c + shadow_cycle[n + i] * (1 - c);
    }
#endif
}

template <typename scalar>
struct thurston_op {
    typedef vec<scalar> InputType;
    typedef vec<scalar> ValueType;

    vec<scalar> family2d;
    context ctx;
    std::string w0;
    std::string w1;
    std::string knead0;
    std::string knead1;

    thurston_op(const vec<scalar> &family2d_, const context &ctx_,
            const std::string &w0_, const std::string &w1_)
        : family2d(family2d_), ctx(ctx_), w0(w0_), w1(w1_)
    {
        knead0 = (w1 + w0).substr(1);
        knead1 = (w0 + w1).substr(1);
    }

    void realization(vec<scalar> &lorenz, const vec<scalar> &shadow_cycle)
    {
        lorenz = family2d;
        lorenz[1] = shadow_cycle[1];
        lorenz[2] = shadow_cycle[w0.size() + w1.size() + 1];
    }

    template <typename t>
    void operator() (const vec<t> &x, vec<t> *py) const
    {
        size_t n0 = w0.size(), n1 = w1.size(), n = n0 + n1;

        // Choose member of 2d family to pull back with
#if 1   // Assume shadow cycle is admissible
        const t &v0 = x[1], &v1 = x[n + 1];
#else   // Sanity check: pick (v0, v1) so that all pull backs are defined
        t v0 = x[1], v1 = x[n + 1];
        for (size_t i = 0, j = 0; j < n - 1; ++i, ++j) {
            if ('L' == knead0[j] && x[i + 1] < v0)
                v0 = x[i + 1];
            if ('R' == knead0[j] && x[i + 1] > v1)
                v1 = x[i + 1];
        }
        for (size_t i = n, j = 0; j < n - 1; ++i, ++j) {
            if ('L' == knead1[j] && x[i + 1] < v0)
                v0 = x[i + 1];
            if ('R' == knead1[j] && x[i + 1] > v1)
                v1 = x[i + 1];
        }
#endif
        vec<t> lorenz(family2d);
        lorenz[1] = v0;
        lorenz[2] = v1;

        vec<t> &y = *py;
        y.resize(x.size());
        for (size_t i = 0, j = 0; j < n - 1; ++i, ++j)
            pull_back(y[i], x[i + 1], 'L' == knead0[j], lorenz, ctx);
        for (size_t i = n, j = 0; j < n - 1; ++i, ++j)
            pull_back(y[i], x[i + 1], 'L' == knead1[j], lorenz, ctx);

        // Choose pull back so that (v0, v1) are fixed under renormalization
        const t &l = x[n1 - 1], &r = x[n + n0 - 1];
        y[n - 1] = l + v0 * (r - l);
        y[2 * n - 1] = l + v1 * (r - l);
    }
};

template <typename scalar>
struct boundary_op {
    typedef vec<scalar> InputType;
    typedef vec<scalar> ValueType;

    vec<scalar> family3d;
    context ctx;
    std::string w0;
    std::string w1;

    boundary_op(const vec<scalar> &family3d_, const context &ctx_,
            const std::string &w0_, const std::string &w1_)
        : family3d(family3d_), ctx(ctx_), w0(w0_), w1(w1_) { }

    void realization(vec<scalar> &lorenz, const vec<scalar> &x)
    {
        lorenz = family3d;
        lorenz.head(3) = x;
    }

    template <typename t>
    void operator() (const vec<t> &x, vec<t> *py) const
    {
        size_t n0 = w0.size(), n1 = w1.size();

        vec<t> lorenz(family3d);
        lorenz.head(3) = x;

        // l = f^{n1 - 1)(0), vl = f^{n0}(l)
        t l(0), vl;
        iterate(l, l, n1 - 1, lorenz, ctx);
        iterate(vl, l, n0, lorenz, ctx);

        // r = f^{n0 - 1)(1), vr = f^{n1}(r)
        t r(1), vr;
        iterate(r, r, n0 - 1, lorenz, ctx);
        iterate(vr, r, n1, lorenz, ctx);

        vec<t> &y = *py;
        y[0] = x[0] * (r - l - 1) + l;
        y[1] = x[1] * (r - l) - vl + l,
        y[2] = x[2] * (r - l) - vr + l;
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
        y[0] = (x[0] - l) / (r - l);
        y[1] = (vl - l) / (r - l);
        y[2] = (vr - l) / (r - l);

        if (x.size() > 5) {
            // Ensure endpoints of diffeos are fixed
            y[3] = 0;
            y[3 + ctx.grid0.size()] = 0;
            y[3 + ctx.grid0.size() - 1] = 1;
            y[3 + ctx.grid0.size() + ctx.grid1.size() - 1] = 1;

            // Left diffeo
            for (size_t k = 1; k < ctx.grid0.size() - 1; ++k) {
                t p(l + (1 - pow(1 - ctx.grid0[k], 1 / ctx.alpha)) *
                        (x[0] - l));
                iterate(p, p, n0, x, ctx);
                y[3 + k] = (p - vl) / (r - vl);
            }

            // Right diffeo
            for (size_t k = 1; k < ctx.grid1.size() - 1; ++k) {
                t p(x[0] + pow(ctx.grid0[k], 1 / ctx.alpha) * (r - x[0]));
                iterate(p, p, n1, x, ctx);
                y[3 + ctx.grid0.size() + k] = (p - l) / (vr - l);
            }
        }
    }
};

int main(int argc, char *argv[])
{
    if (argc != 7) {
        std::cerr << "usage: renorm w0 w1 c alpha n prec\n\n";
        exit(EXIT_FAILURE);
    }

    std::string w0(argv[1]);
    std::string w1(argv[2]);
#if USE_MPFR
    real c(argv[3]);
    real alpha(argv[4]);
    int precision = atoi(argv[6]);

    mpreal::set_default_prec(precision);
    real eps = machine_epsilon(precision);
    std::cerr << "eps = " << eps << std::endl;
#else
    real c = atof(argv[3]);
    real alpha = atof(argv[4]);
    real eps = std::numeric_limits<real>::epsilon();
#endif
    real desired_err2(pow(eps, 1.6));

    int diffeo_size = atoi(argv[5]);
    if (diffeo_size < 3)
        diffeo_size = 0;

    context ctx;
    ctx.alpha = alpha;
    if (diffeo_size > 2) {
        ctx.grid0.setLinSpaced(diffeo_size, 0, 1);
        ctx.grid1.setLinSpaced(diffeo_size, 0, 1);
    }

    vec<real> f0, f1;
    init_lorenz(f0, ctx, c);
    init_lorenz(f1, ctx);

    vec<real> pb0, pb1;
    thurston_guess(pb0, w0, w1);

    renorm_op<real> renorm(ctx, w0, w1);
    AutoDiffJacobian< renorm_op<real> > drenorm(renorm);
    mat<real> jac(f0.size(), f0.size());

    for (size_t i = 0; i < 100; ++i) {
        std::cerr << "i = " << i << std::endl;
        thurston_op<real> thurston(f0, ctx, w0, w1);
        real err2 = 1;
        for (size_t j = 0; j < 1000 && err2 > desired_err2; ++j, pb0 = pb1) {
            thurston(pb0, &pb1);
            err2 = (pb0 - pb1).squaredNorm();
            // std::cerr << "\t[" << j << "] " << pb1.transpose() <<
                // std::endl;
        }
        thurston.realization(f0, pb1);
        std::cerr << "thurston = " << f0.transpose() << std::endl;
        std::cerr << "    err2 = " << err2 << std::endl;
        if (99 == i)
            break;

        if (i % 2) {
            boundary_op<real> boundary(f0, ctx, w0, w1);
            AutoDiffJacobian< boundary_op<real> > dboundary(boundary);
            vec<real> p0(f0.head(3));
            vec<real> p1(p0.size());
            mat<real> bjac(p0.size(), p0.size());

            for (int k = 0; k < 1; ++k) {
                dboundary(p0, &p1, &bjac);
                vec<real> ph = bjac.fullPivLu().solve(p1);
                p0[0] -= ph[0];
                // std::cerr << "\t[" << ph.squaredNorm() << "] " <<
                //     p0.transpose() << std::endl;
            }

            boundary.realization(f0, p0);
            std::cerr << "boundary = " << p0.transpose() << std::endl;
        } else {
            renorm(f0, &f1);
            f0.tail(f0.size() - 3) = f1.tail(f0.size() - 3);
            std::cerr << "renorm   = " << f1.transpose() << std::endl;
        }
    }

    std::cerr << "\n========================\n" << std::endl;

    drenorm(f0, &f1, &jac);
    std::cerr << "R^0(f) = " << f0.transpose() << std::endl;
    for (int i = 1; i <= 1; ++i, f0 = f1) {
        renorm(f0, &f1);
        std::cerr << "R^" << i << "(f) = " << f1.transpose() <<
            std::endl;
    }

    vec<real> evals = jac.eigenvalues().array();
    std::sort(evals.data(), evals.data() + evals.size(),
            [](real a, real b) { return abs(a) > abs(b); } );
    std::cerr << "eigvals = " << evals.transpose().head(3) << std::endl;

    return EXIT_SUCCESS;
}
