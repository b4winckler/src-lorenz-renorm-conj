#pragma once

#include <unsupported/Eigen/MPRealSupport>

using namespace Eigen;
using namespace mpfr;


#define error(msg) \
{ \
    std::cerr << __FILE__ << ':' << __LINE__ << ": " << msg << std::endl; \
    exit(EXIT_FAILURE); \
}

#define clamp(x, x0, x1) ((x) < (x0) ? (x0) : ((x) > (x1) ? (x1) : (x)))


template <typename t> using vec = Matrix<t, Dynamic, 1>;
template <typename t> using mat = Matrix<t, Dynamic, Dynamic>;


struct context {
    mpreal alpha;
    vec<mpreal> grid0;
    vec<mpreal> grid1;
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
        const vec<mpreal> &grid = left_branch ? ctx.grid0 : ctx.grid1;
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
        const vec<mpreal> &grid = left_branch ? ctx.grid0 : ctx.grid1;
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
