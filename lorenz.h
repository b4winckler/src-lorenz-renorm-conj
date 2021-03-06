#pragma once

#include <vector>

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


// A context holds the critical exponent alpha and the x-coordinates used for
// interpolating the diffeomorphisms.
template <typename scalar>
void init_context(vec<scalar> &ctx, size_t n, const scalar &alpha)
{
    ctx.resize(3 + 2 * n);
    ctx[0] = alpha;
    ctx.segment(3, n) = vec<scalar>::LinSpaced(n, 0, 1);
    ctx.tail(n) = vec<scalar>::LinSpaced(n, 0, 1);
}

// A Lorenz map holds the critical point c, boundary values (v0, v1), and the
// y-coordinates used for interpolating the diffeomorphisms.
template <typename scalar>
void init_lorenz(vec<scalar> &lorenz, const vec<scalar> &ctx, scalar c = 0.5,
        scalar v0 = 0, scalar v1 = 1)
{
    lorenz = ctx;
    lorenz[0] = c;
    lorenz[1] = v0;
    lorenz[2] = v1;
}

template <typename derived, typename scalar>
void apply(derived &y, const derived &x, const vec<derived> &lorenz,
        const vec<scalar> &ctx)
{
    // Adjust for critical point, then fold
    size_t k0 = 3, k1 = 3 + (ctx.size() - 3) / 2 - 1;
    bool left_branch = x < lorenz[0];
    bool right_branch = x > lorenz[0];
    if (left_branch) {
        y = x / lorenz[0];
        y = 1 - pow(1 - y, ctx[0]);
    } else if (right_branch) {
        y = (x - lorenz[0]) / (1 - lorenz[0]);
        y = pow(y, ctx[0]);
        k0 = k1 + 1;
        k1 = ctx.size() - 1;
    } else {
        error("can't evaluate lorenz map at x = " << x <<
                " (lorenz = " << lorenz.head(3).transpose() << ")");
    }

    // Interpolate diffeomorphism
    if (lorenz.size() > 5 && y > 0 && y < 1) {
        size_t k;
        // With autodiff library it is not possible (?) to take the floor of y
        // so we need to search to find the interval on which to interpolate.
        while (k1 > k0 + 1) {
            k = (k0 + k1) / 2;
            if (y < ctx[k]) {
                k1 = k;
            } else if (y > ctx[k]) {
                k0 = k;
            } else {
                k0 = k;
                k1 = k0 + 1;
            }
        }
        y = (ctx[k1] - y) * lorenz[k0] + (y - ctx[k0]) * lorenz[k1];
        y /= (ctx[k1] - ctx[k0]);
    }

    // Adjust for boundary value
    y = left_branch ? lorenz[1] + y * (1 - lorenz[1]) : y = y * lorenz[2];
}

template <typename derived, typename scalar>
void pull_back(derived &y, const derived &x, bool left_branch,
        const vec<derived> &lorenz, const vec<scalar> &ctx)
{
    // Adjust for boundary values
    size_t k0 = 3, k1 = 3 + (ctx.size() - 3) / 2 - 1;
    if (left_branch) {
        y = (x - lorenz[1]) / (1 - lorenz[1]);
    } else {
        y = x / lorenz[2];
        k0 = k1 + 1;
        k1 = ctx.size() - 1;
    }

    // HACK! Avoid rounding errors
    y = clamp(y, 0, 1);

    // Interpolate diffeomorphism
    if (lorenz.size() > 5 && y > 0 && y < 1) {
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
        y = (lorenz[k1] - y) * ctx[k0] + (y - lorenz[k0]) * ctx[k1];
        y /= (lorenz[k1] - lorenz[k0]);
    }

    // Fold, then adjust for critical point
    if (left_branch) {
        y = 1 - pow(1 - y, 1 / ctx[0]);
        y = y * lorenz[0];
    } else {
        y = pow(y, 1 / ctx[0]);
        y = lorenz[0] + (1 - lorenz[0]) * y;
    }
}

template <typename derived, typename scalar>
void iterate(derived &y, const derived &x, size_t n,
        const vec<derived> &lorenz, const vec<scalar> &ctx)
{
    y = x;
    for (size_t i = 0; i < n; ++i)
        apply(y, y, lorenz, ctx);
}

template <typename scalar>
void thurston_guess(vec<scalar> &shadow_orbit, const std::string &w0,
        const std::string &w1)
{
    size_t n = w0.size() + w1.size();
    scalar c(0.5);
    shadow_orbit.resize(2 * n);
    shadow_orbit.head(n) = vec<scalar>::LinSpaced(n, 0, c);
    shadow_orbit.tail(n) = vec<scalar>::LinSpaced(n, 1, c);
}

// The strange design of this code is due to interface with AutoDiff lib
template <typename scalar>
struct thurston_op {
    typedef vec<scalar> InputType;
    typedef vec<scalar> ValueType;

    vec<scalar> family2d;
    vec<scalar> ctx;
    std::string w0;
    std::string w1;
    std::string knead0;
    std::string knead1;

    thurston_op(const vec<scalar> &family2d_, const vec<scalar> &ctx_,
            const std::string &w0_, const std::string &w1_)
        : family2d(family2d_), ctx(ctx_), w0(w0_), w1(w1_)
    {
        knead0 = (w1 + w0).substr(1);
        knead1 = (w0 + w1).substr(1);
    }

    void realization(vec<scalar> &lorenz, const vec<scalar> &shadow_orbit)
    {
        lorenz = family2d;
        lorenz[1] = shadow_orbit[1];
        lorenz[2] = shadow_orbit[w0.size() + w1.size() + 1];
    }

    template <typename derived>
    void operator() (const vec<derived> &x, vec<derived> *py) const
    {
        size_t n0 = w0.size(), n1 = w1.size(), n = n0 + n1;

        // Choose member of 2d family to pull back with
        const derived &v0 = x[1], &v1 = x[n + 1];
        vec<derived> lorenz(family2d);
        lorenz[1] = v0;
        lorenz[2] = v1;

        vec<derived> &y = *py;
        y.resize(x.size());
        for (size_t i = 0, j = 0; j < n - 1; ++i, ++j)
            pull_back(y[i], x[i + 1], 'L' == knead0[j], lorenz, ctx);
        for (size_t i = n, j = 0; j < n - 1; ++i, ++j)
            pull_back(y[i], x[i + 1], 'L' == knead1[j], lorenz, ctx);

        // Choose pull back so that (v0, v1) are fixed under renormalization
        const derived &l = x[n1 - 1], &r = x[n + n0 - 1];
        y[n - 1] = l + v0 * (r - l);
        y[2 * n - 1] = l + v1 * (r - l);
    }
};

// The strange design of this code is due to interface with AutoDiff lib
template <typename scalar>
struct modrenorm_op {
    typedef vec<scalar> InputType;
    typedef vec<scalar> ValueType;

    vec<scalar> family3d;
    vec<scalar> ctx;
    std::string w0;
    std::string w1;

    modrenorm_op(const vec<scalar> &family3d_, const vec<scalar> &ctx_,
            const std::string &w0_, const std::string &w1_)
        : family3d(family3d_), ctx(ctx_), w0(w0_), w1(w1_) { }

    void realization(vec<scalar> &lorenz, const vec<scalar> &x)
    {
        lorenz = family3d;
        lorenz.head(3) = x;
    }

    template <typename derived>
    void operator() (const vec<derived> &x, vec<derived> *py) const
    {
        size_t n0 = w0.size(), n1 = w1.size();

        vec<derived> lorenz(family3d);
        lorenz.head(3) = x;

        // l = f^{n1 - 1)(0), vl = f^{n0}(l)
        derived l(0), vl;
        iterate(l, l, n1 - 1, lorenz, ctx);
        iterate(vl, l, n0, lorenz, ctx);

        // r = f^{n0 - 1)(1), vr = f^{n1}(r)
        derived r(1), vr;
        iterate(r, r, n0 - 1, lorenz, ctx);
        iterate(vr, r, n1, lorenz, ctx);

        vec<derived> &y = *py;
        y[0] = x[0] * (r - l - 1) + l;
        y[1] = x[1] * (r - l) - vl + l,
        y[2] = x[2] * (r - l) - vr + l;
    }
};

// The strange design of this code is due to interface with AutoDiff lib
template <typename scalar>
struct renorm_op {
    typedef vec<scalar> InputType;
    typedef vec<scalar> ValueType;

    vec<scalar> ctx;
    std::string w0;
    std::string w1;

    renorm_op(const vec<scalar> &ctx_, const std::string &w0_,
            const std::string &w1_)
        : ctx(ctx_), w0(w0_), w1(w1_) { }

    template <typename derived>
    void operator() (const vec<derived> &x, vec<derived> *py) const
    {
        size_t n0 = w0.size(), n1 = w1.size();

        // l = f^{n1 - 1)(0), vl = f^{n0}(l)
        derived l(0), vl;
        iterate(l, l, n1 - 1, x, ctx);
        iterate(vl, l, n0, x, ctx);

        // r = f^{n0 - 1)(1), vr = f^{n1}(r)
        derived r(1), vr;
        iterate(r, r, n0 - 1, x, ctx);
        iterate(vr, r, n1, x, ctx);

        // The return interval is [l, r] and the image of the first-return map
        // of f on this interval is [vl, vr]
        vec<derived> &y = *py;
        y[0] = (x[0] - l) / (r - l);
        y[1] = (vl - l) / (r - l);
        y[2] = (vr - l) / (r - l);

        if (x.size() > 5) {
            // Ensure endpoints of diffeos are fixed
            y[3] = 0;
            y[3 + (ctx.size() - 3) / 2 - 1] = 1;
            y[3 + (ctx.size() - 3) / 2] = 0;
            y[ctx.size() - 1] = 1;

            // Left diffeo
            size_t k = 4;
            for (; k < 3 + (ctx.size() - 3) / 2 - 1; ++k) {
                derived p(l + (1 - pow(1 - ctx[k], 1 / ctx[0])) * (x[0] - l));
                iterate(p, p, n0, x, ctx);
                y[k] = (p - vl) / (r - vl);
            }
            // Right diffeo
            for (k += 2; k < ctx.size() - 1; ++k) {
                derived p(x[0] + pow(ctx[k], 1 / ctx[0]) * (r - x[0]));
                iterate(p, p, n1, x, ctx);
                y[k] = (p - l) / (vr - l);
            }
        }
    }
};

template <typename scalar>
void print_vec(const vec<scalar> &x)
{
    if (x.size() == 0)
        return;

    std::cout << x[0].toString();
    for (size_t i = 1; i < x.size(); ++i)
        std::cout << '\t' << x[i].toString();
}

template <typename scalar>
void read_vec(vec<scalar> &v)
{
    std::string line;
    std::getline(std::cin, line);
    std::istringstream ss(line);

    std::vector<scalar> xs;
    scalar x;
    while (ss >> x)
        xs.push_back(x);

    v.resize(xs.size());
    for (size_t i = 0; i < xs.size(); ++i)
        v[i] = xs[i];
}
