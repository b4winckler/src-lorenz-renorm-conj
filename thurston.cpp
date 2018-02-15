#include <iostream>
#include <mpreal.h>
#include <vector>

using namespace mpfr;


void error(const char *msg)
{
    std::cerr << msg << std::endl;
    exit(EXIT_FAILURE);
}

void pull_back_affine(
        // output
        mpreal &y,
        // input
        const mpreal &x, const mpreal &c, int left)
{
    y = left ? c - (1 - x) * c : c + x * (1 - c);
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

enum thurston_type {
    thurston_trivial,
    thurston_full,
    thurston_fixed_boundary_values,
    thurston_desired_boundary_values
};

enum { thurston_max_iterations = 1000000 };

static mpreal thurston_sqr_eps = 1e-40;

void thurston(
        // output
        mpreal &v0, mpreal &v1,
        // input
        const mpreal &c, const char *w0, const char *w1,
        thurston_type type, mpreal desired_v0 = 0, mpreal desired_v1 = 0)
{
    size_t n0 = strlen(w0), n1 = strlen(w1), n = n0 + n1;

    // Embed the words w1 + w0 and w0 + w1 in the unit interval by pulling back
    // the critical point with the full affine Lorenz map (nb. a full map has
    // pre-orbits defined for all words, and an affine map is fast to
    // calculate).  This provides an admissible initial guess to the Thurston
    // algorithm.
    std::vector<mpreal> x0(n + 1), x1(n + 1);
    x0[n] = x1[n] = c;
    for (size_t i = n, j = n - 1; i > 0; --i, --j) {
        pull_back_affine(x0[j], x0[i], c, (j >= n1 ? w0[j - n1] : w1[j]) == 'L');
        pull_back_affine(x1[j], x1[i], c, (j >= n0 ? w1[j - n0] : w0[j]) == 'L');
    }

    for (auto x : x0) std::cerr << x << ' '; std::cerr << std::endl;
    for (auto x : x1) std::cerr << x << ' '; std::cerr << std::endl;

    if (thurston_full == type) {
        desired_v0 = 0;
        desired_v1 = 1;
    }

    mpreal sqr_err(1), tmp;
    size_t count = 0;
    for (; sqr_err > thurston_sqr_eps && count < thurston_max_iterations;
            ++count) {
        v0 = x0[2];
        v1 = x1[2];

        switch (type) {
            case thurston_fixed_boundary_values:
                desired_v0 = v0;
                desired_v1 = v1;
                // fall through ...
            case thurston_full:
            case thurston_desired_boundary_values:
                x0[n] = x0[n0] + desired_v0 * (x1[n1] - x0[n0]);
                x1[n] = x0[n0] + desired_v1 * (x1[n1] - x0[n0]);
            case thurston_trivial: // do nothing
                break;
        }

        // This is the Thurston pull back step: pull back each point of the
        // embedded words in such a way that a fixed point of the pull back is
        // an actual orbit of the desired Lorenz family.  The magic of this
        // algorithm is that it always converges to a fixed point (as long as
        // the initial guess is admissible).
        sqr_err = 0;
        for (size_t i = 1; i <= n; ++i) {
            tmp = x0[i - 1];
            // nb. the '<' ensures x0[0] = c
            pull_back_lorenz(x0[i - 1], x0[i], c, v0, v1, x0[i - 1] < c);
            sqr_err += sqr(tmp - x0[i - 1]);

            tmp = x1[i - 1];
            // nb. the '<=' ensures x1[0] = c
            pull_back_lorenz(x1[i - 1], x1[i], c, v0, v1, x1[i - 1] <= c);
            sqr_err += sqr(tmp - x1[i - 1]);
        }

        for (auto x : x0) std::cerr << x << ' '; std::cerr << std::endl;
        for (auto x : x1) std::cerr << x << ' '; std::cerr << std::endl;
        std::cerr << sqr_err << std::endl;
    }

    std::cerr << "convergence after " << count << " iterations (out of max " <<
        thurston_max_iterations << ")" << std::endl;
    std::cerr << "error " << sqrt(sqr_err) << std::endl;
    // for (auto x : x0) std::cerr << x << ' '; std::cerr << std::endl;
    // for (auto x : x1) std::cerr << x << ' '; std::cerr << std::endl;
}

void usage()
{
    std::cerr << "usage: thurston w0 w1 c\n\n";
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
    if (argc != 4)
        usage();

    mpreal::set_default_prec(512);

    char *w0 = argv[1];
    char *w1 = argv[2];
    mpreal c(atof(argv[3]));
    mpreal v0, v1;

    thurston(v0, v1, c, w0, w1, thurston_trivial);

    return EXIT_SUCCESS;
}
