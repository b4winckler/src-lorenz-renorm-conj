#include <iostream>
#include <csignal>

#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <unsupported/Eigen/AutoDiff>

#include "lorenz.h"


static char usage[] =
"usage: fixedpt w0 w1 c alpha dim prec niter\n"
"\n"
"Locates a (w0,w1)-renormalization fixed point; the  critical exponent\n"
"is 'alpha', the diffeomorphisms are truncated to 'dim' dimensions,\n"
"and 'prec' bits of precision are used for the calculations (and output).\n"
"The critical point 'c' is used as an initial guess.  The parameter 'niter'\n"
"determines how many times to iterate (more iterates lead to a better\n"
"approximation).  CTRL-C exists iteration prematurely and prints results.\n"
"\n"
"Output:\n"
"\n"
"    alpha  0  0 x0 ... xn\n"
"    c     v0 v1 y0 ... yn\n"
"\n"
"The coordinates (xk,yk) give piecewise linear approximations of the\n"
"diffeomorphisms.\n";


volatile std::sig_atomic_t signal_status = 0;


void signal_handler(int signal)
{
    signal_status = signal;
}

int main(int argc, char *argv[])
{
    if (argc != 8) {
        std::cerr << usage << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string w0(argv[1]);
    std::string w1(argv[2]);
    mpreal c(argv[3]);
    mpreal alpha(argv[4]);
    int precision = atoi(argv[6]);

    mpreal::set_default_prec(precision);
    mpreal eps = machine_epsilon(precision);
    mpreal desired_err2(pow(eps, 1.6));

    int dim = atoi(argv[5]);
    if (dim < 3)
        dim = 0;

    int niter = atoi(argv[7]);

    vec<mpreal> ctx;
    init_context(ctx, dim, alpha);

    vec<mpreal> f0, f1;
    init_lorenz(f0, ctx, c);
    init_lorenz(f1, ctx);

    vec<mpreal> pb0, pb1;
    thurston_guess(pb0, w0, w1);

    renorm_op<mpreal> renorm(ctx, w0, w1);

    std::signal(SIGINT, signal_handler);

    for (size_t i = 0; i < niter; ++i) {
        std::cerr << "[" << i << "]";
        thurston_op<mpreal> thurston(f0, ctx, w0, w1);
        mpreal err2 = 1;
        for (size_t j = 0; j < 1000 && err2 > desired_err2; ++j, pb0 = pb1) {
            thurston(pb0, &pb1);
            err2 = (pb0 - pb1).squaredNorm();
        }
        thurston.realization(f0, pb1);
        std::cerr << "\tthurston error = " << sqrt(err2);

        // Iteration should end with a Thurston step since it makes the error
        // smaller whereas the other two steps make the error worse.
        if (niter - 1 == i || SIGINT == signal_status) {
            std::cerr << std::endl;
            break;
        }

        if (dim > 2 && i % 2) {
            // Take step with renormalization operator on diffeos
            renorm(f0, &f1);
            mpreal err =
                (f0.tail(f0.size() - 3) - f1.tail(f0.size() - 3)).norm();
            f0.tail(f0.size() - 3) = f1.tail(f0.size() - 3);
            std::cerr << "\tdiffeo error = " << err;
        } else {
            // Take Newton step with modified renorm operator on c
            modrenorm_op<mpreal> modrenorm(f0, ctx, w0, w1);
            AutoDiffJacobian< modrenorm_op<mpreal> > dmodrenorm(modrenorm);
            vec<mpreal> p0(f0.head(3));
            vec<mpreal> p1(p0.size());
            mat<mpreal> bjac(p0.size(), p0.size());

            dmodrenorm(p0, &p1, &bjac);
            vec<mpreal> ph = bjac.fullPivLu().solve(p1);
            p0[0] -= ph[0];

            modrenorm.realization(f0, p0);
            std::cerr << "\tcritpt error = " << abs(ph[0]);
        }
        std::cerr << std::endl;
    }

    // Output result
    print_vec(ctx); std::cout << std::endl;
    print_vec(f0); std::cout << std::endl;

    // Log some more information to confirm fixed point
    std::cerr << "R^0(f) = " << f0.head(3).transpose() << std::endl;
    for (int i = 1; i <= 3; ++i, f0 = f1) {
        renorm(f0, &f1);
        std::cerr << "R^" << i << "(f) = " << f1.head(3).transpose() <<
            std::endl;
    }

    return EXIT_SUCCESS;
}

// TODO:
// - number of thurston iter should depend on precision
// - output in human readable form:
//      * f = F(c, v, phi)
//      * evaluation of f, Rf at the grid points (for plotting)
//      * largest eigenvalues of DR at f, descending
//      * distortion of phi?
//      * ...?
// - use library for faster eigenvalue estimates
