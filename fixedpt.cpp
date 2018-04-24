#include <iostream>

#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <unsupported/Eigen/AutoDiff>

#include "lorenz.h"


int main(int argc, char *argv[])
{
    if (argc != 8) {
        std::cerr << "usage: renorm w0 w1 c alpha ngrid prec niter\n\n";
        exit(EXIT_FAILURE);
    }

    std::string w0(argv[1]);
    std::string w1(argv[2]);
    mpreal c(argv[3]);
    mpreal alpha(argv[4]);
    int precision = atoi(argv[6]);

    mpreal::set_default_prec(precision);
    mpreal eps = machine_epsilon(precision);
    std::cerr << "eps = " << eps << std::endl;
    mpreal desired_err2(pow(eps, 1.6));

    int diffeo_size = atoi(argv[5]);
    if (diffeo_size < 3)
        diffeo_size = 0;

    int niter = atoi(argv[7]);

    context ctx;
    ctx.alpha = alpha;
    if (diffeo_size > 2) {
        ctx.grid0.setLinSpaced(diffeo_size, 0, 1);
        ctx.grid1.setLinSpaced(diffeo_size, 0, 1);
    }

    vec<mpreal> f0, f1;
    init_lorenz(f0, ctx, c);
    init_lorenz(f1, ctx);

    vec<mpreal> pb0, pb1;
    thurston_guess(pb0, w0, w1);

    renorm_op<mpreal> renorm(ctx, w0, w1);
    AutoDiffJacobian< renorm_op<mpreal> > drenorm(renorm);
    mat<mpreal> jac(f0.size(), f0.size());

    for (size_t i = 0; i < niter; ++i) {
        std::cerr << "i = " << i << std::endl;
        thurston_op<mpreal> thurston(f0, ctx, w0, w1);
        mpreal err2 = 1;
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

        if (diffeo_size > 2 && i % 2) {
            renorm(f0, &f1);
            f0.tail(f0.size() - 3) = f1.tail(f0.size() - 3);
            std::cerr << "renorm   = " << f1.transpose() << std::endl;
        } else {
            boundary_op<mpreal> boundary(f0, ctx, w0, w1);
            AutoDiffJacobian< boundary_op<mpreal> > dboundary(boundary);
            vec<mpreal> p0(f0.head(3));
            vec<mpreal> p1(p0.size());
            mat<mpreal> bjac(p0.size(), p0.size());

            for (int k = 0; k < 1; ++k) {
                dboundary(p0, &p1, &bjac);
                vec<mpreal> ph = bjac.fullPivLu().solve(p1);
                p0[0] -= ph[0];
                // std::cerr << "\t[" << ph.squaredNorm() << "] " <<
                //     p0.transpose() << std::endl;
            }

            boundary.realization(f0, p0);
            std::cerr << "boundary = " << p0.transpose() << std::endl;
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

    vec<mpreal> evals = jac.eigenvalues().array();
    std::sort(evals.data(), evals.data() + evals.size(),
            [](mpreal a, mpreal b) { return abs(a) > abs(b); } );
    std::cerr << "eigvals = " << evals.transpose().head(3) << std::endl;

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
