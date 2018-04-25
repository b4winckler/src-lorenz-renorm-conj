#include <iostream>

#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/AutoDiff>

#include "lorenz.h"

static char usage[] =
"usage: deriv w0 w1 prec\n"
"\n"
"A Lorenz map f as output by 'fixedpt' is read on stdin (can be piped).\n"
"Outputs the first 3 eigenvalues of the derivative of the\n"
"(w0,w1)-renormalization at f, calculated with 'prec' bits of precision.\n";


int main(int argc, char *argv[])
{
    if (argc != 4) {
        std::cerr << usage << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string w0(argv[1]);
    std::string w1(argv[2]);
    int precision = atoi(argv[3]);

    mpreal::set_default_prec(precision);

    vec<mpreal> ctx;
    read_vec(ctx);

    vec<mpreal> f;
    read_vec(f);

    vec<mpreal> rf = f;

    renorm_op<mpreal> renorm(ctx, w0, w1);
    AutoDiffJacobian< renorm_op<mpreal> > drenorm(renorm);
    mat<mpreal> jac(f.size(), f.size());

    drenorm(f, &rf, &jac);
    vec<mpreal> evals = jac.topLeftCorner(3, 3).eigenvalues().array();
    std::sort(evals.data(), evals.data() + evals.size(),
            [](mpreal a, mpreal b) { return abs(a) > abs(b); } );

    std::cout << evals[0] << '\t' << evals[1] << '\t' << evals[2] << std::endl;

    return EXIT_SUCCESS;
}
