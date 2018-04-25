#include <iostream>
#include <vector>

#include "lorenz.h"

static char usage[] =
"usage: renorm w0 w1 nrenorm ngrid prec\n"
"\n"
"A Lorenz map as output by 'fixedpt' is read on stdin (can be piped).\n"
"Evaluates 'nrenorm' (w0,w1)-renormalizations on 'ngrid' points for each\n"
"branch, using 'prec' bits of precision.\n"
"Outputs pairs of columns with (x, y) data for each renormalization.\n";


int main(int argc, char *argv[])
{
    if (argc != 6) {
        std::cerr << usage << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string w0(argv[1]);
    std::string w1(argv[2]);
    int nrenorm = atoi(argv[3]);
    int ngrid = atoi(argv[4]);
    int precision = atoi(argv[5]);

    mpreal::set_default_prec(precision);

    vec<mpreal> ctx;
    read_vec(ctx);

    vec<mpreal> f;
    read_vec(f);

    vec<mpreal> rf = f;

    renorm_op<mpreal> renorm(ctx, w0, w1);

    vec<mpreal> x(2 * ngrid), y(2 * ngrid);
    std::vector< vec<mpreal> > xs, ys;

    for (int i = 0; i <= nrenorm; ++i) {
        // Evaluate f at evenly spaced points in interior of branches
        x.head(ngrid) = vec<mpreal>::LinSpaced(ngrid, 0, f[0]);
        x.tail(ngrid) = vec<mpreal>::LinSpaced(ngrid, f[0], 1);
        for (size_t j = 1; j < ngrid - 1; ++j) {
            apply(y[j], x[j], f, ctx);
            apply(y[ngrid + j], x[ngrid + j], f, ctx);
        }

        // Evaluate f at boundary of branches
        y[0] = f[1]; y[ngrid - 1] = 1;
        y[ngrid] = 0; y[2 * ngrid - 1] = f[2];

        xs.push_back(x); ys.push_back(y);

        if (i < nrenorm) {
            renorm(f, &rf);
            f = rf;
        }
    }

    // Output header
    std::cout << "x0\ty0";
    for (int i = 1; i <= nrenorm; ++i)
        std::cout << '\t' << "x" << i << '\t' << "y" << i;
    std::cout << std::endl;

    // Output columns of (x,y) data
    // NB. No need to print with high precision since this is for plotting
    for (size_t j = 0; j < 2 * ngrid; ++j) {
        std::cout << xs[0][j] << '\t' << ys[0][j];
        for (int i = 1; i <= nrenorm; ++i)
            std::cout << '\t' << xs[i][j] << '\t' << ys[i][j];
        std::cout << std::endl;
    }

    return EXIT_SUCCESS;
}
