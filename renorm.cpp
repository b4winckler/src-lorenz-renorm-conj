#include <iostream>
#include <vector>

#include "lorenz.h"


void read_vec(vec<mpreal> &v)
{
    std::string line;
    std::getline(std::cin, line);
    std::istringstream ss(line);

    std::vector<mpreal> xs;
    mpreal x;
    while (ss >> x)
        xs.push_back(x);

    v.resize(xs.size());
    for (size_t i = 0; i < xs.size(); ++i)
        v[i] = xs[i];
}

int main(int argc, char *argv[])
{
    if (argc != 6) {
        std::cerr << "usage: renorm w0 w1 nrenorm ngrid prec\n\n";
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

    vec<mpreal> rf;
    init_lorenz(rf, ctx);

    renorm_op<mpreal> renorm(ctx, w0, w1);

    vec<mpreal> x0(ngrid), y0(ngrid);
    vec<mpreal> x1(ngrid), y1(ngrid);

    for (int i = 0; i <= nrenorm; ++i, f = rf) {
        x0.setLinSpaced(ngrid, 0, f[0]);
        x1.setLinSpaced(ngrid, f[0], 1);

        for (size_t j = 1; j < ngrid - 1; ++j) {
            apply(y0[j], x0[j], f, ctx);
            apply(y1[j], x1[j], f, ctx);
        }

        y0[0] = f[1];
        y0[ngrid - 1] = 1;
        y1[0] = 0;
        y1[ngrid - 1] = f[2];

        print_vec(x0); std::cout << '\t'; print_vec(x1); std::cout << std::endl;
        print_vec(y0); std::cout << '\t'; print_vec(y1); std::cout << std::endl;

        renorm(f, &rf);
    }

    return EXIT_SUCCESS;
}
