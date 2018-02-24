#include <unsupported/Eigen/MPRealSupport>
#include <unsupported/Eigen/AutoDiff>

#include <iostream>

using namespace mpfr;
using namespace Eigen;

template<typename R>
struct test_operator
{
    typedef Matrix<R, 2, 1> InputType;
    typedef Matrix<R, 2, 1> ValueType;

    template <typename T1, typename T2>
    void operator() (const T1 &x, T2 *yp) const
    {
        T2 &y = *yp;
#if 0
        // linear example
        Matrix<R, 2, 2> m; m <<
            1, 2,
            1, 0;
        y = m * x;
#elif 0
        // same as above
        y(0) = x(0) + 2 * x(1);
        y(1) = x(0);
#else
        // nonlinear example
        y(0) = exp(-0.5 * x(0) + pow(x(1), 3));
        y(1) = sin(x(0)) - cos(x(1));
#endif
    }
};

void usage()
{
    std::cerr << "usage: test\n\n";
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
    if (argc != 1)
        usage();

    mpreal::set_default_prec(512);

    test_operator<mpreal> op;
    AutoDiffJacobian< test_operator<mpreal> > ad_op(op);

    Matrix<mpreal, 2, 2> dy;
    Matrix<mpreal, 2, 1> x, y;
    x << 1, -1;
    ad_op(x, &y, &dy);

    std::cout << "x =\n" << x << std::endl;
    std::cout << "y =\n" << y << std::endl;
    std::cout << "dy =\n" << dy << std::endl;

    return EXIT_SUCCESS;
}
