#include <iostream>
#include <mpreal.h>

using namespace mpfr;


void usage()
{
    std::cerr << "usage: lorenz\n\n";
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
    if (argc != 1)
        usage();

    mpreal::set_default_prec(512);

    return EXIT_SUCCESS;
}
