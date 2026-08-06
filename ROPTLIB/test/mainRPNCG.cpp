#include "test/TestRPNCGSPCA.h"
#include "randgen.h"

#include <ctime>
#include <iostream>

using namespace ROPTLIB;

int main(void)
{
    unsigned seed = 2;
    std::cout << "seed:" << seed << std::endl;
    genrandseed(seed);
    testRPNCGSPCA();
    return 0;
}
