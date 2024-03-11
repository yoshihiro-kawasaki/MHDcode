#ifndef DEFS_HPP
#define DEFS_HPP

// #######################################################################
// #######################################################################
// #######################################################################

#define SQR(x) ( (x)*(x) )
#define CUB(x) ( (x)*(x)*(x) )
#define SIGN(x) ( (x >= 0.0 ? 1.0 : -1.0))

#define FAILED_TO_OPEN(file_stream, file_name) if(!file_stream) {std::cerr << "##ERROR : Failed to open " << file_name << std::endl; exit(1); }

// #######################################################################
// #######################################################################
// #######################################################################

constexpr double GAMMA = 5.0 / 3.0;
constexpr double CR    = 0.18;
constexpr double CFL   = 0.3;

// #######################################################################
// #######################################################################
// #######################################################################

enum PrimitiveVariablesIndex
{
    IRHO,
    IV1,
    IV2,
    IV3,
    IPR,
    IB1,
    IB2,
    IB3,
    IPSI,
    NVAR
};


enum ConservativeVariablesIndex
{
    IM1 = 1,
    IM2,
    IM3,
    IEN
};

// #######################################################################
// #######################################################################
// #######################################################################


// #######################################################################
// #######################################################################
// #######################################################################

#endif /* DEFS_HPP */