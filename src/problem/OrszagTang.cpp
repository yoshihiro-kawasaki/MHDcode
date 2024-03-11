#include "../data.hpp"
#include "../driver.hpp"

void GridInfo::SetGridInfo()
{
    nx1_   = 64;
    x1min_ = -0.5;
    x1max_ = 0.5;
    nx2_   = 1;
    x2min_ = 0.0;
    x2max_ = 0.0;
    nx3_   = 1;
    x3min_ = 0.0;
    x3max_ = 0.0;
    ngh_   = 2;
    is_set_gridinfo_ = true;
}