#pragma once
#include <limits>

namespace rad{
  namespace constant{

    constexpr double InvalidEntry(){return std::numeric_limits<double>::quiet_NaN();};
    constexpr int InvalidIndex(){return -1;};


  }
}
