#pragma once
#include <limits>

namespace rad{
  namespace constant{
    /**
     *    InvalidEntry
     */
    // Generic template declaration
    template<typename T>
    constexpr T InvalidEntry();

    template<typename T>
    constexpr T InvalidEntry() {
      return std::numeric_limits<T>::quiet_NaN();
    }
  
    // Specialization for int
    template<>
    constexpr int InvalidEntry<int>() {
      return std::numeric_limits<int>::max();
    }

    // Specialization for unsigned int
    template<>
    constexpr unsigned int InvalidEntry<unsigned int>() {
      return std::numeric_limits<unsigned int>::max();
    }

    template<typename T>
    inline bool IsInvalidEntry(const T& entry){
      return entry==InvalidEntry<T>();
    }
    // constexpr double InvalidEntry(){return std::numeric_limits<double>::quiet_NaN();};
    constexpr int InvalidIndex(){return -1;};


  }
}
