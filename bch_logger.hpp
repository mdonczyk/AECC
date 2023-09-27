#ifndef BCH_LOGGER_HPP
#define BCH_LOGGER_HPP

#include <iostream>

namespace bch_logger {
    inline bool enable_logging{}; // off by default

    template<typename... Args>
    inline void log(Args&&... args) {
        if (enable_logging) {
            ((std::cout << args), ...);
        }
    }
}

#endif /* BCH_LOGGER_HPP */