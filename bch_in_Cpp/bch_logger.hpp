#ifndef BCH_LOGGER_HPP
#define BCH_LOGGER_HPP

#include <iostream>

class BchLogger {
    public:
        BchLogger() = default;
        template<typename... Args>
        inline void log(Args&&... args) {
            if (enable_logging_) {
                ((std::cout << args), ...);
            }
        }
        bool enable_logging_;
};

#endif /* BCH_LOGGER_HPP */