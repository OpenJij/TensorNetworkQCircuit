#pragma once

#include <exception>

namespace qcircuit {
    class QCircuitException : public std::exception {
    private:
        std::string message;

    public:
        QCircuitException(const std::string& message) : message(message) {}

        const char* what() const noexcept override {
            return message.c_str();
        }
    };
}
