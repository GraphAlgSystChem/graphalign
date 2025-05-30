#ifndef FATAL_EXCEPTION_H
#define FATAL_EXCEPTION_H

#include <iostream>
#include <string>

// template from https://rollbar.com/blog/cpp-custom-exceptions/

class FatalException : public std::exception {
private:
    const char *message;

public:
    FatalException(const char *msg) : message(msg) {}
    const char * what() const noexcept override {
        return message;
    }
};

class TimeoutException : public std::exception {
private:
    const char *message;

public:
    TimeoutException(const char *msg) : message(msg) {}
    const char * what() const noexcept override {
        return message;
    }
};

#endif