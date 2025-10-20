#pragma once
#ifndef MATH_CORE_IO_HPP
#define MATH_CORE_IO_HPP

#include <iostream>

namespace math::io {

template<typename... Args>
void print(const Args&... args) {
    (std::cout << ... << args);
}

template<typename... Args>
void println(const Args&... args) {
    (std::cout << ... << args);
    std::cout << '\n';
}

}

#endif
