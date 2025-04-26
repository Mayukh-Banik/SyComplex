#pragma once

#include "complex.hpp"
#include "constants.hpp"

namespace sycl
{



    template <typename T, typename U = T>
    inline constexpr Complex<T> pow(const Complex<T> &val, const Complex<U> &power)
    {
        return Complex<T>();
    }

    template <typename T, typename U = T>
    inline constexpr Complex<T> pow(const T val, const Complex<U> &power)
    {
        return Complex<T>();
    }

    template <typename T>
    inline constexpr Complex<T> sin(const Complex<T> &val)
    {
        return Complex<T>(
            sycl::sin(val.real()) * sycl::cosh(val.imag()),
            sycl::cos(val.real()) * sycl::sinh(val.imag()));
    }


    template <typename T>
    inline constexpr Complex<T> cos(const Complex<T> &val)
    {
        return Complex<T>(
            sycl::cos(val.real()) * sycl::cosh(val.imag()),
            -1 * (sycl::sin(val.real()) * sycl::sinh(val.imag())));
    }


    template <typename T>
    inline constexpr Complex<T> tan(const Complex<T> &val)
    {
        Complex<T> sin_val = sycl::sin(val);
        Complex<T> cos_val = sycl::cos(val);
        return sin_val / cos_val;
    }


    template <typename T>
    inline constexpr Complex<T> sinh(const Complex<T> &val)
    {
        return Complex<T>(
            sycl::sinh(val.real()) * sycl::cos(val.imag()),
            sycl::cosh(val.real()) * sycl::sin(val.imag()));
    }


    template <typename T>
    inline constexpr Complex<T> cosh(const Complex<T> &val)
    {
        return Complex<T>(
            sycl::cosh(val.real()) * sycl::cos(val.imag()),
            sycl::sinh(val.real()) * sycl::sin(val.imag()));
    }


    template <typename T>
    inline constexpr Complex<T> tanh(const Complex<T> &val)
    {
        Complex<T> sinh_val = sinh(val);
        Complex<T> cosh_val = cosh(val);
        return sinh_val / cosh_val;
    }


}