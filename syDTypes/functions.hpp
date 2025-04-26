#pragma once

#include "complex.hpp"
#include "constants.hpp"

namespace sycl
{
    template <typename T>
    inline Complex<T> acospi(const Complex<T> &val)
    {
        return sycl::acos(val) / syConstants::PI;
    }

    template <typename T>
    inline Complex<T> asinpi(const Complex<T> &val)
    {
        return sycl::asin(val) / syConstants::PI;
    }

    template <typename T>
    inline  Complex<T> atan2(const Complex<T> &z)
    {
        if (z.real() == 0 && z.imag() == 0)
        {
            return Complex<T>(0, 0);
        }
        Complex<T> one(1, 0);
        Complex<T> z_squared = sycl::pow(z, 2);
        Complex<T> denominator = one + sycl::sqrt(one + z_squared);
        return T(2) * sycl::atan(z / denominator);
    }

    template <typename T>
    inline Complex<T> atanpi(const Complex<T> &val)
    {
        return sycl::atan(val) / syConstants::PI;
    }

    template <typename T>
    inline Complex<T> atan2pi(const Complex<T> &val)
    {
        return sycl::atan2(val) / syConstants::PI;
    }

    template <typename T>
    inline Complex<T> cbrt(const Complex<T>& val)
    {
        return sycl::pow(val, 1/3);
    }

    template <typename T>
    inline Complex<T> ceil(const Complex<T>& val)
    {
        Complex<T> a(sycl::ceil(val.real()), sycl::ceil(val.imag()));
    }

}