#pragma once

#include <type_traits>
#include <complex>
#include <sycl/sycl.hpp>
#include "constants.hpp"

#define CONSTEXPR_IF_UT(x, y)           \
    if constexpr (std::is_same_v<T, U>) \
    {                                   \
        return x;                       \
    }                                   \
    else                                \
    {                                   \
        return y;                       \
    }

namespace sycl
{

    template <typename T>
    struct Complex
    {
        T _real;
        T _imag;
#define a this->_real
#define b this->_imag
#define c rhs._real
#define d rhs._imag
        constexpr Complex(T real = 0, T imag = 0) : _real(real), _imag(imag) {}

        operator std::complex<T> const() { return std::complex<T>(a, b); }

        constexpr Complex(const std::complex<T> &complex)
        {
            a = complex.real;
            b = complex.imag;
        }

        inline constexpr T real() const { return a; }
        inline constexpr T imag() const { return b; }
        inline constexpr void real(T val) { a = val; }
        inline constexpr void imag(T val) { b = val; }

        Complex<T> &operator++()
        {
            a++;
            b++;
            return *this;
        }

        Complex<T> operator++(int)
        {
            Complex<T> temp = *this;
            a++;
            b++;
            return temp;
        }

        Complex<T> &operator--()
        {
            a--;
            b--;
            return *this;
        }

        Complex<T> operator--(int)
        {
            Complex<T> temp = *this;
            a--;
            b--;
            return temp;
        }

        template <typename U = T>
        constexpr Complex<T> operator+(const Complex<U> &rhs)
        {
            CONSTEXPR_IF_UT(
                Complex<T>(a + c, b + d),
                Complex<T>(a + (T)c, b + (T)d))
        }

        template <typename U = T>
        constexpr Complex<T> operator+(const U &rhs)
        {
            CONSTEXPR_IF_UT(
                Complex<T>(a + rhs, b),
                Complex<T>(a + (T)rhs, b))
        }

        template <typename U = T>
        constexpr Complex<T> operator-(const Complex<U> &rhs)
        {
            CONSTEXPR_IF_UT(
                Complex<T>(a - c, b - d),
                Complex<T>(a - (T)c, b - (T)d))
        }

        template <typename U = T>
        constexpr Complex<T> operator-(const U &rhs)
        {
            CONSTEXPR_IF_UT(
                Complex<T>(a - rhs, b),
                Complex<T>(a - (T)rhs, b))
        }

        template <typename U = T>
        constexpr Complex<T> operator*(const Complex<U> &rhs)
        {
            CONSTEXPR_IF_UT(
                Complex<T>(
                    a * c - b * d,
                    a * d - b * c),
                Complex<T>(
                    a * (T)c - b * (T)d,
                    a * (T)d - b * (T)c))
        }

        template <typename U = T>
        constexpr Complex<T> operator*(const U &rhs)
        {
            CONSTEXPR_IF_UT(
                Complex<T>(a * rhs, b * rhs),
                Complex<T>(a * (T)rhs, b * (T)rhs))
        }

        template <typename U = T>
        constexpr Complex<T> operator/(const Complex<U> &rhs)
        {
            T num1, num2, denom;
            if constexpr (std::is_same_v<T, U>)
            {
                num1 = a * c + b * d;
                num2 = b * c - a * d;
                denom = c * c + d * d;
            }
            else
            {
                num1 = a * (T)c + b * (T)d;
                num2 = b * (T)c - a * (T)d;
                denom = (T)(c * c + d * d);
            }
            if (denom == 0)
            {
                return Complex<T>(0, 0);
            }
            else
            {
                return Complex<T>(num1 / denom, num2 / denom);
            }
        }

        template <typename U = T>
        constexpr Complex<T> operator/(const U &rhs)
        {
            CONSTEXPR_IF_UT(
                Complex<T>(a / rhs, b / rhs),
                Complex<T>(a / (T)rhs, b / (T)rhs))
        }

        constexpr Complex<T> operator+() const
        {
            return Complex<T>(+a, +b);
        }

        constexpr Complex<T> operator-() const
        {
            return Complex<T>(-a, -b);
        }

        constexpr bool operator!() const
        {
            return !a && !b;
        }

        constexpr Complex<T> operator~() const
        {
            return Complex<T>(~a, ~b);
        }

        template <typename U = T>
        inline constexpr bool operator==(const Complex<U> &rhs) const
        {
            CONSTEXPR_IF_UT(
                a == c && b == d,
                a == (T)c && b == (T)d)
        }

        template <typename U = T>
        inline constexpr bool operator==(U &rhs) const
        {
            CONSTEXPR_IF_UT(
                false,
                false)
        }

        template <typename U = T>
        inline constexpr bool operator!=(const Complex<U> &rhs) const
        {
            return !(*this == rhs);
        }

        template <typename U = T>
        inline constexpr bool operator!=(U &rhs) const
        {
            return !(*this == rhs);
        }

        template <typename U = T>
        inline constexpr Complex<T> &operator+=(const Complex<U> &rhs) const
        {
            if constexpr (std::is_same_v<U, T>)
            {
                a = a + c;
                b = b + d;
            }
            else
            {
                a = a + (T)c;
                b = b + (T)d;
            }
            return *this;
        }

        template <typename U = T>
        inline constexpr Complex<T> &operator+=(const U &rhs) const
        {
            if constexpr (std::is_same_v<U, T>)
            {
                a = a + rhs;
            }
            else
            {
                a = a + (T)rhs;
            }
            return *this;
        }

        template <typename U = T>
        inline constexpr Complex<T> &operator-=(const Complex<U> &rhs) const
        {
            if constexpr (std::is_same_v<U, T>)
            {
                a = a - c;
                b = b - d;
            }
            else
            {
                a = a - (T)c;
                b = b - (T)d;
            }
            return *this;
        }

        template <typename U = T>
        inline constexpr Complex<T> &operator-=(const U &rhs) const
        {
            if constexpr (std::is_same_v<U, T>)
            {
                a = a - rhs;
            }
            else
            {
                a = a - (T)rhs;
            }
            return *this;
        }

        template <typename U = T>
        inline constexpr Complex<T> &operator*=(const Complex<U> &rhs) const
        {
            T t1, t2;
            if constexpr (std::is_same_v<U, T>)
            {
                t1 = a * c - b * d;
                t2 = a * d - b * c;
            }
            else
            {
                t1 = a * (T)c - b * (T)d;
                t2 = a * (T)d - b * (T)c;
            }
            a = t1;
            b = t2;
            return *this;
        }

        template <typename U = T>
        inline constexpr Complex<T> &operator*=(const U &rhs) const
        {
            if constexpr (std::is_same_v<U, T>)
            {
                a = a * rhs;
            }
            else
            {
                a = a * (T)rhs;
            }
            return *this;
        }

        template <typename U = T>
        inline constexpr Complex<T> &operator/=(const Complex<U> &rhs) const
        {
            T num1, num2, denom;
            if constexpr (std::is_same_v<T, U>)
            {
                num1 = a * c + b * d;
                num2 = b * c - a * d;
                denom = c * c + d * d;
            }
            else
            {
                num1 = a * (T)c + b * (T)d;
                num2 = b * (T)c - a * (T)d;
                denom = (T)(c * c + d * d);
            }
            if (denom == 0)
            {
                a = 0;
                b = 0;
                return *this;
            }
            else
            {
                a = num1 / denom;
                b = num2 / denom;
                return *this;
            }
        }

        template <typename U = T>
        inline constexpr Complex<T> &operator/=(const U &rhs) const
        {
            if constexpr (std::is_same_v<U, T>)
            {
                a = a / rhs;
                b = b / rhs;
            }
            else
            {
                a = a / (T)rhs;
                b = b / (T)rhs;
            }
            return *this;
        }

        friend std::ostream &operator<<(std::ostream &os, const Complex<T> &rhs)
        {
            os << '(' << c << ',' << d << ')';
            return os;
        }
#undef a
#undef b
#undef c
#undef d
    };

    template <typename T>
    Complex<T> I_V{0, 1};

    template <typename T, typename U>
    Complex<U> operator+(const U &lhs, const Complex<T> &rhs)
    {
        Complex<U> temp(lhs, 0);
        return temp + rhs;
    }

    template <typename T, typename U>
    Complex<U> operator-(const U &lhs, const Complex<T> &rhs)
    {
        Complex<U> temp(lhs, 0);
        return temp - rhs;
    }

    template <typename T, typename U>
    Complex<U> operator*(const U &lhs, const Complex<T> &rhs)
    {
        Complex<U> temp(lhs, 0);
        return temp * rhs;
    }

    template <typename T, typename U>
    Complex<U> operator/(const U &lhs, const Complex<T> &rhs)
    {
        Complex<U> temp(lhs, 0);
        return temp / rhs;
    }

    template <typename T>
    inline T abs(const Complex<T> &val)
    {
        return (T)sycl::sqrt(val.real() * val.real() + val.imag() * val.imag());
    }

    template <typename T>
    inline double arg(const Complex<T> &val)
    {
        return sycl::atan2<double, double>(val.imag(), val.real());
    }

    template <typename T>
    inline T norm(const Complex<T> &val)
    {
        return val.real() * val.real() + val.imag() * val.imag();
    }

    template <typename T>
    inline Complex<T> conj(const Complex<T> &val)
    {
        return Complex<T>(val.real(), -val.imag());
    }

    template <typename T>
    inline Complex<T> polar(const T &r, const T theta = 0)
    {
        if (r < 0)
        {
            return Complex<T>(0, 0);
        }
        return Complex<T>(r * sycl::cos(theta), r * sycl::sin(theta));
    }

    template <typename T>
    inline Complex<T> exp(const Complex<T> &z)
    {
        return Complex<T>(sycl::exp(z.real()) * sycl::cos(z.imag()), sycl::sin(z.imag()));
    }

    template <typename T>
    inline Complex<T> log(const Complex<T> &z)
    {
        T r = sycl::abs(z);
        T theta = sycl::arg(z);
        return Complex<T>(sycl::log(r), theta);
    }

    template <typename T>
    inline Complex<T> log10(const Complex<T> &z)
    {
        T r = sycl::abs(z);
        T theta = sycl::arg(z);
        return Complex<T>(sycl::log10(r), sycl::log10(theta));
    }

    template <typename T>
    inline Complex<T> pow(const Complex<T> &val, const T &power)
    {
        switch (power)
        {
        case 0:
            return Complex<T>(1, 0);
        case 1:
            return val;
        case 2:
            return val * val;
        }
        T r = sycl::abs(val);
        r = sycl::pow(r, power);
        T theta = sycl::arg(val);
        theta = power > 0 ? theta * power : theta / sycl::abs(power);
        return Complex<T>(r * sycl::cos(theta), r * sin(theta));
    }

    template <typename T>
    inline Complex<T> pow(const Complex<T> &val, const Complex<T> &power)
    {
        T r = sycl::abs(val);
        T theta = sycl::pow(r, power);
        T a = sycl::pow(syConstants::e, val.real() * sycl::log(r) - val.imag() * theta);
        T b = val.imag() * sycl::log(r) * val.real() * theta;
        return Complex<T>(a * sycl::cos(b), a * sycl::sin(b));
    }

    template <typename T>
    inline Complex<T> pow(const T &val, const Complex<T> &power)
    {
        T a = sycl::pow(val, power.real());
        T b = power.imag() * sycl::log(val);
        return Complex<T>(a * sycl::cos(b), a * sycl::sin(b));
    }

    template <typename T>
    inline Complex<T> sqrt(const Complex<T> &val)
    {
        T r = sycl::abs(val);
        T a = sycl::sqrt((r + val.real()) / 2);
        T b = sycl::sqrt((r - val.real()) / 2);
        return Complex<T>(a, val.imag() < 0 ? b : -b);
    }

    template <typename T>
    inline Complex<T> sin(const Complex<T> &val)
    {
        return Complex<T>(
            sycl::sin(val.real()) * sycl::cosh(val.imag()),
            sycl::cos(val.real()) * sycl::sinh(val.imag()));
    }

    template <typename T>
    inline Complex<T> asin(const Complex<T> &val)
    {
        Complex<T> i(0, 1);
        // -i * ln(i * z + sqrt(1 - z^2))
        return -i * sycl::log(i * val + sycl::sqrt(Complex<T>(1, 0) - sycl::pow(val, 2)));
    }

    template <typename T>
    inline Complex<T> cos(const Complex<T> &val)
    {
        return Complex<T>(
            sycl::cos(val.real()) * sycl::cosh(val.imag()),
            -1 * (sycl::sin(val.real()) * sycl::sinh(val.imag())));
    }

    template <typename T>
    inline Complex<T> acos(const Complex<T> &val)
    {
        Complex<T> one(1, 0);
        Complex<T> i(0, -1);
        // -i * ln(z + i * sqrt(1 - z^2))
        return i * sycl::log(val + i * sycl::sqrt(one - sycl::pow(val, 2)));
    }

    template <typename T>
    inline Complex<T> tan(const Complex<T> &val)
    {
        Complex<T> sin_val = sycl::sin(val);
        Complex<T> cos_val = sycl::cos(val);
        return sin_val / cos_val;
    }

    template <typename T>
    inline Complex<T> atan(const Complex<T> &val)
    {
        Complex<T> i(0, 1);
        // (i/2) * ln((i + z)/(i - z))
        return (i / T(2)) * sycl::log((i + val) / (i - val));
    }

    template <typename T>
    inline Complex<T> sinh(const Complex<T> &val)
    {
        return Complex<T>(
            sycl::sinh(val.real()) * sycl::cos(val.imag()),
            sycl::cosh(val.real()) * sycl::sin(val.imag()));
    }

    template <typename T>
    inline Complex<T> asinh(const Complex<T> &val)
    {
        // ln(z + sqrt(z^2 + 1))
        return sycl::log(val + sycl::sqrt(sycl::pow(val, 2) + Complex<T>(1, 0)));
    }

    template <typename T>
    inline Complex<T> cosh(const Complex<T> &val)
    {
        return Complex<T>(
            sycl::cosh(val.real()) * sycl::cos(val.imag()),
            sycl::sinh(val.real()) * sycl::sin(val.imag()));
    }

    template <typename T>
    inline Complex<T> acosh(const Complex<T> &val)
    {
        // // ln(z + sqrt(z + 1) * sqrt(z - 1))
        // return sycl::log(val + sycl::sqrt(val + Complex<T>(1, 0)) * sycl::sqrt(val - Complex<T>(1, 0)));
        // Alternative form: ln(z + sqrt(z^2 - 1))
        return sycl::log(val + sycl::sqrt(sycl::pow(val, 2) - Complex<T>(1, 0)));
    }

    template <typename T>
    inline Complex<T> tanh(const Complex<T> &val)
    {
        Complex<T> sinh_val = sinh(val);
        Complex<T> cosh_val = cosh(val);
        return sinh_val / cosh_val;
    }

    template <typename T>
    inline Complex<T> atanh(const Complex<T> &val)
    {
        // (1/2) * ln((1 + z)/(1 - z))
        return (Complex<T>(0.5, 0) * sycl::log((Complex<T>(1, 0) + val) / (Complex<T>(1, 0) - val)));
    }
}

#undef CONSTEXPR_IF_UT