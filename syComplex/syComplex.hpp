#pragma once

#include <sycl/sycl.hpp>

namespace syConstants
{
    const double PI = 3.14159265358979323846264338327950288419716939937510582097494459230781640628;
    const double PIsqrt = 1.77245385090551602729816748334114518279754945612238712821380778985291128459;
    const double PIcbrt = 1.46459188756152326302014252726379039173859685562793717435725593713839364979;
    const double PIsquared = 9.86960440108935861883449099987615113531369940724079062641334937622004482241;
    const double PIcubed = 31.0062766802998201754763150671013952022252885658851076941445381038063949174;
    // 1 / sqrt(pi)
    const double PIrsqrt = 0.56418958354775628694807945156077258584405062932899885684408572171064246844;
    // 2 / sqrt(pi)
    const double PIdblrsqrt = 2 * PIrsqrt;
    const double PIdouble = 2 * PI;
    const double PIhalf = PI / 2;
    const double PIthird = PI / 3;
    const double PIfourth = PI / 4;
    const double PIfifth = PI / 5;
    const double PIsixth = PI / 6;
    const double PIseventh = PI / 7;
    const double PIeight = PI / 8;
    const double PIninth = PI / 9;
    const double PItenth = PI / 10;
    const double PIeleventh = PI / 11;
    const double PItwelvth = PI / 12;

    const double e = 2.71828182845904523536028747135266249775724709369995957496696762772407663035;
    // const double esquared = e * e;
    // const double esqrt = e ^ 0.5;

}

namespace sycl
{
    template <typename T>
    struct Complex
    {
        T _real;
        T _imag;

        constexpr Complex(T real = 0, T imag = 0) : _real(real), _imag(imag) {}

        inline constexpr T real() const { return this->_real; }
        inline constexpr T imag() const { return this->_imag; }
        inline constexpr void real(T val) { this->_real = val; }
        inline constexpr void imag(T val) { this->_imag = val; }

        inline constexpr Complex<T> operator+(const Complex<T> &rhs)
        {
            return Complex<T>(this->_real + rhs._real, this->_imag + rhs._imag);
        }

        template <typename U>
        inline constexpr Complex<T> operator+(const U &rhs)
        {
            return Complex<T>(this->_real + (T)rhs, this->_imag);
        }

        inline constexpr Complex<T> operator-(const Complex<T> &rhs)
        {
            return Complex<T>(this->_real - rhs._real, this->_imag - rhs._imag);
        }

        inline constexpr Complex<T> operator-(const T &rhs)
        {
            return Complex<T>(this->_real - (T)rhs, this->_imag);
        }

        constexpr Complex<T> operator*(const Complex<T> &rhs)
        {
            return Complex<T>(this->_real * rhs._real - this->_imag * rhs._imag, this->_real * rhs._imag + this->_imag * rhs._real);
        }

        constexpr Complex<T> operator*(const T &rhs)
        {
            return Complex<T>(this->_real * rhs, this->_imag * rhs);
        }

        constexpr Complex<T> operator/(const Complex<T> &rhs)
        {
            T num1, num2, denom;
            num1 = this->_real * rhs._real + this->_imag * rhs._imag;
            num2 = this->_imag * rhs._real - this->_real * rhs._imag;
            denom = rhs._real * rhs._real + rhs._imag * rhs._imag;
            if (denom == 0)
            {
                return Complex<T>(0, 0);
            }
            else
            {
                return Complex<T>(num1 / denom, num2 / denom);
            }
        }

        constexpr Complex<T> operator/(const T &rhs)
        {
            return Complex<T>(this->_real / rhs, this->_imag / rhs);
        }

        constexpr Complex<T> operator+() const
        {
            return Complex<T>(+this->_real, +this->_imag);
        }

        constexpr Complex<T> operator-() const
        {
            return Complex<T>(-this->_real, -this->_imag);
        }

        inline constexpr bool operator==(const Complex<T> &rhs) const
        {
            return this->_real == rhs._real && this->_imag == rhs._imag;
        }

        inline constexpr bool operator==(T &rhs) const
        {
            return this->_real == rhs && this->_imag == 0;
        }

        inline constexpr bool operator!=(const Complex<T> &rhs) const
        {
            return !(*this == rhs);
        }

        template <typename U>
        inline constexpr bool operator!=(const Complex<U> &rhs) const
        {
            return !(*this == rhs);
        }

        template <typename U = T>
        inline constexpr bool operator!=(U &rhs) const
        {
            return !(*this == rhs);
        }

        inline constexpr Complex<T> &operator+=(const Complex<T> &rhs)
        {
            *this = *this + rhs;
            return *this;
        }

        inline constexpr Complex<T> &operator-=(const Complex<T> &rhs)
        {
            *this = *this - rhs;
            return *this;
        }

        inline constexpr Complex<T> &operator*=(const Complex<T> &rhs)
        {
            *this = *this * rhs;
            return *this;
        }

        inline constexpr Complex<T> &operator/=(const Complex<T> &rhs)
        {
            *this = *this / rhs;
            return *this;
        }

        friend std::ostream &operator<<(std::ostream &os, const Complex<T> &rhs)
        {
            os << '(' << rhs._real << ',' << rhs._imag << ')';
            return os;
        }
    };

    namespace literal
    {
        template <typename T>
        Complex<T> i{0, 1};
    }

    template <typename T>
    Complex<T> operator+(const T &lhs, const Complex<T> &rhs)
    {
        return Complex<T>(lhs + rhs.real(), rhs.imag());
    }

    template <typename T>
    Complex<T> operator-(const T &lhs, const Complex<T> &rhs)
    {
        return Complex<T>(lhs - rhs.real(), -rhs.imag());
    }

    template <typename T>
    Complex<T> operator*(const T &lhs, const Complex<T> &rhs)
    {
        return Complex<T>(lhs * rhs.real(), lhs * rhs.imag());
    }

    template <typename T>
    Complex<T> operator/(const T &lhs, const Complex<T> &rhs)
    {
        T a = lhs;
        T b = rhs.real();
        T c = rhs.imag();
        T ab = a * b;
        T ac = a * c;
        T b2c2 = b * b + c * c;
        return Complex<T>(ab / b2c2, -ac / b2c2);
    }

    template <typename T>
    inline T abs(const Complex<T> &val)
    {
        return (T)sycl::sqrt(val.real() * val.real() + val.imag() * val.imag());
    }

    template <typename T>
    inline T arg(const Complex<T> &val)
    {
        return sycl::atan2(val.imag(), val.real());
    }

    template <typename T>
    inline T norm(const Complex<T> &val)
    {
        return val.real() * val.real() + val.imag() * val.imag();
    }
    template <typename T>
    bool operator<(const Complex<T> &lhs, const Complex<T> &rhs)
    {
        return sycl::norm(lhs) < sycl::norm(rhs) ? true : false;
    }

    template <typename T>
    bool operator<=(const Complex<T> &lhs, const Complex<T> &rhs)
    {
        return sycl::norm(lhs) <= sycl::norm(rhs) ? true : false;
    }

    template <typename T>
    bool operator>(const Complex<T> &lhs, const Complex<T> &rhs)
    {
        return sycl::norm(lhs) < sycl::norm(rhs) ? true : false;
    }

    template <typename T>
    bool operator>=(const Complex<T> &lhs, const Complex<T> &rhs)
    {
        return sycl::norm(lhs) <= sycl::norm(rhs) ? true : false;
    }

    template <typename T>
    inline Complex<T> conj(const Complex<T> &val)
    {
        return Complex<T>(val.real(), -val.imag());
    }

    template <typename T>
    inline Complex<T> polar(const T &r, const T theta = 0)
    {
        return Complex<T>(r * sycl::cos(theta), r * sycl::sin(theta));
    }

    template <typename T>
    inline Complex<T> exp(const Complex<T> &z)
    {
        const T a = sycl::exp(z.real());
        return Complex<T>(a * sycl::cos(z.imag()), a * sycl::sin(z.imag()));
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
        return (Complex<T>(0.5, 0) * sycl::log((Complex<T>(1, 0) + val) / (Complex<T>(1, 0) - val)));
    }

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
    inline Complex<T> atan2(const Complex<T> &z)
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
    inline Complex<T> cbrt(const Complex<T> &val)
    {
        return sycl::pow(val, 1 / 3);
    }

    // template <typename T>
    // inline Complex<T> ceil(const Complex<T> &val)
    // {
    //     return val;
    // }

    template <typename T>
    inline Complex<T> copysign(const Complex<T> &lhs, const Complex<T> &rhs)
    {
        T theta1 = sycl::arg(lhs);
        T theta2 = sycl::arg(rhs);
        theta1 = sycl::copysign(theta1, theta2);
        return sycl::polar(sycl::abs(lhs), theta1);
    }

    template <typename T>
    inline Complex<T> cospi(const Complex<T> &val)
    {
        return sycl::cos(val * syConstants::PI);
    }

    template <typename T>
    inline Complex<T> erf(const Complex<T> &val)
    {
        auto factorial = [](int n) -> int
        {
            if (n <= 1)
            {
                return 1;
            }
            int result = 1;
            for (int i = 2; i <= n; ++i)
            {
                result *= i;
            }
            return result;
        };
        Complex<T> sum(0, 0);
        Complex<T> term = val;
        Complex<T> z_squared = val * val;
        for (int n = 0; n < 20; n++)
        {
            T coef = static_cast<T>(1.0) / (factorial(n) * (2 * n + 1));
            if (n > 0)
            {
                term = -term * z_squared;
            }
            sum += term * coef;
        }

        return sum * syConstants::PIdblrsqrt;
    }

    template <typename T>
    inline Complex<T> erfc(const Complex<T> &val)
    {
        return Complex<T>(1, 0) - sycl::erf(val);
    }

    template <typename T>
    inline Complex<T> exp2(const Complex<T> &val)
    {
        return sycl::pow(2, val);
    }

    template <typename T>
    inline Complex<T> exp10(const Complex<T> &val)
    {
        return sycl::pow(10, val);
    }

    template <typename T>
    inline Complex<T> expm1(const Complex<T> &val)
    {
        return sycl::exp(val) - Complex<T>(1, 0);
    }

    template <typename T>
    inline Complex<T> fabs(const Complex<T> &val)
    {
        return Complex<T>(sycl::abs(val.real()), sycl::abs(val.imag()));
    }

    template <typename T>
    inline Complex<T> fdim(const Complex<T> &lhs, const Complex<T> &rhs)
    {
        if (sycl::abs(lhs) > sycl::abs(rhs))
        {
            return Complex<T>();
        }
        return lhs - rhs;
    }

    template <typename T>
    inline Complex<T> fmax(const Complex<T> &lhs, const Complex<T> &rhs)
    {
        return sycl::abs(lhs) < sycl::abs(rhs) ? rhs : lhs;
    }

    template <typename T>
    inline Complex<T> fmin(const Complex<T> &lhs, const Complex<T> &rhs)
    {
        return sycl::abs(lhs) > sycl::abs(rhs) ? rhs : lhs;
    }

    template <typename T>
    inline Complex<T> hypot(const Complex<T> &lhs, const Complex<T> &rhs)
    {
        return sycl::sqrt(lhs * lhs + rhs * rhs);
    }

    // do ilogb later

    template <typename T>
    inline Complex<T> ldexp(const Complex<T> &lhs, int k)
    {
        return lhs * sycl::pow(2, k);
    }

    // do lgamma later

    // do lgamma_r later

    template <typename T>
    inline Complex<T> log2(const Complex<T> &z)
    {
        T r = sycl::abs(z);
        T theta = sycl::arg(z);
        return Complex<T>(sycl::log2(r), sycl::log2(theta));
    }

    template <typename T>
    inline Complex<T> log1p(const Complex<T> &z)
    {
        T r = sycl::abs(z + Complex<T>(1, 0));
        T theta = sycl::arg(z + Complex<T>(1, 0));
        return Complex<T>(sycl::log(r), sycl::log2(theta));
    }

    // do logb later

    template <typename T>
    Complex<T> mad(const Complex<T> &a, const Complex<T> &b, const Complex<T> &c)
    {
        return a * b + c;
    }

    template <typename T>
    Complex<T> maxmag(const Complex<T> &a, const Complex<T> &b)
    {
        T bb = sycl::norm(a);
        T cc = sycl::norm(b);
        if (bb > cc)
        {
            return a;
        }
        else if (cc > bb)
        {
            return b;
        }
        else
        {
            return sycl::fmax(a, b);
        }
    }

    template <typename T>
    Complex<T> minmag(const Complex<T> &a, const Complex<T> &b)
    {
        T bb = sycl::norm(a);
        T cc = sycl::norm(b);
        if (bb < cc)
        {
            return a;
        }
        else if (cc < bb)
        {
            return b;
        }
        else
        {
            return sycl::fmin(a, b);
        }
    }

    // do modf later

    template <typename T>
    Complex<T> nextafter(const Complex<T> &a, const Complex<T> &b)
    {
        return Complex<T>(sycl::nextafter(a.real(), b.real()), sycl::nextafter(a.imag(), b.imag()));
    }

    template <typename T>
    Complex<T> pown(const Complex<T> &a, int power)
    {
        return sycl::pow(a, power);
    }

    template <typename T>
    Complex<T> powr(const Complex<T> &a, T power)
    {
        if (power < 0)
        {
            return Complex<T>();
        }
        return sycl::pow(a, power);
    }

    template <typename T>
    Complex<T> rootn(const Complex<T> &a, int power)
    {
        double v = 1 / power;
        return sycl::pow(a, power);
    }

    template <typename T>
    Complex<T> rsqrt(const Complex<T> &a)
    {
        return 1 / sycl::sqrt(a);
    }

    template <typename T>
    Complex<T> sinpi(const Complex<T> &a)
    {
        return sycl::sin(a * syConstants::PI);
    }

    template <typename T>
    Complex<T> tanpi(const Complex<T> &a)
    {
        return sycl::tan(a * syConstants::PI);
    }

    template <typename T>
    Complex<T> tgamma(const Complex<T> &z)
    {
        // alpha        beta        theta                           omega
        // sqrt(2pi) * sqrt(1/z) * [ (1 / (12z - 1/10z)) + z] ^ z * e^-z
        const T alpha = sycl::sqrt(2 * syConstants::PI);
        const Complex<T> beta = sycl::rsqrt(z);
        const Complex<T> theta = sycl::pow(((1 / (12 * z - (1 / (10 * z)))) + z), z);
        const Complex<T> omega = sycl::exp(-z);
        return alpha * beta * theta * omega;
    }

    template <typename T>
    T lgamma(const Complex<T> &z)
    {
        return sycl::log(sycl::abs(sycl::tgamma(z)));
    }

    template <typename T>
    struct plus
    {
        sycl::Complex<T> operator()(const sycl::Complex<T> &x, const sycl::Complex<T> &y) const
        {
            return x + y;
        }
    };

    template <typename T>
    struct multiplies
    {
        sycl::Complex<T> operator()(const sycl::Complex<T> &x, const sycl::Complex<T> &y) const
        {
            return x * y;
        }
    };

    template <typename T>
    struct bit_and
    {
        sycl::Complex<T> operator()(const sycl::Complex<T> &x, const sycl::Complex<T> &y) const
        {
            return Complex<T>(x.real() & y.real(), x.imag() & y.imag());
        }
    };

    template <typename T>
    struct bit_or
    {
        sycl::Complex<T> operator()(const sycl::Complex<T> &x, const sycl::Complex<T> &y) const
        {
            return Complex<T>(x.real() | y.real(), x.imag() | y.imag());
        }
    };

    template <typename T>
    struct bit_xor
    {
        sycl::Complex<T> operator()(const sycl::Complex<T> &x, const sycl::Complex<T> &y) const
        {
            return Complex<T>(x.real() ^ y.real(), x.imag() ^ y.imag());
        }
    };

    template <typename T>
    struct logical_and
    {
        sycl::Complex<T> operator()(const sycl::Complex<T> &x, const sycl::Complex<T> &y) const
        {
            return Complex<T>(x.real() && y.real(), x.imag() && y.imag());
        }
    };

    template <typename T>
    struct logical_or
    {
        sycl::Complex<T> operator()(const sycl::Complex<T> &x, const sycl::Complex<T> &y) const
        {
            return Complex<T>(x.real() || y.real(), x.imag() || y.imag());
        }
    };

    template <typename T>
    struct minimum
    {
        sycl::Complex<T> operator()(const sycl::Complex<T> &x, const sycl::Complex<T> &y) const
        {
            return (std::abs(x) < std::abs(y)) ? x : y;
        }
    };

    template <typename T>
    struct maximum
    {
        sycl::Complex<T> operator()(const sycl::Complex<T> &x, const sycl::Complex<T> &y) const
        {
            return (std::abs(x) > std::abs(y)) ? x : y;
        }
    };

    using Complex128 = Complex<double>;
    using Complex64 = Complex<float>;
    using Complex32 = Complex<sycl::half>;
}
