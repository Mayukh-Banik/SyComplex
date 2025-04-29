#include <gtest/gtest.h>
#include "syComplex.hpp"
#include <complex>
#include <random>
#include <type_traits>

using TestTypes = ::testing::Types<float, double>;

template <typename T>
T generateRandomNumber()
{
    static std::mt19937 rng(std::random_device{}());
    if constexpr (std::is_integral_v<T>)
        return std::uniform_int_distribution<T>(-100, 100)(rng);
    else
        return std::uniform_real_distribution<T>(-100.0, 100.0)(rng);
}

template <typename T>
bool complexEquals(const std::complex<T> &a, const sycl::Complex<T> &b, T tol = static_cast<T>(1e-6))
{
    return (std::abs(a.real() - b.real()) <= tol) && (std::abs(a.imag() - b.imag()) <= tol);
}

template <typename T>
class ComplexComparisonTest : public ::testing::Test
{
protected:
    using ValueType = T;

    T r1, i1, r2, i2;
    std::complex<T> std_c1, std_c2;
    sycl::Complex<T> sycl_c1, sycl_c2;

    void SetUp() override
    {
        r1 = generateRandomNumber<T>();
        i1 = generateRandomNumber<T>();
        r2 = generateRandomNumber<T>();
        i2 = generateRandomNumber<T>();

        std_c1 = std::complex<T>(r1, i1);
        std_c2 = std::complex<T>(r2, i2);
        sycl_c1 = sycl::Complex<T>(r1, i1);
        sycl_c2 = sycl::Complex<T>(r2, i2);
    }
};

TYPED_TEST_SUITE(ComplexComparisonTest, TestTypes);

TYPED_TEST(ComplexComparisonTest, RandomInitialization)
{
    using T = typename TestFixture::ValueType;
    ASSERT_TRUE(complexEquals(this->std_c1, this->sycl_c1))
        << "Mismatch on initialization. "
        << "std::complex: (" << this->std_c1.real() << ", " << this->std_c1.imag()
        << "), sycl::Complex: (" << this->sycl_c1.real() << ", " << this->sycl_c1.imag() << ")";
}

#define PIPE_OUTPUT << "std::complex: (" << STD.real() << ", " << STD.imag()                                          \
                    << "), sycl::Complex: (" << SYCL.real() << ", " << SYCL.imag() << ")" << "\nRandom values used: " \
                    << "r1: " << this->r1 << ", i1: " << this->i1                                                     \
                    << ", r2: " << this->r2 << ", i2: " << this->i2;

TYPED_TEST(ComplexComparisonTest, AdditionCC)
{
    auto STD = this->std_c1 + this->std_c2;
    auto SYCL = this->sycl_c1 + this->sycl_c2;
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}

TYPED_TEST(ComplexComparisonTest, AdditionCSReal)
{
    auto STD = this->std_c1 + this->r2;
    auto SYCL = this->sycl_c1 + this->r2;
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}

TYPED_TEST(ComplexComparisonTest, SubtractionCC)
{
    auto STD = this->std_c1 - this->std_c2;
    auto SYCL = this->sycl_c1 - this->sycl_c2;
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}

TYPED_TEST(ComplexComparisonTest, SubtractionCSReal)
{
    auto STD = this->std_c1 - this->r2;
    auto SYCL = this->sycl_c1 - this->r2;
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}

TYPED_TEST(ComplexComparisonTest, MultiplicationCC)
{
    auto STD = this->std_c1 * this->std_c2;
    auto SYCL = this->sycl_c1 * this->sycl_c2;
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}

TYPED_TEST(ComplexComparisonTest, MultiplicationCSReal)
{
    auto STD = this->std_c1 * this->r2;
    auto SYCL = this->sycl_c1 * this->r2;
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}

TYPED_TEST(ComplexComparisonTest, DivisionCC)
{
    auto STD = this->std_c1 / this->std_c2;
    auto SYCL = this->sycl_c1 / this->sycl_c2;
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}

TYPED_TEST(ComplexComparisonTest, DivisionCSReal)
{
    auto STD = this->std_c1 / this->r2;
    auto SYCL = this->sycl_c1 / this->r2;
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}

TYPED_TEST(ComplexComparisonTest, AdditionItself)
{
    this->sycl_c1 += this->sycl_c2;
    this->std_c1 += this->std_c2;
    auto STD = this->std_c1;
    auto SYCL = this->sycl_c1;
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}

TYPED_TEST(ComplexComparisonTest, SubtractionItself)
{
    this->sycl_c1 -= this->sycl_c2;
    this->std_c1 -= this->std_c2;
    auto STD = this->std_c1;
    auto SYCL = this->sycl_c1;
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}

TYPED_TEST(ComplexComparisonTest, MultiplicationItself)
{
    this->sycl_c1 *= this->sycl_c2;
    this->std_c1 *= this->std_c2;
    auto STD = this->std_c1;
    auto SYCL = this->sycl_c1;
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}

TYPED_TEST(ComplexComparisonTest, DivisionItself)
{
    this->sycl_c1 /= this->sycl_c2;
    this->std_c1 /= this->std_c2;
    auto STD = this->std_c1;
    auto SYCL = this->sycl_c1;
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}

TYPED_TEST(ComplexComparisonTest, AdditionSCReal)
{
    auto STD = this->r2 + this->std_c1;
    auto SYCL = this->r2 + this->sycl_c1;
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}

TYPED_TEST(ComplexComparisonTest, SubtractionSCReal)
{
    auto STD = this->r2 - this->std_c1;
    auto SYCL = this->r2 - this->sycl_c1;
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}

TYPED_TEST(ComplexComparisonTest, MultiplicationSCReal)
{
    auto STD = this->r2 * this->std_c1;
    auto SYCL = this->r2 * this->sycl_c1;
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}

TYPED_TEST(ComplexComparisonTest, DivisionSCReal)
{
    auto STD = this->r2 / this->std_c1;
    auto SYCL = this->r2 / this->sycl_c1;
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}

TYPED_TEST(ComplexComparisonTest, ABS)
{
    auto STD = std::abs(this->std_c1);
    auto SYCL = sycl::abs(this->sycl_c1);
    ASSERT_TRUE(STD - SYCL <= 1e-06)
        << "STD: " << STD << " SYCL: " << SYCL << std::endl;
    ;
}

TYPED_TEST(ComplexComparisonTest, ARG)
{
    auto STD = std::arg(this->std_c1);
    auto SYCL = sycl::arg(this->sycl_c1);
    ASSERT_TRUE(STD - SYCL <= 1e-06)
        << "STD: " << STD << " SYCL: " << SYCL << std::endl;
    ;
}

TYPED_TEST(ComplexComparisonTest, NORM)
{
    auto STD = std::norm(this->std_c1);
    auto SYCL = sycl::norm(this->sycl_c1);
    ASSERT_TRUE(STD - SYCL <= 1e-06)
        << "STD: " << STD << " SYCL: " << SYCL << std::endl;
    ;
}

TYPED_TEST(ComplexComparisonTest, CONJ)
{
    auto STD = std::conj(this->std_c1);
    auto SYCL = sycl::conj(this->sycl_c1);
    ASSERT_TRUE(complexEquals(STD, SYCL))
        << "STD: " << STD << " SYCL: " << SYCL << std::endl;
    ;
}

TYPED_TEST(ComplexComparisonTest, LN)
{
    auto SYCL = sycl::log(this->sycl_c1);
    auto STD = std::log(this->std_c1);
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}

TYPED_TEST(ComplexComparisonTest, EXP)
{
    auto SYCL = sycl::exp(this->sycl_c1);
    auto STD = std::exp(this->std_c1);
    ASSERT_TRUE(complexEquals(STD, SYCL))
    PIPE_OUTPUT;
}