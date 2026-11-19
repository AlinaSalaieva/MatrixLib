#include "matrixlib/Matrix.hpp"
#include "ExceptionFlow.hpp"

#include <gtest/gtest.h>

#include <limits>
#include <sstream>
#include <string_view>
#include <utility>
#include <vector>

namespace {
using namespace matrixlib;

constexpr double tolerance = 1e-10;
constexpr double maxDouble = std::numeric_limits<double>::max();
constexpr double infinity = std::numeric_limits<double>::infinity();
constexpr double nan = std::numeric_limits<double>::quiet_NaN();

void expectSameState(const Matrix& current, const Matrix& before) {
    EXPECT_EQ(current.rows(), before.rows());
    EXPECT_EQ(current.cols(), before.cols());
    EXPECT_EQ(current.data(), before.data());
}

TEST(MatrixConstruction, InitializesFlatRowMajorBuffer) {
    const Matrix matrix{{1, 2, 3}, {4, 5, 6}};
    EXPECT_EQ(matrix.rows(), 2U);
    EXPECT_EQ(matrix.cols(), 3U);
    EXPECT_EQ(matrix.data(), (std::vector<double>{1, 2, 3, 4, 5, 6}));
}

TEST(MatrixConstruction, CreatesMatrixAtDimensionLimit) {
    const Matrix matrix(128, 128, 2.0);
    EXPECT_EQ(matrix.data().size(), 128U * 128U);
    EXPECT_EQ(matrix.at(127, 127), 2.0);
}

TEST(MatrixConstruction, RejectsDimensionsOutsideContract) {
    EXPECT_THROW((Matrix(0, 2)), std::invalid_argument);
    EXPECT_THROW((Matrix(2, 0)), std::invalid_argument);
    EXPECT_THROW((Matrix(129, 1)), std::invalid_argument);
    EXPECT_THROW((Matrix(1, 129)), std::invalid_argument);
    EXPECT_THROW((Matrix(Matrix::Rows{})), std::invalid_argument);
}

TEST(MatrixConstruction, RejectsJaggedInput) {
    EXPECT_THROW((Matrix{{1, 2}, {3}}), MatrixShapeException);
}

TEST(MatrixConstruction, RejectsNonFiniteElements) {
    EXPECT_THROW((Matrix{{1, nan}}), InvalidElementException);
    EXPECT_THROW((Matrix(1, 1, infinity)), InvalidElementException);
}

TEST(MatrixAccess, ReadsAndWritesAnElement) {
    Matrix matrix{{1, 2}, {3, 4}};
    EXPECT_EQ(matrix.at(1, 0), 3.0);
    matrix.set(1, 0, 7.0);
    EXPECT_EQ(matrix.data(), (std::vector<double>{1, 2, 7, 4}));
}

TEST(MatrixAccess, RejectsOutOfRangeRead) {
    const Matrix matrix{{1, 2}};
    EXPECT_THROW(static_cast<void>(matrix.at(1, 0)), std::out_of_range);
    EXPECT_THROW(static_cast<void>(matrix.at(0, 2)), std::out_of_range);
}

TEST(MatrixAccess, RejectsOutOfRangeWriteWithoutMutation) {
    Matrix matrix{{1, 2}};
    const Matrix before = matrix;
    EXPECT_THROW(matrix.set(0, 2, 7.0), std::out_of_range);
    expectSameState(matrix, before);
}

TEST(MatrixAccess, RejectsNonFiniteWriteWithoutMutation) {
    Matrix matrix{{1, 2}};
    const Matrix before = matrix;
    EXPECT_THROW(matrix.set(0, 0, -infinity), InvalidElementException);
    expectSameState(matrix, before);
}

TEST(MatrixAlgorithms, FindsDeterminantWithPivoting) {
    const Matrix matrix{{0, 1, 0}, {1, 1, 0}, {0, 0, 1}};
    EXPECT_NEAR(matrix.determinant(), -1.0, tolerance);
}

TEST(MatrixAlgorithms, FindsFullRank) {
    const Matrix matrix{{0, 1, 0}, {1, 1, 0}, {0, 0, 1}};
    EXPECT_EQ(matrix.rank(), 3U);
}

TEST(MatrixAlgorithms, FindsRankOfRectangularMatrix) {
    const Matrix matrix{{0, 1, 2}, {0, 2, 4}, {0, 0, 1}, {0, 3, 7}};
    EXPECT_EQ(matrix.rank(), 2U);
}

TEST(MatrixAlgorithms, RejectsNonSquareInputForSquareOperations) {
    const Matrix matrix{{1, 2}};
    EXPECT_THROW(static_cast<void>(matrix.determinant()), MatrixShapeException);
    EXPECT_THROW(static_cast<void>(matrix.inverse()), MatrixShapeException);
    EXPECT_THROW(static_cast<void>(matrix.isSingular()), MatrixShapeException);
}

TEST(MatrixAlgorithms, DetectsSingularityAndRejectsInverse) {
    const Matrix matrix{{1, 2}, {2, 4}};
    EXPECT_TRUE(matrix.isSingular());
    EXPECT_NEAR(matrix.determinant(), 0.0, tolerance);
    EXPECT_THROW(static_cast<void>(matrix.inverse()), SingularMatrixException);
}

TEST(MatrixAlgorithms, InverseProductIsIdentity) {
    const Matrix matrix{{0, 2, 1}, {2, 1, 0}, {1, 0, 3}};
    EXPECT_FALSE(matrix.isSingular());
    const Matrix product = matrix * matrix.inverse();
    const Matrix reverseProduct = matrix.inverse() * matrix;
    const Matrix expected = Matrix::identity(3);
    for (std::size_t i = 0; i < expected.data().size(); ++i) {
        EXPECT_NEAR(product.data()[i], expected.data()[i], tolerance);
        EXPECT_NEAR(reverseProduct.data()[i], expected.data()[i], tolerance);
    }
}

TEST(MatrixAlgorithms, ChecksAbsolutePivotTolerance) {
    const Matrix matrix{{1e-14, 0}, {0, 1}};
    EXPECT_TRUE(matrix.isSingular());
    EXPECT_FALSE(matrix.isSingular(1e-15));
    EXPECT_THROW(static_cast<void>(matrix.rank(-1.0)), std::invalid_argument);
    EXPECT_THROW(static_cast<void>(matrix.inverse(nan)), std::invalid_argument);
}

TEST(MatrixAlgorithms, ReportsNonFiniteDeterminant) {
    const Matrix matrix{{maxDouble, 0}, {0, 2}};
    EXPECT_THROW(static_cast<void>(matrix.determinant()), std::overflow_error);
}

TEST(MatrixFunctions, AddsMatrices) {
    const Matrix a{{1, 2}, {3, 4}};
    const Matrix b{{2, 4}, {6, 8}};
    const Matrix expected{{3, 6}, {9, 12}};
    expectSameState(a + b, expected);
    EXPECT_THROW(static_cast<void>(a + Matrix(1, 2)), MatrixShapeException);
}

TEST(MatrixFunctions, MultipliesRectangularMatrices) {
    const Matrix a{{1, 2, 3}, {4, 5, 6}};
    const Matrix b{{7, 8}, {9, 10}, {11, 12}};
    const Matrix expected{{58, 64}, {139, 154}};
    expectSameState(a * b, expected);
    EXPECT_THROW(static_cast<void>(a * Matrix(2, 2)), MatrixShapeException);
}

TEST(MatrixFunctions, TransposesRectangularMatrix) {
    const Matrix a{{1, 2, 3}, {4, 5, 6}};
    const Matrix expected{{1, 4}, {2, 5}, {3, 6}};
    expectSameState(a.transpose(), expected);
}

TEST(MatrixFunctions, CreatesIdentityAndPrintsMatrix) {
    const Matrix matrix = Matrix::identity(2);
    EXPECT_EQ(matrix.data(), (std::vector<double>{1, 0, 0, 1}));
    std::ostringstream out;
    out << matrix;
    EXPECT_EQ(out.str(), "[1, 0]\n[0, 1]\n");
}

TEST(MatrixFunctions, RejectsNonFiniteArithmeticResults) {
    const Matrix matrix{{maxDouble}};
    EXPECT_THROW(static_cast<void>(matrix + matrix), std::overflow_error);
    EXPECT_THROW(static_cast<void>(matrix * Matrix{{2}}), std::overflow_error);
    EXPECT_EQ(matrix.at(0, 0), maxDouble);
}

TEST(MatrixStrong, BoundsFailurePreservesState) {
    Matrix matrix{{1, 2}, {3, 4}};
    const Matrix before = matrix;
    const Matrix::Rows block{{5, 6}, {7, 8}};
    EXPECT_THROW(matrix.replaceSubmatrix(1, 1, block), std::out_of_range);
    expectSameState(matrix, before);
}

TEST(MatrixStrong, HugeIndexFailurePreservesState) {
    Matrix matrix{{1, 2}, {3, 4}};
    const Matrix before = matrix;
    const auto huge = std::numeric_limits<std::size_t>::max();
    EXPECT_THROW(matrix.replaceSubmatrix(huge, 0, {{5}}), std::out_of_range);
    EXPECT_THROW(matrix.replaceSubmatrix(0, huge, {{5}}), std::out_of_range);
    expectSameState(matrix, before);
}

TEST(MatrixStrong, EmptyBlockFailurePreservesState) {
    Matrix matrix{{1, 2}, {3, 4}};
    const Matrix before = matrix;
    EXPECT_THROW(matrix.replaceSubmatrix(0, 0, {}), std::invalid_argument);
    expectSameState(matrix, before);
}

TEST(MatrixStrong, JaggedBlockFailurePreservesState) {
    Matrix matrix{{1, 2}, {3, 4}};
    const Matrix before = matrix;
    EXPECT_THROW(matrix.replaceSubmatrix(0, 0, {{5, 6}, {7}}), MatrixShapeException);
    expectSameState(matrix, before);
}

TEST(MatrixStrong, NonFiniteBlockFailurePreservesState) {
    Matrix matrix{{1, 2}, {3, 4}};
    const Matrix before = matrix;
    EXPECT_THROW(matrix.replaceSubmatrix(0, 0, {{5, 6}, {7, nan}}), InvalidElementException);
    EXPECT_THROW(matrix.replaceSubmatrix(0, 0, {{infinity}}), InvalidElementException);
    expectSameState(matrix, before);
}

TEST(MatrixStrong, PartlyReplacesMatrixWithoutTouchingOutsideBlock) {
    Matrix matrix{{1, 2}, {3, 4}};
    const Matrix::Rows block{{5}, {6}};
    const auto blockBefore = block;
    matrix.replaceSubmatrix(0, 1, block);
    const Matrix expected{{1, 5}, {3, 6}};
    expectSameState(matrix, expected);
    EXPECT_EQ(block, blockBefore);
}

TEST(MatrixStrong, ReplacesLastElementAtDimensionLimit) {
    Matrix matrix(128, 128, 1.0);
    matrix.replaceSubmatrix(127, 127, {{6}});
    Matrix expected(128, 128, 1.0);
    expected.set(127, 127, 6.0);
    expectSameState(matrix, expected);
}

TEST(MatrixStrong, ReplacesEntireMatrix) {
    Matrix matrix{{1, 2}, {3, 4}};
    matrix.replaceSubmatrix(0, 0, {{5, 6}, {7, 8}});
    const Matrix expected{{5, 6}, {7, 8}};
    expectSameState(matrix, expected);
}

TEST(MatrixCopySwap, CopyHasIndependentData) {
    const Matrix original{{1, 2}, {3, 4}};
    Matrix copy(original);
    copy.set(0, 0, 8.0);
    EXPECT_EQ(original.at(0, 0), 1.0);
    EXPECT_EQ(copy.at(0, 0), 8.0);
}

TEST(MatrixCopySwap, AssignmentAndSelfAssignmentPreserveInvariants) {
    Matrix matrix{{1}};
    const Matrix other{{1, 2, 3}, {4, 5, 6}};
    matrix = other;
    expectSameState(matrix, other);
    const auto& alias = matrix;
    matrix = alias;
    expectSameState(matrix, other);
}

TEST(MatrixCopySwap, SwapExchangesDimensionsAndDataWithoutThrowing) {
    Matrix a{{1, 2}};
    Matrix b{{3}, {4}};
    const Matrix beforeA = a;
    const Matrix beforeB = b;
    static_assert(noexcept(swap(a, b)));
    swap(a, b);
    expectSameState(a, beforeB);
    expectSameState(b, beforeA);
}

TEST(ExceptionFlow, ScopeLoggersUnwindInReverseOrder) {
    demo::Trace trace;
    EXPECT_THROW(demo::level1(trace), SingularMatrixException);
    const auto events = trace.events();
    const std::vector<std::string_view> expected{
        "enter level1", "enter level2", "enter level3", "leave level3",
        "level2 caught MatrixException and rethrows", "leave level2", "leave level1"
    };
    EXPECT_EQ((std::vector<std::string_view>{events.begin(), events.end()}), expected);
}

TEST(ExceptionFlow, RethrowPreservesDynamicTypeAndMessage) {
    demo::Trace trace;
    try {
        demo::level1(trace);
        FAIL() << "Expected SingularMatrixException";
    } catch (const SingularMatrixException& error) {
        EXPECT_STREQ(error.what(), "Cannot invert a singular matrix");
    } catch (...) {
        FAIL() << "Exception dynamic type was not preserved";
    }
}

}  // namespace
