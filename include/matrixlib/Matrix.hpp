#pragma once

#include <cstddef>
#include <initializer_list>
#include <iosfwd>
#include <stdexcept>
#include <vector>

namespace matrixlib {

class MatrixException : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class MatrixShapeException : public MatrixException {
public:
    using MatrixException::MatrixException;
};

class InvalidElementException : public MatrixException {
public:
    using MatrixException::MatrixException;
};

class SingularMatrixException : public MatrixException {
public:
    using MatrixException::MatrixException;
};

class Matrix {
public:
    using Row = std::vector<double>;
    using Rows = std::vector<Row>;
    static constexpr std::size_t maxDimension = 128;
    static constexpr double defaultTolerance = 1e-12;

    Matrix(std::size_t rows, std::size_t cols, double initialValue = 0.0);
    explicit Matrix(const Rows& values);
    Matrix(std::initializer_list<std::initializer_list<double>> values);

    Matrix(const Matrix&) = default;
    Matrix& operator=(const Matrix& other);

    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t cols() const noexcept { return cols_; }
    [[nodiscard]] const std::vector<double>& data() const noexcept { return data_; }

    [[nodiscard]] double at(std::size_t row, std::size_t col) const;
    void set(std::size_t row, std::size_t col, double value);

    // Zero-based coordinates. On failure, the matrix is unchanged.
    void replaceSubmatrix(std::size_t row, std::size_t col, const Rows& values);

    // Pivot tests use an absolute tolerance chosen by the caller.
    [[nodiscard]] double determinant(double tolerance = defaultTolerance) const;
    [[nodiscard]] std::size_t rank(double tolerance = defaultTolerance) const;
    [[nodiscard]] bool isSingular(double tolerance = defaultTolerance) const;
    [[nodiscard]] Matrix inverse(double tolerance = defaultTolerance) const;
    [[nodiscard]] Matrix transpose() const;
    [[nodiscard]] static Matrix identity(std::size_t size);

    void swap(Matrix& other) noexcept;

    friend void swap(Matrix& lhs, Matrix& rhs) noexcept { lhs.swap(rhs); }
    friend Matrix operator+(const Matrix& lhs, const Matrix& rhs);
    friend Matrix operator*(const Matrix& lhs, const Matrix& rhs);

private:
    std::size_t rows_;
    std::size_t cols_;
    std::vector<double> data_;

    [[nodiscard]] std::size_t checkedIndex(std::size_t row, std::size_t col) const;
    void requireSquare() const;
    void swapRows(std::size_t first, std::size_t second) noexcept;
    [[nodiscard]] std::size_t pivotRow(std::size_t first, std::size_t col) const noexcept;
    void eliminateBelow(std::size_t row, std::size_t col);
};

Matrix operator+(const Matrix& lhs, const Matrix& rhs);
Matrix operator*(const Matrix& lhs, const Matrix& rhs);
std::ostream& operator<<(std::ostream& out, const Matrix& matrix);

}  // namespace matrixlib
