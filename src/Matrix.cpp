#include "matrixlib/Matrix.hpp"

#include <algorithm>
#include <cmath>
#include <ostream>
#include <utility>

namespace matrixlib {
namespace {

std::size_t checkedArea(std::size_t rows, std::size_t cols) {
    if (rows == 0 || cols == 0 || rows > Matrix::maxDimension ||
        cols > Matrix::maxDimension) {
        throw std::invalid_argument("Matrix dimensions must be in [1, 128]");
    }
    return rows * cols;  // The dimension limits make this product safe.
}

double checkedElement(double value) {
    if (!std::isfinite(value)) {
        throw InvalidElementException("Matrix elements must be finite");
    }
    return value;
}

double checkedResult(double value) {
    if (!std::isfinite(value)) {
        throw std::overflow_error("A matrix computation produced a non-finite value");
    }
    return value;
}

void checkTolerance(double tolerance) {
    if (!std::isfinite(tolerance) || tolerance < 0.0) {
        throw std::invalid_argument("Tolerance must be finite and non-negative");
    }
}

template<class Input>
std::vector<double> flatten(const Input& input, std::size_t rows, std::size_t cols) {
    const auto area = checkedArea(rows, cols);
    for (const auto& row : input) {
        if (row.size() != cols) {
            throw MatrixShapeException("Matrix rows must have the same length");
        }
    }

    std::vector<double> result;
    result.reserve(area);
    for (const auto& row : input) {
        for (double value : row) {
            result.push_back(checkedElement(value));
        }
    }
    return result;
}

}  // namespace

Matrix::Matrix(std::size_t rows, std::size_t cols, double initialValue)
    : rows_(rows), cols_(cols) {
    const auto area = checkedArea(rows_, cols_);
    data_.assign(area, checkedElement(initialValue));
}

Matrix::Matrix(const Rows& values)
    : rows_(values.size()), cols_(values.empty() ? 0 : values.front().size()),
      data_(flatten(values, rows_, cols_)) {}

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> values)
    : rows_(values.size()), cols_(values.size() == 0 ? 0 : values.begin()->size()),
      data_(flatten(values, rows_, cols_)) {}

Matrix& Matrix::operator=(const Matrix& other) {
    if (this != &other) {
        Matrix copy(other);
        swap(copy);
    }
    return *this;
}

std::size_t Matrix::checkedIndex(std::size_t row, std::size_t col) const {
    if (row >= rows_ || col >= cols_) {
        throw std::out_of_range("Matrix index is outside its bounds");
    }
    return row * cols_ + col;
}

double Matrix::at(std::size_t row, std::size_t col) const {
    return data_[checkedIndex(row, col)];
}

void Matrix::set(std::size_t row, std::size_t col, double value) {
    const auto index = checkedIndex(row, col);
    data_[index] = checkedElement(value);
}

void Matrix::replaceSubmatrix(std::size_t row, std::size_t col, const Rows& values) {
#if defined(LAB7_INJECT_BUG) && LAB7_INJECT_BUG
    // Negative control only: this deliberately breaks the strong guarantee.
    std::fill(data_.begin(), data_.end(), 0.0);
#endif

    // Preparation: any exception leaves the original buffer untouched.
    const Matrix block(values);
    if (row >= rows_ || col >= cols_ ||
        block.rows_ > rows_ - row || block.cols_ > cols_ - col) {
        throw std::out_of_range("Replacement block is outside matrix bounds");
    }

    auto prepared = data_;
    for (std::size_t r = 0; r < block.rows_; ++r) {
        for (std::size_t c = 0; c < block.cols_; ++c) {
            prepared[(row + r) * cols_ + col + c] = block.data_[r * block.cols_ + c];
        }
    }

    // Commit: vector<double> with the default allocator can swap without throwing.
    static_assert(noexcept(data_.swap(prepared)));
    data_.swap(prepared);
}

void Matrix::swap(Matrix& other) noexcept {
    using std::swap;
    swap(rows_, other.rows_);
    swap(cols_, other.cols_);
    data_.swap(other.data_);
}

void Matrix::requireSquare() const {
    if (rows_ != cols_) {
        throw MatrixShapeException("Operation requires a square matrix");
    }
}

void Matrix::swapRows(std::size_t first, std::size_t second) noexcept {
    if (first != second) {
        for (std::size_t c = 0; c < cols_; ++c) {
            std::swap(data_[first * cols_ + c], data_[second * cols_ + c]);
        }
    }
}

std::size_t Matrix::pivotRow(std::size_t first, std::size_t col) const noexcept {
    auto pivot = first;
    for (std::size_t r = first + 1; r < rows_; ++r) {
        if (std::abs(data_[r * cols_ + col]) > std::abs(data_[pivot * cols_ + col])) {
            pivot = r;
        }
    }
    return pivot;
}

void Matrix::eliminateBelow(std::size_t row, std::size_t col) {
    for (std::size_t r = row + 1; r < rows_; ++r) {
        const auto factor = checkedResult(data_[r * cols_ + col] / data_[row * cols_ + col]);
        data_[r * cols_ + col] = 0.0;
        for (std::size_t c = col + 1; c < cols_; ++c) {
            const auto term = checkedResult(factor * data_[row * cols_ + c]);
            auto& value = data_[r * cols_ + c];
            value = checkedResult(value - term);
        }
    }
}

double Matrix::determinant(double tolerance) const {
    requireSquare();
    checkTolerance(tolerance);
    Matrix work(*this);
    double sign = 1.0;
    for (std::size_t c = 0; c < cols_; ++c) {
        const auto pivot = work.pivotRow(c, c);
        if (std::abs(work.data_[pivot * cols_ + c]) <= tolerance) {
            return 0.0;
        }
        if (pivot != c) {
            work.swapRows(c, pivot);
            sign = -sign;
        }
        work.eliminateBelow(c, c);
    }
    for (std::size_t i = 0; i < rows_; ++i) {
        sign = checkedResult(sign * work.data_[i * cols_ + i]);
    }
    return sign;
}

std::size_t Matrix::rank(double tolerance) const {
    checkTolerance(tolerance);
    Matrix work(*this);
    std::size_t row = 0;
    for (std::size_t c = 0; c < cols_ && row < rows_; ++c) {
        const auto pivot = work.pivotRow(row, c);
        if (std::abs(work.data_[pivot * cols_ + c]) <= tolerance) {
            continue;
        }
        work.swapRows(row, pivot);
        work.eliminateBelow(row, c);
        ++row;
    }
    return row;
}

bool Matrix::isSingular(double tolerance) const {
    requireSquare();
    return rank(tolerance) < rows_;
}

Matrix Matrix::inverse(double tolerance) const {
    requireSquare();
    checkTolerance(tolerance);
    Matrix work(*this);
    Matrix result = identity(rows_);

    for (std::size_t c = 0; c < cols_; ++c) {
        const auto pivot = work.pivotRow(c, c);
        const auto divisor = work.data_[pivot * cols_ + c];
        if (std::abs(divisor) <= tolerance) {
            throw SingularMatrixException("Cannot invert a singular matrix");
        }

        work.swapRows(c, pivot);
        result.swapRows(c, pivot);
        for (std::size_t j = 0; j < cols_; ++j) {
            // Earlier pivot columns in the current working row are already zero.
            if (j > c) {
                work.data_[c * cols_ + j] = checkedResult(work.data_[c * cols_ + j] / divisor);
            }
            result.data_[c * cols_ + j] = checkedResult(result.data_[c * cols_ + j] / divisor);
        }
        work.data_[c * cols_ + c] = 1.0;

        for (std::size_t r = 0; r < rows_; ++r) {
            if (r == c) {
                continue;
            }
            const auto factor = work.data_[r * cols_ + c];
            work.data_[r * cols_ + c] = 0.0;
            for (std::size_t j = 0; j < cols_; ++j) {
                if (j > c) {
                    const auto term = checkedResult(factor * work.data_[c * cols_ + j]);
                    auto& value = work.data_[r * cols_ + j];
                    value = checkedResult(value - term);
                }
                const auto term = checkedResult(factor * result.data_[c * cols_ + j]);
                auto& value = result.data_[r * cols_ + j];
                value = checkedResult(value - term);
            }
        }
    }
    return result;
}

Matrix Matrix::transpose() const {
    Matrix result(cols_, rows_);
    for (std::size_t r = 0; r < rows_; ++r) {
        for (std::size_t c = 0; c < cols_; ++c) {
            result.data_[c * rows_ + r] = data_[r * cols_ + c];
        }
    }
    return result;
}

Matrix Matrix::identity(std::size_t size) {
    Matrix result(size, size);
    for (std::size_t i = 0; i < size; ++i) {
        result.data_[i * size + i] = 1.0;
    }
    return result;
}

Matrix operator+(const Matrix& lhs, const Matrix& rhs) {
    if (lhs.rows_ != rhs.rows_ || lhs.cols_ != rhs.cols_) {
        throw MatrixShapeException("Addition requires identical matrix dimensions");
    }
    Matrix result(lhs.rows_, lhs.cols_);
    for (std::size_t i = 0; i < result.data_.size(); ++i) {
        result.data_[i] = checkedResult(lhs.data_[i] + rhs.data_[i]);
    }
    return result;
}

Matrix operator*(const Matrix& lhs, const Matrix& rhs) {
    if (lhs.cols_ != rhs.rows_) {
        throw MatrixShapeException("Incompatible dimensions for matrix multiplication");
    }
    Matrix result(lhs.rows_, rhs.cols_);
    for (std::size_t r = 0; r < lhs.rows_; ++r) {
        for (std::size_t c = 0; c < rhs.cols_; ++c) {
            double sum = 0.0;
            for (std::size_t k = 0; k < lhs.cols_; ++k) {
                const auto term = checkedResult(lhs.data_[r * lhs.cols_ + k] *
                                                rhs.data_[k * rhs.cols_ + c]);
                sum = checkedResult(sum + term);
            }
            result.data_[r * rhs.cols_ + c] = sum;
        }
    }
    return result;
}

std::ostream& operator<<(std::ostream& out, const Matrix& matrix) {
    for (std::size_t r = 0; r < matrix.rows(); ++r) {
        out << '[';
        for (std::size_t c = 0; c < matrix.cols(); ++c) {
            if (c > 0) {
                out << ", ";
            }
            out << matrix.at(r, c);
        }
        out << "]\n";
    }
    return out;
}

}  // namespace matrixlib
