#include "ExceptionFlow.hpp"

#include <cstdio>
#include <exception>
#include <iostream>
#include <stdexcept>

using matrixlib::Matrix;

int main() {
    std::cout << "Strong exception safety in MatrixLib\n\n";
    Matrix matrix{{1, 2}, {3, 4}};
    const auto before = matrix;
    std::cout << "Before an invalid replacement:\n" << matrix;

    bool expectedException = false;
    try {
        matrix.replaceSubmatrix(1, 1, {{5, 6}, {7, 8}});
    } catch (const std::out_of_range& error) {
        expectedException = true;
        std::cout << "Caught std::out_of_range: " << error.what() << '\n';
    }

    const bool unchanged = before.rows() == matrix.rows() && before.cols() == matrix.cols() &&
                           before.data() == matrix.data();
    std::cout << "After an invalid replacement:\n" << matrix;
    std::cout << "State unchanged: " << (unchanged ? "yes" : "NO") << '\n';

    Matrix result = before;
    result.replaceSubmatrix(0, 1, {{9}, {8}});
    std::cout << "\nValid replacement result:\n" << result;
    std::cout << "Determinant: " << result.determinant() << '\n';
    std::cout << "Rank: " << result.rank() << '\n';
    if (result.isSingular()) {
        std::cout << "Inverse unavailable: matrix is singular.\n";
    } else {
        std::cout << "Inverse:\n" << result.inverse();
    }

    std::cout << "\nRAII and exception propagation\n" << std::flush;
    matrixlib::demo::Trace trace;
    trace.output = stdout;
    bool correctDynamicType = false;
    try {
        matrixlib::demo::level1(trace);
    } catch (const matrixlib::SingularMatrixException&) {
        correctDynamicType = true;
        trace.add("main caught SingularMatrixException");
    } catch (const std::exception& error) {
        std::cerr << "Unexpected exception: " << error.what() << '\n';
    }
    return expectedException && unchanged && correctDynamicType ? 0 : 1;
}
