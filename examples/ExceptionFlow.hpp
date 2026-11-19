#pragma once

#include "matrixlib/Matrix.hpp"

#include <array>
#include <cstdio>
#include <memory>
#include <span>
#include <string_view>

namespace matrixlib::demo {

class Trace {
public:
    void add(std::string_view event) noexcept {
        if (size_ < events_.size()) {
            events_[size_++] = event;
        }
        if (output != nullptr) {
            std::fwrite(event.data(), 1, event.size(), output);
            std::fputc('\n', output);
        }
    }

    [[nodiscard]] std::span<const std::string_view> events() const noexcept {
        return {events_.data(), size_};
    }

    void recordCause(const char* message) const noexcept {
        if (output != nullptr) {
            std::fprintf(output, "Cause recorded by level2: %s\n", message);
        }
    }

    std::FILE* output = nullptr;

private:
    std::array<std::string_view, 10> events_{};
    std::size_t size_ = 0;
};

class ScopeLogger {
public:
    ScopeLogger(Trace& trace, std::string_view enter, std::string_view leave) noexcept
        : trace_(trace), leave_(leave) {
        trace_.add(enter);
    }

    ScopeLogger(const ScopeLogger&) = delete;
    ScopeLogger& operator=(const ScopeLogger&) = delete;

    ~ScopeLogger() noexcept { trace_.add(leave_); }

private:
    Trace& trace_;
    std::string_view leave_;
};

inline void level3(Trace& trace) {
    const ScopeLogger scope(trace, "enter level3", "leave level3");
    const auto matrix = std::make_unique<Matrix>(Matrix::Rows{{1, 2}, {2, 4}});
    static_cast<void>(matrix->inverse());
}

inline void level2(Trace& trace) {
    const ScopeLogger scope(trace, "enter level2", "leave level2");
    try {
        level3(trace);
    } catch (const MatrixException& error) {
        trace.add("level2 caught MatrixException and rethrows");
        trace.recordCause(error.what());
        throw;  // Preserve the original dynamic type.
    }
}

inline void level1(Trace& trace) {
    const ScopeLogger scope(trace, "enter level1", "leave level1");
    level2(trace);
}

}  // namespace matrixlib::demo
