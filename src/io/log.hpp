#pragma once

// ------------------------------------------------------------
// This silences a known **false-positive** — Wstringop-overflow —
// triggered by inlined std::atomic<bool>::operator bool() and
// __atomic_load_n in GCC's <atomic> headers under optimization
// or sanitizers. It's a compiler diagnostic issue, not a real bug.
//
// See GCC Bugzilla #113775: warning from __atomic_load_n as a
// false positive due to sanitizer-influenced optimization paths.
// Link: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=113775
// ------------------------------------------------------------
#if defined(__GNUC__) && !defined(__clang__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#if defined(__GNUC__) && !defined(__clang__)
#  pragma GCC diagnostic pop
#endif

#include <memory>

#ifndef SPDLOG_FUNCTION
#define SPDLOG_FUNCTION __func__
#endif

// Macros that pass source location so %! (func), %s (file), %# (line) work.
#define CIRCA_LOG(level, ...)                            \
    (spdlog::get("circa"))->log(                                               \
        spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, \
        level, __VA_ARGS__)

#define CIRCA_TRACE(...) CIRCA_LOG(spdlog::level::trace, __VA_ARGS__)
#define CIRCA_DEBUG(...) CIRCA_LOG(spdlog::level::debug, __VA_ARGS__)
#define CIRCA_INFO(...) CIRCA_LOG(spdlog::level::info, __VA_ARGS__)
#define CIRCA_WARN(...) CIRCA_LOG(spdlog::level::warn, __VA_ARGS__)
#define CIRCA_ERROR(...) CIRCA_LOG(spdlog::level::err, __VA_ARGS__)
#define CIRCA_CRITICAL(...) CIRCA_LOG(spdlog::level::critical, __VA_ARGS__)

namespace circa {
namespace log {

// Create or reuse a named console logger with pattern + level.
inline void init(const std::string& name = "circa",
                 const std::string& pattern = "[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %s:%# %! | %v",
                 spdlog::level::level_enum level = spdlog::level::info) {
    auto existing = spdlog::get(name);
    if (existing) {
        existing->set_level(level);
        spdlog::set_pattern(pattern);
    }
    auto logger = spdlog::stdout_color_mt(name);
    logger->set_level(level);
    spdlog::set_pattern(pattern);
}

}  // namespace log
}  // namespace circa
