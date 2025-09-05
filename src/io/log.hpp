// io/log.hpp
#pragma once
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <memory>

#ifndef SPDLOG_FUNCTION
#define SPDLOG_FUNCTION __func__
#endif

namespace circa {
namespace log {

// Create once, reuse forever, and set as default so SPDLOG_* are safe.
inline std::shared_ptr<spdlog::logger>
init_and_get(const std::string& name = "circa",
             const std::string& pattern = "[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %s:%# %! | %v",
             spdlog::level::level_enum level = spdlog::level::info) {
    static std::once_flag once;
    static std::shared_ptr<spdlog::logger> logger;

    std::call_once(once, [&] {
        if(auto existing = spdlog::get(name)) {
            logger = existing;
        } 
        else {
            logger = spdlog::stdout_color_mt(name);
        }
        spdlog::set_default_logger(logger);
    });

    if(logger) {
        logger->set_level(level);
        logger->set_pattern(pattern);
    }

    return logger;
}

}  // namespace log
}  // namespace circa

#define CIRCA_LOG(level, ...)                                                \
    do {                                                                     \
        auto lg = circa::log::init_and_get();                              \
        if(lg && lg->should_log(level))                                     \
            lg->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, \
                    level, __VA_ARGS__);                                     \
    } while(0)

#define CIRCA_TRACE(...) CIRCA_LOG(spdlog::level::trace, __VA_ARGS__)
#define CIRCA_DEBUG(...) CIRCA_LOG(spdlog::level::debug, __VA_ARGS__)
#define CIRCA_INFO(...) CIRCA_LOG(spdlog::level::info, __VA_ARGS__)
#define CIRCA_WARN(...) CIRCA_LOG(spdlog::level::warn, __VA_ARGS__)
#define CIRCA_ERROR(...) CIRCA_LOG(spdlog::level::err, __VA_ARGS__)
#define CIRCA_CRITICAL(...) CIRCA_LOG(spdlog::level::critical, __VA_ARGS__)
