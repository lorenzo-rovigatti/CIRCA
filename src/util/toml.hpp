#pragma once

#include "../io/log.hpp"

#define TOML_EXCEPTIONS 0
#define TOML_ENABLE_FORMATTERS 0
// the two define's that follow are require to work around a known toml++ bug (see https://github.com/marzer/tomlplusplus/issues/213)
#define TOML_RETURN_BOOL_FROM_FOR_EACH_BROKEN 1
#define TOML_RETURN_BOOL_FROM_FOR_EACH_BROKEN_ACKNOWLEDGED 1
#include <toml++/toml.hpp>

#include <string_view>

namespace circa {

inline const toml::table* as_table_ptr(const toml::node_view<const toml::node>& nv) {
    if (!nv) return nullptr;
    return nv.as_table();
}

template <typename T, size_t N>
std::array<T, N> array_from_toml(const toml::array& a, const char* key) {
    if (a.size() != N) {
        throw std::runtime_error(fmt::format("Expected {} to have {} elements", key, N));
    }
    std::array<T, N> out{};
    for (size_t i = 0; i < N; ++i) {
        out[i] = a[i].value<T>().value_or(T{});
    }
    return out;
}

template <class T>
std::optional<T> value_or_die(const toml::table& tbl, std::string_view key_path) {
    // Resolve dotted path (e.g. "server.port")
    toml::node_view<const toml::node> nv = tbl.at_path(key_path);

    if (!nv) {
        CIRCA_CRITICAL(fmt::format("Missing key '{}'", key_path));
        throw std::runtime_error("");
    }

    // Try to extract as T (like node.value<T>())
    if (auto v = nv.value<T>()) {
        return *v;
    }

    CIRCA_CRITICAL(fmt::format("Key '{}' has incompatible type (expected something convertible to {})", key_path, typeid(T).name()));
    throw std::runtime_error("");

    return std::nullopt;
}

// Convenience: returns T or a default; logs when default is used (optional).
template <class T>
T value_or(const toml::table& tbl, std::string_view key_path, T default_value) {
    if(auto v = tbl[key_path].value<T>()) {
        return *v;
    }

    CIRCA_INFO("Using default for {} ({})", key_path, default_value);
    return default_value;
}

template <typename T>
T value_or(const toml::table* tp, std::string_view key_path, T default_value) {
    if(!tp) {
        return default_value;
    }

    return value_or<T>(*tp, key_path, default_value);
}

}  // namespace circa
