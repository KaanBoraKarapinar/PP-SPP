#pragma once

#include <stdexcept>
#include <tuple>
#include <cstdint>

template <typename T> class Vector2d{
public:

    Vector2d(){
        contents = std::tuple<T, T>{T(), T()};
    }

    Vector2d(T arg_first, T arg_second){
        contents = std::tuple<T, T>{arg_first, arg_second};
    }

    void set(T arg_first, T arg_second){
        std::get<0>(contents) = arg_first;
        std::get<1>(contents) = arg_second;
    }

    bool operator==(Vector2d other) const {
        return (std::get<0>(contents) == other[0]) && (std::get<1>(contents) == other[1]);
    }

    Vector2d operator+(Vector2d other) const {
        return Vector2d(std::get<0>(contents) + other[0], std::get<1>(contents) + other[1]);
    }

    Vector2d operator-(Vector2d other) const {
        return Vector2d(std::get<0>(contents) - other[0], std::get<1>(contents) - other[1]);
    }

    Vector2d operator*(T scalar) const {
        return Vector2d(std::get<0>(contents) * scalar, std::get<1>(contents) * scalar);
    }

    Vector2d operator/(T scalar) const {
        return Vector2d(std::get<0>(contents) / scalar, std::get<1>(contents) / scalar);
    }

    T operator[](std::int32_t position) const {
        if (position == 0){
            return std::get<0>(contents);
        }
        if (position == 1){
            return std::get<1>(contents);
        }
        throw std::invalid_argument("Out of bounds access to Vector2d");
    }

private:
    std::tuple<T, T> contents;
};