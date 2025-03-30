#pragma once

#include <stdexcept>

template <typename T>

class Vector2d {

	T x;
	T y;

public: 

	using value_type = T;

	// Ensure T is a numeric type.
    static_assert(std::is_arithmetic<T>::value, "T must be a numeric type");

	Vector2d(): x(static_cast<T>(0)), y(static_cast<T>(0)) {} //Standardconstructor
	
	Vector2d(T v1, T v2) : x( v1 ), y( v2 )  {} //Constructor

	void set(T v1, T v2) {

		x = v1;
		y = v2;
	}

	T operator[] (size_t index ) const{
		
		if (index == 0) {
			return x;
		}
		else if (index == 1) {
			return y;
		}
		else
			throw std::out_of_range("Invalid index!");
	}

	// Operators for Math Operations (Exercise: 2b)
	Vector2d<T> operator+(const Vector2d<T> vector2) const { // Addition: Vector + Vector
		return Vector2d<T>(x + vector2.x, y + vector2.y);
	}

	Vector2d<T> operator-(const Vector2d<T> vector2) const { // Subtraction: Vector - Vector
		return Vector2d<T>(x - vector2.x, y - vector2.y);
	}

	Vector2d<T> operator*(T scalar) const { // Multiplication: Scalar * Vector (Not the dot product!)
		return Vector2d<T>(x * scalar, y * scalar);
	}

	Vector2d<T> operator/(T scalar) const { // Division: Vector / Scalar
		if (scalar == 0) {
			throw std::invalid_argument("Cannot divide by zero!");
		}
		return Vector2d<T>(x / scalar, y / scalar);
	}

	bool operator==(const Vector2d<T> vector2) const { // Equality
		return x == vector2.x && y == vector2.y;
	}

};
