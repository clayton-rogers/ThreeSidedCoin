#pragma once
#include <cmath>
#include <iostream>

class Vector {

public:
	double x, y;

	Vector(double x_ = 0.0, double y_ = 0.0) :
		x(x_), y(y_)
	{}

	Vector(const Vector& v) :
		x(v.x), y(v.y)
	{}



	Vector operator+ (const Vector& v) const {
		return Vector(x + v.x, y + v.y);
	}

	Vector flip() const {
		return Vector(-x, -y);
	}

	Vector operator- (const Vector& v) const {
		return *this + v.flip();
	}

	double mag() const {
		return std::sqrt(x*x + y*y);
	}

	double dot(const Vector& v) const {
		return x*v.x + y*v.y;
	}

	double cross(const Vector& v) const {
		return x*v.y - y*v.x;
	}

	Vector normal() const {
		return Vector(y, -x);
	}

	Vector rotate(double alpha) const {
		double new_x = x * std::cos(alpha) + y * std::sin(alpha);
		double new_y = -x * std::sin(alpha) + y * std::cos(alpha);

		return Vector(new_x, new_y);
	}

	static Vector fromAngleAndMag(double angle, double mag) {
		Vector v(0.0, mag);
		return v.rotate(angle);
	}

	void print() const {
		std::cout << "Vector: " << x << " " << y << std::endl;
	}
};

Vector operator* (double value, const Vector& vector) {
	return Vector(vector.x * value, vector.y * value);
}