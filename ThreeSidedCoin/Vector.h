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

// Assumes first vector only has a z component and the second
// has x,y. Then the resultant cross product will be in the 
// x,y plane.
// https://en.wikipedia.org/wiki/Cross_product#Coordinate_notation
Vector cross(double z_value, const Vector& v) {
	//double u1 = 0.0, u2 = 0.0, u3 = z_value;
	//double v1 = v.x, v2 = v.y, v3 = 0;

	//double new_x = u2*v3 - u3*v2;
	//double new_y = u3*v1 - u1*v3;
	//double new_z = u1*v2 - u2*v1; // Always zero

	//return Vector(new_x, new_y);

	return Vector(
		(-z_value*v.y),
		(z_value*v.x)
	);
}