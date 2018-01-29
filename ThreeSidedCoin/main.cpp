#include <iostream>
#include <cmath>

const double PI = 3.141592653589793238462643383279502884197169399375105820974944592;
const double DENSITY = 1225.0; // kg/m^3
const double RADIUS = 0.03; // m = 3 cm

double calculateMass(double radius, double height) {
	// V = pi r^2 h
	double volume = PI * radius * radius * height;
	// m = rho V
	double mass = volume * DENSITY;

	return mass;
}

double calculateMomentOfInertial(double mass, double radius, double height) {
	return 1 / 12 * mass * (3 * radius * radius + height* height);
	// 1/2 m (3r^2 + h^2)
}

class ThreeSidedCoim {
private:


public:
	const double radius;
	const double theta;
	const double height;
	const double mass;
	const double momentI;

	double KE = 0;
	double PE = 0;
	double RE = 0;

	ThreeSidedCoim(double theta_) :
		radius (RADIUS),
		theta(theta_),
		height (radius * 2 * std::tan(theta)),
		mass (calculateMass(radius, height)),
		momentI (calculateMomentOfInertial(mass, radius, height))
	{
	}


	void collide() {

	}

	
};

int main() {

	std::cout << "Hello world" << std::endl;

	return 0;
}

