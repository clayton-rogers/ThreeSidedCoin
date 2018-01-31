#include <iostream>
#include <cmath>
#include <random>

const double PI = 3.141592653589793238462643383279502884197169399375105820974944592;
const double G = 9.80665; // m/s^2
const double DENSITY = 1225.0; // kg/m^3
const double RADIUS = 0.03; // m = 3 cm
const double INITIAL_HEIGHT = 0.2; // m = 20 cm
const bool VERBOSE = true;

double radians(double degrees) {
	return degrees / 180.0 * PI;
}

double degrees(double radians) {
	return radians / PI * 180.0;
}

double calculateMass(double radius, double height) {
	// V = pi r^2 h
	double volume = PI * radius * radius * height;
	// m = rho V
	double mass = volume * DENSITY;

	return mass;
}

double calculateMomentOfInertial(double mass, double radius, double height) {
	double mI = 1.0 / 12.0 * mass * (3 * radius * radius + height* height);
	return mI;
	// 1/2 m (3r^2 + h^2)
}

class ThreeSidedCoin {
private:
public:
	const double radius;
	const double theta;
	const double height;
	double centerToEdge;
	const double mass;
	const double momentI;

	double KE = 0;
	double PE = 0;
	double RE = 0;

	ThreeSidedCoin(double theta_) :
		radius (RADIUS),
		theta(theta_),
		height (radius * 2 * std::tan(theta)),
		centerToEdge (radius / std::cos(theta)),
		mass (calculateMass(radius, height)),
		momentI (calculateMomentOfInertial(mass, radius, height))
	{
		PE = G * INITIAL_HEIGHT * mass;
	}

	double getRandomAlpha() {
		std::random_device rd;
		std::mt19937 mt(rd());
		std::uniform_real_distribution<double> dist(0.0, 1.0);
		double x = dist(mt);
		double angle = std::asin(x);
		// probability distribution is cos(theta) , 0 .. 90
		// cumulative is therefore int(cos(theta)) = sin(theta)
		// to get a value then you take the inverse: asin(x) = theta

		// angle can actually be plus or minus this amout
		std::uniform_int_distribution<> intdist(0, 1);
		int sign = intdist(mt);
		if (sign = 1) {
			angle *= -1.0;
		}
		return angle;
	}

	double getHeight(double alpha) {
		// TODO 
		return RADIUS;
	}

	void collide() {
		// randomly choose the angle that the coin impacts at
		double alpha = getRandomAlpha();

		// position the coin at the required height for the edge to touch
		// calculate the new KE and PE that this would entail
		double height = centerToEdge * std::cos(alpha - theta);
		double newPE = G * height * mass;
		double deltaPE = PE - newPE;
		KE += deltaPE;
		PE = newPE;

		// for now, elastic collision

	}

	
};

int main() {

	ThreeSidedCoin c (radians(30));

	c.collide();

	return 0;
}

