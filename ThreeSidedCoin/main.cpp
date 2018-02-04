#include <iostream>
#include <cmath>
#include <random>
#include <map>
#include <string>
#include "Vector.h"

const double PI = 3.141592653589793238462643383279502884197169399375105820974944592;
const double G = 9.80665; // m/s^2
const double DENSITY = 1225.0; // kg/m^3
const double COEFFICIENT_RESTITUTION = 1.0; // ND 
const double RADIUS = 0.03; // m = 3 cm
const double INITIAL_HEIGHT = 0.2; // m = 20 cm
const bool   VERBOSE = true;

void log(std::string message) {
	if (VERBOSE) {
		std::cout << message << std::endl;
	}
}

double radians(double degrees) {
	return degrees / 180.0 * PI;
}

double degrees(double radians) {
	return radians / PI * 180.0;
}

double calculateMass(double radius, double thickness) {
	// V = pi r^2 h
	const double volume = PI * radius * radius * thickness;
	// m = rho V
	const double mass = volume * DENSITY;

	return mass;
}

double calculateMomentOfInertial(double mass, double radius, double thickness) {
	const double mI = 1.0 / 12.0 * mass * (3 * radius * radius + thickness* thickness);
	return mI;
	// 1/2 m (3r^2 + h^2)
}

class ThreeSidedCoin {
private:
public:
	const double radius;
	const double theta;
	const double thickness;
	double centerToEdge;
	const double mass;
	const double momentI;

	double KE = 0;
	double PE = 0;
	double RE = 0;
	double omegaSign = 1.0; // either 1.0 or -1.0

	ThreeSidedCoin(double theta_) :
		radius (RADIUS),
		theta(theta_),
		thickness (radius * 2 * std::tan(theta)),
		centerToEdge (radius / std::cos(theta)),
		mass (calculateMass(radius, thickness)),
		momentI (calculateMomentOfInertial(mass, radius, thickness))
	{
		PE = G * INITIAL_HEIGHT * mass;
	}

	double getRandomAlpha() const {
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
		if (sign == 1) {
			angle *= -1.0;
		}
		return angle;
	}

	void logTotalEnergy() const {
		log(" PE:" + std::to_string(PE) + 
			" KE:" + std::to_string(KE) +
			" RE:" + std::to_string(RE) +
			" Total energy : " + std::to_string(PE + KE + RE));
	}

	void setPE(const double height) {
		PE = G * height * mass;
	}

	void setKE(const double velocity) {
		KE = 0.5 * mass * velocity * velocity;
	}

	void setRE(const double omega) {
		if (omega < 0) {
			omegaSign = -1.0;
		} else {
			omegaSign = 1.0;
		}

		RE = 0.5 * momentI * omega * omega;
	}

	double getHeight() const {
		return PE / (G * mass);
	}

	double getVelocity() const {
		return std::sqrt(2 * KE / mass);
	}

	double getOmega() const {
		return omegaSign * std::sqrt(2 * RE / momentI);
	}

	void collide() {
		// randomly choose the angle that the coin impacts at
		const double alpha = getRandomAlpha();

		logTotalEnergy();
		// position the coin at the required height for the edge to touch
		// calculate the new KE and PE that this would entail
		const double height = centerToEdge * std::cos(std::abs(alpha) - theta);
		const double newPE = G * height * mass;
		const double deltaPE = PE - newPE;
		KE += deltaPE;
		PE = newPE;
		logTotalEnergy();
		

		// Find impulse
		const Vector velocityCG (0.0, -getVelocity());
		const double R_angle = (alpha > 0.0) ? alpha - theta : alpha + theta;
		const Vector R = Vector::fromAngleAndMag(R_angle, centerToEdge); // In world coordinates
		const double omega = getOmega();
		const Vector velocity_contact = velocityCG + cross(omega, R);
		// Velocity along the normal
		const double v_n = velocity_contact.y;
		if (v_n > 0) {
			return; // The point is actually moving away; no collision
		}
		
		const Vector impulse(
			0.0,
			-(1 + COEFFICIENT_RESTITUTION) * v_n/
			(1/mass + (R.x*R.x)/momentI)
		);

		const double new_omega = omega + impulse.cross(R)/momentI;
		const Vector new_vel = velocityCG + impulse / mass;

		setKE(new_vel.y);
		setRE(new_omega);

		logTotalEnergy();
	}

	
};

int main() {

	ThreeSidedCoin c (radians(30));

	c.collide();

	return 0;
}

