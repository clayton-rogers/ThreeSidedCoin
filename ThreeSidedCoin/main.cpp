#include "Vector.h"


#include <iostream>
#include <cmath>
#include <random>
#include <map>
#include <string>

// TODO remove
#include <thread>
#include <chrono>

// Actual constants
const double PI = 3.141592653589793238462643383279502884197169399375105820974944592;
const double G = 9.80665; // m/s^2
const double DENSITY = 1225.0; // kg/m^3

// Experimental setup
const double COEFFICIENT_RESTITUTION = 0.9; // ND 
const double RADIUS = 0.03; // m = 3 cm
const double INITIAL_HEIGHT = 0.2; // m = 20 cm
// NOTE: theta is set in the program body
const int    NUMBER_OF_THROWS = 1000;
using namespace std::chrono_literals;
const auto   SLEEP_TIME = 0.0s;

// Program constants
const bool   VERBOSE = false;

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

	// Critical energy is the potential energy when the coin
	// is at rest (linear and rotational) and balancing on its edge
	const double criticalEnergy;
	double KE = 0;
	double PE = 0;
	double RE = 0;
	double omegaSign = 1.0; // either 1.0 or -1.0
	double last_alpha = 0;

	enum class Side {
		Heads,
		Tails,
		Edge
	};

	ThreeSidedCoin(double theta_) :
		radius(RADIUS),
		theta(theta_),
		thickness(radius * 2 * std::tan(theta)),
		centerToEdge(radius / std::cos(theta)),
		mass(calculateMass(radius, thickness)),
		momentI(calculateMomentOfInertial(mass, radius, thickness)),
		criticalEnergy(G * centerToEdge * mass)
	{
		PE = G * INITIAL_HEIGHT * mass;
	}

	double getRandomAlpha() const {
		std::random_device rd;
		std::mt19937_64 mt(rd());
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
		if (VERBOSE) {
			log(" PE:" + std::to_string(PE) +
				" KE:" + std::to_string(KE) +
				" RE:" + std::to_string(RE) +
				" Total energy : " + std::to_string(getTotalEnergy()));
		}
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

	double getTotalEnergy() const {
		return PE + KE + RE;
	}

	void collide() {
		// randomly choose the angle that the coin impacts at
		const double alpha = getRandomAlpha();
		last_alpha = alpha;

		// position the coin at the required height for the edge to touch
		// calculate the new KE and PE that this would entail
		const double height = centerToEdge * std::cos(std::abs(alpha) - theta);
		const double newPE = G * height * mass;
		const double deltaPE = PE - newPE;

		// If we could not possibly borrow 
if (newPE > (KE + PE)) {
	log("Impossible alpha!");
	return;
}

KE += deltaPE;
PE = newPE;

// Find impulse
const Vector velocityCG(0.0, -getVelocity());
const double R_angle = (alpha > 0.0) ? alpha - theta : alpha + theta;
const Vector R = Vector::fromAngleAndMag(R_angle, centerToEdge); // In world coordinates
const double omega = getOmega();
const Vector velocity_contact = velocityCG + cross(omega, R);
// Velocity along the normal
const double v_n = velocity_contact.y;
if (v_n > 0) {
	log("Point is moving away from collision!");
	return; // The point is actually moving away; no collision
}

const Vector impulse(
	0.0,
	-(1 + COEFFICIENT_RESTITUTION) * v_n /
	(1 / mass + (R.x*R.x) / momentI)
);

const double new_omega = omega - impulse.cross(R) / momentI;
const Vector new_vel = velocityCG + impulse / mass;

setKE(new_vel.y);
setRE(new_omega);

logTotalEnergy();
	}

	Side getSideFromAlpha(double alpha) {
		if (alpha > theta) {
			return Side::Heads;
		} else if (alpha < (-theta)) {
			return Side::Tails;
		} else {
			return Side::Edge;
		}
	}

	Side throwCoin() {
		while (true) {
			// Do a collision
			std::this_thread::sleep_for(SLEEP_TIME);
			collide();

			// If the total energy is smaller than the critical energy
			// then the coin has no way of getting out of its current gravity well
			if (getTotalEnergy() < criticalEnergy) {
				log("No more energy!");
				return getSideFromAlpha(last_alpha);
			}

			// If the kinetic energy and potential energy are greater than 
			// the critical energy then the coin has enough speed to escape
			// the gravity well by speed alone, so the next collision is a normal one
			if ((KE + PE) > criticalEnergy) {
				log("Bounced clear!");
				continue;
			}

			// Otherwise the coin has enough energy to escape its gravity well but only
			// if it uses some of its rotational energy (i.e. the rotating edge hits and
			// pushes the coin upwards). This will then make if more likely to hit at 
			// certain angles. For now we do nothing special in this case.
			log("Needed to use rotational energy!");
		}
	}

};

int main() {



	double theta = 30.0;

	do {
		std::map<ThreeSidedCoin::Side, int> results;
		results[ThreeSidedCoin::Side::Heads] = 0;
		results[ThreeSidedCoin::Side::Tails] = 0;
		results[ThreeSidedCoin::Side::Edge] = 0;
		for (int i = 0; i < NUMBER_OF_THROWS; ++i) {
			ThreeSidedCoin c(radians(theta));
			ThreeSidedCoin::Side result = c.throwCoin();
			results.at(result) += 1;
		}

		std::cout << "Results:"
			<< " Heads: " + std::to_string(results.at(ThreeSidedCoin::Side::Heads))
			<< " Tails: " + std::to_string(results.at(ThreeSidedCoin::Side::Tails))
			<< " Edge:  " + std::to_string(results.at(ThreeSidedCoin::Side::Edge))
			<< " Theta: " + std::to_string(theta)
			<< std::endl;

		if (results.at(ThreeSidedCoin::Side::Tails) > results.at(ThreeSidedCoin::Side::Edge)) {
			theta += 0.1;
		} else {
			theta -= 0.1;
		}

	} while (true);

	return 0;
}
