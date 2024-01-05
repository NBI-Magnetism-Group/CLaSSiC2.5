#include "integrator.hpp"
#include <iostream>
#include <cmath>
hey
void Integrator::calculateEffectiveField(std::vector<std::vector<int>> &neighbours, std::vector<double> &spin)
{
	/*
	Calculates the effective by combining: external field, anisotropy and nearest neighbour interaction
	*/
	double jx;
	double jy;
	double jz;

	for (int i = 0; i < constants::nAtoms; i++)
	{
		// External field
		memcpy(&effectiveField[3 * i], &constants::magneticField, sizeof(constants::magneticField)); // set to external field

		// Anisotropy
		effectiveField[3 * i + 2] += constants::anisotropyAxis * spin[3 * i + 2];
		effectiveField[3 * i + 2] -= constants::anisotropyPlane * spin[3 * i + 2];

		jx = 0;
		jy = 0;
		jz = 0;
		// Nearest neighbours
		for (int j : neighbours[i])
		{
			jx += spin[3 * j];	
			jy += spin[3 * j + 1];
			jz += spin[3 * j + 2];
		}
		effectiveField[3 * i]     += constants::exchangePrefactor * jx;	 // x
		effectiveField[3 * i + 1] += constants::exchangePrefactor * jy; // y
		effectiveField[3 * i + 2] += constants::exchangePrefactor * jz; // z

		//Stabilizer for kagome lattice by adding anisotropy (1 T) alligned with GS
		if (constants::geometry==4 && constants::stabilize){
			effectiveField[3*i] += stabilizerField[3*i];
			effectiveField[3*i+1] += stabilizerField[3*i+1];
		}
	}
}

void Integrator::evaluate(std::vector<std::vector<int>> &neighbours, std::vector<double> &spin)
{
	/*
	function evaluation
	*/
	calculateEffectiveField(neighbours, spin);
	for (int i = 0; i < constants::nAtoms; i++)
	{
		rk[3 * i] = -constants::gamma * constants::dt * (spin[3 * i + 1] * effectiveField[3 * i + 2] - spin[3 * i + 2] * effectiveField[3 * i + 1]
		- constants::lambda * 
		 (effectiveField[3 * i + 1] * spin[3 * i] * spin[3 * i + 1]
		  - effectiveField[3 * i] * spin[3 * i + 1] * spin[3 * i + 1]
		  + effectiveField[3 * i + 2] * spin[3 * i] * spin[3 * i + 2]
		   - effectiveField[3 * i] * spin[3 * i + 2] * spin[3 * i + 2]));
		
		rk[3 * i + 1] = -constants::gamma * constants::dt * (spin[3 * i + 2] * effectiveField[3 * i] - spin[3 * i] * effectiveField[3 * i + 2]
		- constants::lambda * 
		 (-effectiveField[3 * i + 1] * spin[3 * i] * spin[3 * i]
		  + effectiveField[3 * i] * spin[3 * i] * spin[3 * i + 1]
		  + effectiveField[3 * i + 2] * spin[3 * i + 1] * spin[3 * i + 2]
		   - effectiveField[3 * i + 1] * spin[3 * i + 2] * spin[3 * i + 2]));


		rk[3 * i + 2] = -constants::gamma * constants::dt * (spin[3 * i] * effectiveField[3 * i + 1] - spin[3 * i + 1] * effectiveField[3 * i] 
		- constants::lambda * 
		(-effectiveField[3 * i + 2] * spin[3 * i] * spin[3 * i]
		 - effectiveField[3 * i + 2] * spin[3 * i + 1] * spin[3 * i + 1]
		 + effectiveField[3 * i] * spin[3 * i] * spin[3 * i + 2] 
		 + effectiveField[3 * i + 1] * spin[3 * i + 1] * spin[3 * i + 2]));
	}
}

void Integrator::integrate(std::vector<std::vector<int>> &neighbours, std::vector<double> &spin, std::vector<double> &randomField)
{
	//More hardcore? Does Euler work?
	evaluate(neighbours, spin);
	for (int i = 0; i < constants::nAtoms * 3; i++)
	{
		rkPos[i] = spin[i] + rk[i] / 2.;
	}
	evaluate(neighbours, rkPos);

	double temperatureSpin[3];
	for (int i = 0; i < constants::nAtoms; i++)
	{
		temperatureSpin[0] = constants::gamma * (spin[3 * i + 1] * randomField[3 * i + 2] - spin[3 * i + 2] * randomField[3 * i + 1]);
		temperatureSpin[1] = constants::gamma * (spin[3 * i + 2] * randomField[3 * i] - spin[3 * i] * randomField[3 * i + 2]);
		temperatureSpin[2] = constants::gamma * (spin[3 * i] * randomField[3 * i + 1] - spin[3 * i + 1] * randomField[3 * i]);

		spin[3 * i]     += rk[3 * i + 0] + temperatureSpin[0];
		spin[3 * i + 1] += rk[3 * i + 1] + temperatureSpin[1];
		spin[3 * i + 2] += rk[3 * i + 2] + temperatureSpin[2];
	}
}

double Integrator::calculateEnergy(std::vector<double> &spin){
	double totalEnergy = 0;
	for (int i=0;i<3*constants::nAtoms;i++){
		totalEnergy += spin[i] * effectiveField[i]; 
	}
	return -constants::gFactor*constants::bohrMagneton*totalEnergy;
}
double Integrator::dotProduct(const double a[3], double b[3])
{
	double result = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	return result;
}

std::vector<double>* Integrator::getEffectiveField(){
	return &effectiveField;
};
void Integrator::setStabilizerField(){
	double iterAngle;
	switch (constants::spinInit) {
		case 7: 
			for (int i = 0; i < constants::nAtoms; i++){
				iterAngle = 1./2.*constants::pi-2./3.*constants::pi*(i/3)+2./3.*constants::pi*i;
				stabilizerField[3*i] = std::cos(iterAngle);
				stabilizerField[3*i+1] = std::sin(iterAngle);
				stabilizerField[3*i+2] = 0;
			}
			break;
		case 8:
			for (double i = 0; i < constants::nAtoms; i++){
			iterAngle = constants::pi*7./6.+2./3.*constants::pi*i;
			stabilizerField[3*i] = std::cos(iterAngle);
			stabilizerField[3*i+1] = std::sin(iterAngle);
			stabilizerField[3*i+2] = 0;
		}
	}
	for (int i=0;i<9;i++){
		std::cout << stabilizerField[3*i] << ", " << stabilizerField[3*i+1] << std::endl;
	}

}