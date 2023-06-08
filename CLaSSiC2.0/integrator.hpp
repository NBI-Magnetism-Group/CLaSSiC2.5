#include "constants.hpp"
#include <vector>

class Integrator
{
private:
	std::vector<double> rk = std::vector<double>(3 * constants::nAtoms);
	std::vector<double> rkTot = std::vector<double>(3 * constants::nAtoms);
	std::vector<double> rkPos = std::vector<double>(3 * constants::nAtoms);
	std::vector<double> effectiveField = std::vector<double>(3 * constants::nAtoms);
	std::vector<double> stabilizerField = std::vector<double>(3 * constants::nAtoms);

public:
	void calculateEffectiveField(std::vector<std::vector<int>> &neighbours, std::vector<double> &spin);
	void evaluate(std::vector<std::vector<int>> &neighbours, std::vector<double> &spin);
	void integrate(std::vector<std::vector<int>> &neighbours, std::vector<double> &spin, std::vector<double> &randomField);
	double calculateEnergy(std::vector<double> &spin);
	double dotProduct(const double a[3], double b[3]);

	std::vector<double>* getEffectiveField();
	void setStabilizerField();
};