#include <vector>
#include "integrator.hpp"

class Simulation
{
private:
	std::vector<double> spin = std::vector<double>(constants::nAtoms * 3);
	std::vector<double> position;
	std::vector<double> randomField = std::vector<double>(constants::nAtoms * 3, 0);
	std::vector<std::vector<int>> neighbours = std::vector<std::vector<int>>(constants::nAtoms, std::vector<int>());
	std::vector<double> totalEnergy = std::vector<double>((int)constants::steps/100);
	Integrator integrator;

public:
	Simulation();
	void load();
	void initialize();
	void run();
	void normalize();
	std::string addFileNumber(std::string input);
	void writeConstants(std::ofstream &file);
	void writePositions();
	double distance(std::vector<double>::iterator a, std::vector<double>::iterator b);
	void addNeighbours(std::vector<double> a, std::vector<double> b, int i, int j);

	std::vector<double>* getSpin();
	std::vector<double>* getPosition();
	std::vector<double>* getRandomField();
	std::vector<std::vector<int>>* getNeighbours();
};