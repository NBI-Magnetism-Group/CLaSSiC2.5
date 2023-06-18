#include "simulation.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <filesystem>
#include <chrono>
#include <vector>
#include <sstream>
#include <string>

Simulation::Simulation()
{
	/*
	Constructer for the simulation class.
		- make sure the data directory exists
		- Calls load to build the crystal
		- Calls initialize to set to start orientation of the spins
	*/
	if (!std::filesystem::is_directory("data") || !std::filesystem::exists("data")) {
    	std::filesystem::create_directory("data"); // create src folder
	}

	load();
	initialize();
}

void Simulation::load()
{
	/*
	Builds the crystal from the ground up.
		1. Makes an unit cell
		2. Extends the unit cell in specified dimensions
		3. Determine nearest neighbours
		4. Apply periodic boundary conditions
	*/

	// Create the basis
	for (int i=0; i < constants::basisPosition.size(); i++){
		for (int j = 0; j < 3; j++){
			position.push_back(constants::basisPosition[i][j]);
		}
	}


	// Add more unit cells
	int offset = 0;
	for (int dim = 0; dim < constants::nDimensions; dim++)
	{
		for (int cell = 1; cell < constants::nUnitCells; cell++)
		{
			if (constants::nDimensions>1){
				if (position.size()%(3*constants::basisPosition.size()*2*constants::nUnitCells)==0 && constants::unitVectors[1][0]!=0){
					offset++;
					// std::cout << "Position size: " << position.size()/3 << std::endl;
				}
			}
			// std::cout <<  std::pow(constants::nUnitCells, dim) << std::endl;
			for (int i = 0; i < std::pow(constants::nUnitCells, dim); i++)
			{
				// std::cout <<  constants::basisPosition.size() << std::endl;
				for (int k=0; k<constants::basisPosition.size(); k++){
					for (int j = 0; j < 3; j++)
					{
						// std::cout << "vec: "<<  constants::unitVectors[dim][j] << std::endl;
						// position.push_back(position[3 * i * dim + j + 3*(k*(1+i))] + (constants::unitVectors[dim][j]) * cell-(constants::unitVectors[0][j])*offset);
						// std::cout << position.size()/3 << ": "<<position[3 * k + 3 * i + j] << " | " <<  (constants::unitVectors[dim][j]) * cell << " | " << -(constants::unitVectors[1][j])*offset << std::endl;
						if (j==0){
							position.push_back(position[3 * k + constants::basisPosition.size() * 3 * i + j]+constants::unitVectors[dim][j] * cell-(constants::unitVectors[0][j])*offset);
						} else {
							position.push_back(position[3 * k + constants::basisPosition.size() * 3 * i + j]+constants::unitVectors[dim][j] * cell);
						}
					}
				}
			}
		}
	}
	writePositions();
	// std::cout << "positions complete!\n";

	// Calculate nearest neighbours
	for (int i = 0; i < constants::nAtoms; i++)
	{
		for (int j = 0; j < constants::nAtoms; j++)
		{
			if (i != j)
			{
				// if (i==4){
				// 	std::cout << "distance: " << distance(position.begin()+3*i,position.begin()+3*j) << std::endl;
				// }
				if (distance(position.begin()+3*i,position.begin()+3*j) < constants::minDistance)
				{
					neighbours[i].push_back(j);
				}
			}
		}
	}

	// std::cout << "nearest neighbours complete!\n";
	// Periodic boundary conditions
	bool applyBoundry = true;
	if (constants::periodicBoundary){
		for (int i = 0; i < constants::nAtoms; i++)
		{ 
			// if (neighbours[i].size()!=constants::nNeighbours)
			if (true)
			{
				std::vector<double> xBoundary = std::vector<double>(3);
				std::vector<double> yBoundary = std::vector<double>(3);
				std::vector<double> zBoundary = std::vector<double>(3);
				std::vector<double> xyBoundary = std::vector<double>(3);
				std::vector<double> xyOrthoBoundary = std::vector<double>(3);
				std::vector<double> xzBoundary = std::vector<double>(3);
				std::vector<double> xzOrthoBoundary = std::vector<double>(3);
				std::vector<double> yzBoundary = std::vector<double>(3);
				std::vector<double> yzOrthoBoundary = std::vector<double>(3);

				xBoundary[0] = position[3*i] + constants::nUnitCells * constants::unitVectors[0][0];
				xBoundary[1] = position[3*i+1];
				xBoundary[2] = position[3*i+2];

				if (constants::nDimensions > 1){
					yBoundary[0] = position[3*i];
					yBoundary[1] = position[3*i+1] + constants::nUnitCells * constants::unitVectors[1][1];
					yBoundary[2] = position[3*i+2];

					xyBoundary[0] = position[3*i] + constants::nUnitCells * constants::unitVectors[0][0];
					xyBoundary[1] = position[3*i+1] + constants::nUnitCells * constants::unitVectors[1][1];
					xyBoundary[2] = position[3*i+2];

					xyOrthoBoundary[0] = position[3*i] - constants::nUnitCells * constants::unitVectors[0][0];
					xyOrthoBoundary[1] = position[3*i+1] + constants::nUnitCells * constants::unitVectors[1][1];
					xyOrthoBoundary[2] = position[3*i+2];	
				}

				if (constants::nDimensions > 2) {
					zBoundary[0] = position[3*i];
					zBoundary[1] = position[3*i+1];
					zBoundary[2] = position[3*i+2] + constants::nUnitCells * constants::unitVectors[2][2];

					xzBoundary[0] = position[3*i] + constants::nUnitCells * constants::unitVectors[0][0];
					xzBoundary[1] = position[3*i+1];
					xzBoundary[2] = position[3*i+2] + constants::nUnitCells * constants::unitVectors[2][2];

					xzOrthoBoundary[0] = position[3*i] - constants::nUnitCells * constants::unitVectors[0][0];
					xzOrthoBoundary[1] = position[3*i+1];
					xzOrthoBoundary[2] = position[3*i+2] + constants::nUnitCells * constants::unitVectors[2][2];

					yzBoundary[0] = position[3*i];
					yzBoundary[1] = position[3*i+1] + constants::nUnitCells * constants::unitVectors[1][1];
					yzBoundary[2] = position[3*i+2] + constants::nUnitCells * constants::unitVectors[2][2];

					yzOrthoBoundary[0] = position[3*i];
					yzOrthoBoundary[1] = position[3*i+1] - constants::nUnitCells * constants::unitVectors[1][1];
					yzOrthoBoundary[2] = position[3*i+2] + constants::nUnitCells * constants::unitVectors[2][2];
				}


				for (int j = 0; j < constants::nAtoms; j++)
				{
					if (i!=j){
						// x-boundary conditions
						addNeighbours(xBoundary, position, i, j);
						// y-boundary conditions
						if (constants::nDimensions > 1){
							addNeighbours(yBoundary, position, i, j);
							// Connecting over the diagonal
							addNeighbours(xyBoundary, position, i, j);
							addNeighbours(xyOrthoBoundary, position, i, j);
						}
						// z-boundary conditions
						if (constants::nDimensions > 2) {
							addNeighbours(zBoundary, position, i, j);
							addNeighbours(xzBoundary, position, i ,j);
							addNeighbours(xzOrthoBoundary, position, i, j);
							addNeighbours(yzBoundary, position, i ,j);
							addNeighbours(yzOrthoBoundary, position, i ,j);
						}
					}
				}
			}	
		}	
	}
	// std::cout << "periodic boundaries complete!\n";

	// Printing

	// std::cout << "Atom positions: \n";
	// for (int i = 0; i < position.size() / 3; i++)
	// {
	// 	std::cout << "{" << position[3 * i] << " , " << position[3 * i + 1];
	// 	if (constants::nDimensions > 2) {
	// 		std::cout <<" , " << position[3 * i + 2];
	// 	}
	// 	std::cout <<"}, ";
	// }


	// std::cout << "Neighbours: \n";
	// for (int i = 0; i < constants::nAtoms; i++)
	// {
	// 	if (neighbours[i].size() != constants::nNeighbours){
	// 		std::cout << i << "(" << neighbours[i].size() << ")"<<  ": ";
	// 		for (int j = 0; j < neighbours[i].size(); j++)
	// 		{
	// 			std::cout << neighbours[i][j] << " ";
	// 		}
	// 		std::cout << std::endl;
	// 	}
	// }

}

void Simulation::initialize()
{
	/*
	Initializing the spin orientations
	*/

	int swap = 1;
	double delta = 1;
	double iterAngle;
	switch (constants::spinInit)
	{
	case 0:
		// Alligned with the z-axis
		for (int i = 0; i < constants::nAtoms; i++)
		{
			spin[3 * i] = 0;
			spin[3 * i + 1] = 0;
			spin[3 * i + 2] = 1;
		}
		break;
	case 1:
		// Angle with z-axis
		for (int i = 0; i < constants::nAtoms; i++)
		{
			spin[3 * i] = std::cos((float)i / constants::nAtoms * 2 * constants::pi) * std::sin(constants::angle);
			spin[3 * i + 1] = std::sin((float)i / constants::nAtoms * 2 * constants::pi) * std::sin(constants::angle);
			spin[3 * i + 2] = std::cos(constants::angle);

			// spin[6 * i + 3] = std::cos((float)i / constants::nAtoms * 2 * constants::pi) * std::sin(constants::angle);
			// spin[6 * i + 4] = std::sin((float)i / constants::nAtoms * 2 * constants::pi) * std::sin(constants::angle);
			// spin[6 * i + 5] = -std::cos(constants::angle);
		}
		break;
	case 2:
		// Small z angle for spin waves
		for (int i = 0; i < (int)constants::nAtoms; i++)
		{
				spin[3 * i] = 0.001 * std::cos((float)i / constants::nAtoms * constants::mode * 2 * constants::pi);
				spin[3 * i + 1] = 0.001 * std::sin((float)i / constants::nAtoms * constants::mode * 2 * constants::pi);
				spin[3 * i + 2] = 1;
			if (constants::J < 0)
				{
					// spin[3 * i + 3] = 0.001 * std::cos((float)i / constants::nAtoms * constants::mode * 2 * constants::pi);
					// spin[3 * i + 4] = -0.001 * std::sin((float)i / constants::nAtoms * constants::mode * 2 * constants::pi);
					// spin[3 * i + 5] = -1;
					i++;
					spin[3 * i ] = -0.001 * std::cos((float)i / constants::nAtoms * constants::mode * 2 * constants::pi);
					spin[3 * i + 1] = -0.001 * std::sin((float)i / constants::nAtoms * constants::mode * 2 * constants::pi);
					spin[3 * i + 2] = -1;
				}
		}
		normalize();
		break;
	case 3:
		if (constants::J < 0)
		{
			swap = -1;
		}
		// Rotor mode for two spins
		spin[0] = std::cos(constants::angle);
		spin[1] = std::sin(constants::angle);
		spin[2] = std::cos(constants::angle);
		spin[3] = -std::cos(constants::angle);
		spin[4] = -std::sin(constants::angle); 
		spin[5] = std::cos(constants::angle);
		break;
	case 4:
		// 2D spin waves
		for (int i = 0; i < constants::nAtoms; i++)
		{
			spin[3 * i] = 0;
			spin[3 * i + 1] = 0;
			if ((i + i / constants::nUnitCells) % 2 == 0)
			{
				spin[3 * i + 2] = 1;
			}
			else
			{
				spin[3 * i + 2] = -1;
			}
		}
	break;
	case 5:
		// Triangle
		std::cout << constants::spinInit;	
		std::cout << "Triangle:\n";
		spin[0] = delta*std::cos(30./180.*constants::pi);
		spin[1] = delta*std::sin(30./180.*constants::pi);
		spin[2] = 1;

		spin[3] = -delta*std::cos(30./180.*constants::pi);
		spin[4] = delta*std::sin(30./180.*constants::pi);
		spin[5] = 1;
		
		spin[6] = delta*0;
		spin[7] = -delta*1;
		spin[8] = -1;
	break;
		case 6:
		// Kagome antiferromagnetic sqrt(3) x sqrt(3)
		double offset;
		for (int i = 0; i < constants::nAtoms; i++){
			offset = 2./3.*constants::pi*(i/3/constants::nUnitCells%2);
			iterAngle = 1./2.*constants::pi-2./3.*constants::pi*(i/3)+2./3.*constants::pi*i + offset;
			spin[3*i] = std::cos(iterAngle);
			spin[3*i+1] = std::sin(iterAngle);
			// Pertubation
			if (std::abs(spin[3*i]) > 1e-6){
				if (spin[3*i] > 1e-6){
					spin[3*i] += std::cos(-constants::angle+iterAngle)-std::cos(iterAngle);
					spin[3*i+1] += std::sin(-constants::angle+iterAngle)-std::sin(iterAngle);
				}
				else {
					spin[3*i] += std::cos(constants::angle+iterAngle)-std::cos(iterAngle);
					spin[3*i+1] += std::sin(constants::angle+iterAngle)-std::sin(iterAngle);
				}
			spin[3*i+2] = 0;
			}
		}
	break;
		case 7:
		// Kagome Q0
		for (double i = 0; i < constants::nAtoms; i++){
			iterAngle = constants::pi*7./6.+2./3.*constants::pi*i;
			spin[3*i] = std::cos(iterAngle);
			spin[3*i+1] = std::sin(iterAngle);

			// Pertubation
			if (std::abs(spin[3*i]) > 1e-6){
				if (spin[3*i] > 1e-6){
					spin[3*i] += std::cos(-constants::angle+iterAngle)-std::cos(iterAngle);
					spin[3*i+1] += std::sin(-constants::angle+iterAngle)-std::sin(iterAngle);
				}
				else {
					spin[3*i] += std::cos(constants::angle+iterAngle)-std::cos(iterAngle);
					spin[3*i+1] += std::sin(constants::angle+iterAngle)-std::sin(iterAngle);
				}
			}
			spin[3*i+2] = 0;
		}
	break;
		case 8:
		// Uses spin orientations in the Saved_spin.dat file to initialize 
		// Load file
		std::string file = "Saved_spin.dat";
		std::ifstream saved_spins;
		saved_spins.open(file);

		// Get the data from the file
		int n = constants::nAtoms*3;

    	int position = 0;

		// Pulling data out into loaded_spin
		while(! saved_spins.eof() && position < n)
		{
    		saved_spins >> spin[position];
			position++;
		}
		std::cout << "Initial state set" << std::endl;
	break;
	}
	// Printing

	// std::cout << "Spin initialization: \n";
	// for (int i = 0; i < constants::nAtoms; i++)
	// {
	// 	std::cout << "Atom: " << i / constants::nUnitCells << ", " << i % constants::nUnitCells << "|x: " << spin[3 * i] << " y: " << spin[3 * i + 1] << " z: " << spin[3 * i + 2] << std::endl;
	// }
}

void Simulation::run()
{
	/*
	Actually runs the simulation
		1. Write constants to file
		2. Determine temperature fluctuations
		3. Stopwatch to track time
		4. Actual run
		5. close file and print time
	*/

	// Open and write to file
	std::ofstream fileSpin;
	std::ofstream fileEnergy;
	fileSpin.open(addFileNumber(constants::outputFile), std::ios::binary);
	writeConstants(fileSpin);

	// Temperature fluctatuations
	std::random_device seed;
	std::default_random_engine engine{seed()};
	std::normal_distribution<double> temperatureDistribution;
	if (constants::temperatureSigma != 0 && constants::temperature != 0)
	{
		temperatureDistribution = std::normal_distribution<double>(0, constants::temperatureSigma);
	}

	// Burn in with a higher temperature:
	std::normal_distribution<double> burnInDistribution;
	if (constants::burnInSteps>0){
		burnInDistribution = std::normal_distribution<double>(0, 10.*constants::temperatureSigma);
		
		for (int i = 0; i < constants::burnInSteps; i++)
		{
			if (constants::temperatureSigma != 0 && constants::temperature != 0)
			{
				for (int i = 0; i < constants::nAtoms * 3; i++)
				{
					randomField[i] = burnInDistribution(engine);
				}
				integrator.integrate(neighbours, spin, randomField);
				normalize();
			}
		}
	}


	//Set the stabilizerField
	if (constants::geometry==4 && constants::stabilize){
		integrator.setStabilizerField();
	}

	// Stopwatch
	auto start = std::chrono::high_resolution_clock::now();

	// Actual run
	for (int i = 0; i < constants::steps; i++)
	{
		if (constants::temperatureSigma != 0 && constants::temperature != 0)
		{
			for (int i = 0; i < constants::nAtoms * 3; i++)
			{
				randomField[i] = temperatureDistribution(engine);
			}
		}


		if (i % 100 == 0)
		{
			std::cout << "Progress: " << (float)i / constants::steps * 100 << "%\r";
		}
		if (i % 100 == 0)
		{
			fileSpin.write((char *)&spin[0], sizeof(double) * constants::nAtoms * 3);
			totalEnergy[i/100]=integrator.calculateEnergy(spin);
		}
		integrator.integrate(neighbours, spin, randomField);
		normalize();
	}

	// Wrap up
	fileSpin.close();
	fileEnergy.open(addFileNumber(constants::energyFile), std::ios::binary);
	fileEnergy.write((char *) &totalEnergy[0], sizeof(double)*((int)constants::steps/100));
	fileEnergy.close();
	auto msInt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
	std::cout << "\nDuration: " << (float)msInt.count() / 1000 << "seconds\n"
			  << std::endl;
}

void Simulation::normalize()
{
	/*
	Helper function to normalize all of the spins
	*/
	double invLength;
	for (int i = 0; i < constants::nAtoms; i++)
	{
		invLength = 1 / std::sqrt(spin[3 * i] * spin[3 * i] + spin[3 * i + 1] * spin[3 * i + 1] + spin[3 * i + 2] * spin[3 * i + 2]);
		spin[3 * i] *= invLength;
		spin[3 * i + 1] *= invLength;
		spin[3 * i + 2] *= invLength;
	}
}

std::string Simulation::addFileNumber(std::string input){
	// Manipulates input to add the next number for the file that will be generated
	int i=0;

	std::string fileString = input;
	fileString.insert(input.find("."), std::to_string(i));
	while (std::filesystem::exists(fileString))
	{
		i++;
		fileString = input;
		fileString.insert(input.find("."), std::to_string(i));
	}
	return fileString;
}

void Simulation::writeConstants(std::ofstream &f)
{
	/*
	Writes all constants and system parameters to the start of the file
	*/
	double saveAtoms = constants::nAtoms;
	double saveSteps = constants::steps;
	double saveGeometry = constants::geometry;
	double saveNUnitcells = constants::nUnitCells;
	f.write((char *)&constants::offset, sizeof(double));
	f.write((char *)&saveAtoms, sizeof(double));
	f.write((char *)&constants::dt, sizeof(double));
	f.write((char *)&saveSteps, sizeof(double));
	f.write((char *)&constants::J, sizeof(double));
	f.write((char *)&constants::lambda, sizeof(double));
	f.write((char *)&constants::magneticField[0], 3 * sizeof(double));
	f.write((char *)&constants::anisotropy[0], 3 * sizeof(double));
	f.write((char *)&constants::temperature, sizeof(double));
	f.write((char *)&constants::length, sizeof(double));
	f.write((char *)&saveGeometry, sizeof(double));
	f.write((char *)&saveNUnitcells, sizeof(double));
	f.write((char *)&constants::anisotropyStrength, sizeof(double));
}

void Simulation::writePositions() {
	// Writes the physical positions of the atoms to file
	std::ofstream positionFile;
	positionFile.open(addFileNumber(constants::positionFile));

	std::string seperator;
	for (int i=0; i < position.size();i++){
		positionFile << position[i];
		seperator = ((i+1)%3!=0) ? ", " : "\n";
		positionFile << seperator;
	}
	positionFile << "\n";
	positionFile.close();
}

double Simulation::distance(std::vector<double>::iterator a, std::vector<double>::iterator b){
	//Returns the distance between a and b. Takes the first three elements as Carthesian coordinates
	double distance;
	distance = std::sqrt((*a-*b)*(*a-*b)
						+(*(a+1)-*(b+1))*(*(a+1)-*(b+1))
						+(*(a+2)-*(b+2))*(*(a+2)-*(b+2)));
	return distance;
}

void Simulation::addNeighbours(std::vector<double> a, std::vector<double> b, int i, int j){
	//Checks if the point i and j are neighbours and adds them if true
		if (distance(a.begin(), b.begin()+3*j) < constants::minDistance) {
			if (std::find(neighbours[i].begin(), neighbours[i].end(), j) == neighbours[i].end()) {
				neighbours[i].push_back(j);
				neighbours[j].push_back(i);
		}
	}
}

std::vector<double>* Simulation::getSpin(){
	return &spin;
}

std::vector<double>* Simulation::getPosition(){
	return &position;
}

std::vector<double>* Simulation::getRandomField(){
	return &randomField;
}

std::vector<std::vector<int>>* Simulation::getNeighbours(){
	return &neighbours;
}