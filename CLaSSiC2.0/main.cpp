#include <iostream>
#include <string>
#include <cmath>
#include "simulation.hpp"
#include "constants.hpp"

//HI

int main(int argc, char *argv[])
{
	/*
	Get all the parameters from the program call
	Simulation is initialized here and run.
	*/
	for (int i = 1; i < argc; i++)
	{
		if (std::string(argv[i]) == "-dt"){
			constants::dt = std::stod(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-steps"){
			constants::steps = (int)std::stod(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-burnInSteps"){
			constants::burnInSteps = (int)std::stod(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-J"){
			constants::J = std::stod(argv[i+1])*constants::boltzmann;
		}
		else if (std::string(argv[i]) == "-lambda"){
			constants::lambda = std::stod(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-B"){
			constants::magneticField[2] = std::stod(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-anisotropyAxis"){
			constants::anisotropyAxis = std::stod(argv[i+1]);
			std::cout << "anisotropyAxis: " << constants::anisotropyAxis << std::endl;
		}
		else if (std::string(argv[i]) == "-anisotropyPlane"){
			constants::anisotropyPlane = std::stod(argv[i+1]);
			std::cout << "anisotropyPlane: " << constants::anisotropyPlane << std::endl;
		}
		else if (std::string(argv[i]) == "-T"){
			constants::temperature = std::stod(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-init"){
			constants::spinInit = std::stoi(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-angle"){
			constants::angle = std::stod(argv[i+1])*constants::pi/180.;
			std::cout << "angle: " << constants::angle << std::endl;
		}
		else if (std::string(argv[i]) == "-mode"){
			constants::mode = std::stoi(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-nCellsX"){
			constants::nUnitCells = (int)std::stod(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-periodicBoundary") {
			if (std::string(argv[i+1]) == "false") {
				constants::periodicBoundary = false;
			} else {
				constants::periodicBoundary = true;
			}
		}
		else if (std::string(argv[i]) == "-stabilize"){
			if (std::string(argv[i+1]) == "true") {
				constants::stabilize = true;
			} else {
				constants::stabilize = false;
			}
		}
		else if (std::string(argv[i]) == "-structure"){
			if (std::string(argv[i+1]) == "single"){
				constants::nDimensions = 1;
				constants::nNeighbours = 2;
				constants::nUnitCells = 1;
				constants::basisPosition = {{0., 0., 0.}};
				constants::geometry = 0;
			}
			else if (std::string(argv[i+1]) == "chain"){
				constants::nDimensions = 1;
				constants::nNeighbours = 2;
				constants::unitVectors = {{1., 0., 0.}};
				constants::basisPosition = {{0., 0., 0.}};
				constants::geometry = 1;
			}
			else if (std::string(argv[i+1]) == "square") {
				constants::nDimensions = 2;
				constants::nNeighbours = 4;
				constants::unitVectors = {{1., 0., 0.}, {0., 1., 0.}};
				constants::basisPosition = {{0., 0., 0.}};
				constants::geometry = 2;
			}
			else if (std::string(argv[i+1]) == "triangle") {
				constants::nDimensions = 2;
				constants::nNeighbours = 6;
				constants::unitVectors = {{1., 0., 0.}, {std::cos(constants::pi/3.), std::sin(constants::pi/3.), 0.}};
				constants::basisPosition = {{0., 0., 0.}};
				constants::geometry = 3;
				std::cout << "Unit vectors:\n";
				std::cout << "a0: " << constants::unitVectors[0][0] << ", " << constants::unitVectors[0][1] << ", " << constants::unitVectors[0][2] << std::endl;
				std::cout << "a1: " << constants::unitVectors[1][0] << ", " << constants::unitVectors[1][1] << ", " << constants::unitVectors[1][2] << std::endl;
			} 
			else if (std::string(argv[i+1]) == "kagome") {
				constants::nDimensions = 2;
				constants::nNeighbours = 4;
				constants::unitVectors = {{2., 0., 0.}, {2. * std::cos(constants::pi/3.), 2. * std::sin(constants::pi/3.), 0.}};
				constants::basisPosition = {{0., 0., 0.}, {1., 0., 0.}, {1. * std::cos(constants::pi/3.), 1. * std::sin(constants::pi/3.), 0.}};
				constants::geometry = 4;
			}

			else if (std::string(argv[i+1])== "hexagonal") {
				constants::nDimensions = 2;
				constants::nNeighbours = 3;
				constants::unitVectors = {{std::sqrt(3), 0., 0.}, {std::sqrt(3)*std::cos(constants::pi/3.), std::sqrt(3)*std::sin(constants::pi/3.), 0.}};
				constants::basisPosition = {{0., 0., 0.}, {std::cos(constants::pi/6.), -std::sin(constants::pi/6.), 0.}};
				constants::geometry = 5;
			}

			else if (std::string(argv[i+1])== "cubic") {
				constants::nDimensions = 3;
				constants::nNeighbours = 6;
				constants::unitVectors = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
				constants::basisPosition = {{0., 0., 0.}};
				constants::geometry = 6;
			}
			else if (std::string(argv[i+1]) == "hyperkagome") {
				bool save_time = true;
				constants::nDimensions = 3;
				constants::nNeighbours = 4;
                constants::unitVectors = {{8., 0., 0.}, {0., 8., 0.},{0, 0, 8.}};
                if (save_time) {
					constants::basisPosition = {{6., 7., 0.}, {4., 6., 1.}, {2., 5., 0.}, {3., 4., 2.}, {4., 2., 3.}, {2., 3., 4.}, {0., 2., 5.}, {1., 4., 6.}, {0., 6., 7.}, {5., 0., 2.}, {6., 1., 4.}, {7., 0., 6.}};
				}
				else{
					constants::basisPosition = {{6., 7., 0.}, {4., 6., 1.}, {2., 5., 0.}, {3., 4., 2.}, {4., 2., 3.}, {2., 3., 4.}, {0., 2., 5.}, {1., 4., 6.}, {0., 6., 7.}, {5., 0., 2.}, {6., 1., 4.}, {7., 0., 6.}, {6., 3., 0.}, {8., 2., 1.}, {7., 4., 2.}, {6., 5., 4.}, {5., 4., 6.}, {4., 6., 5.}, {2., 7., 4.}, {0., 6., 3.}, {1., 8., 2.}, {4., 2., 7.}, {2., 1., 8.}, {3., 0., 6.}};
				}
                constants::geometry = 7;
				constants::minDistance = 2.45;
			}
			else {
				std::cerr << "-structure only has the following options: single, chain, square, triangle, kagome, hexagonal, cubic, hyperkagome\n";
			}
		}

		i++;
	}
	setDerivatives();
	Simulation sim;
	sim.run();
	
	return 0;
}