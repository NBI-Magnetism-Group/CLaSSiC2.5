#include <vector>
#include <string>

namespace constants
{
    // Physical constants
    extern double boltzmann;
    extern double gFactor;
    extern double bohrMagneton;
    extern double hBar;
    extern double gamma;
    extern double pi;

    // System constants
    extern double spinSize;
    extern int nUnitCells;
    extern int nDimensions;
    extern int nAtoms;
    extern std::vector<std::vector<double>> unitVectors;
    extern std::vector<std::vector<double>> basisPosition;

    // Simulation settings;
    extern int spinInit;
    extern int geometry;
    extern double angle;
    extern int mode;
    extern int nNeighbours;
    extern double minDistance;
    extern bool periodicBoundary;
    extern bool stabilize;

    // Simulation parameters
    extern double dt;
    extern int steps;
    extern int burnInSteps;
    extern double anisotropy[3];
    extern double anisotropyStrength;
    extern double anisotropyAxis;
	extern double anisotropyPlane;
    extern double magneticField[3];
    extern double J;
    extern double lambda;
    extern double temperature;

    // File constants
    extern std::string positionFile;
    extern std::string outputFile;
    extern std::string energyFile;
    extern double offset;
    extern double length;

    // Prefactors
    extern double exchangePrefactor;
    extern double temperatureSigma;
    extern double anisotropyMatrix[9];
}

void setDerivatives();
// Enables printing for debug purposes
// bool printPositions = false;
// bool printNeighbours = false;
// bool printInitialize = false;