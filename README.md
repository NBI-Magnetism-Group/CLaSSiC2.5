# CLaSSiC2.5

## Installation Guide
### GitHub installation
The easiest way to get CLaSSiC is to git clone from this repository.
*  Install GIT from: https://git-scm.com/downloads 
*  Go to the folder you want to clone to
*  Git clone https://github.com/NBI-Magnetism-Group/CLaSSiC2.5.git

### C++ Setup
To run CLaSSiC a C++ environment is needed.
* Install MSYS2 from: https://www.msys2.org/
* Open the MSYS2 UCRT64 environment
* Run the command 'pacman -S mingw-w64-ucrt-x86_64-gcc'
* In system variables -> environment variables -> Path: Add "C:\msys64\ucrt64\bin" 

### Code Editor
To run CLaSSiC a code editor is needed that both can run Python and C++.
My personal preference is Visual Studio Code, which then needs the following extensions:
* C/C++
* python
* pylance
* github: pull requests and issues

### Python Environment
To run the plotting tools in the simulator a python environment is needed.

Packages needed:
* pyfftw
* matplotlib
* numpy
* scipy


### Task.json structure
When everything is downloaded and located on your computer all the above steps are complete. Then open the 'run.bash' and press "ctrl/command"+"shift"+"b" to create the 'task.json' file.

If this returns an error, then in the root of CLaSSiC folder, create '.vscode'-folder and create a 'tasks.json' file in that folder.

Else the 'tasks.json' file needs to contain the following:

```
{
  "version":"2.0.0",
  "tasks": [
    {
			"type": "cppbuild",
			"label": "release",
			"command": "C:/msys64/ucrt64/bin/g++.exe",
			"args": [
				"-fdiagnostics-color=always",
				"-O3",
				"'Your location to CLaSSiC'/CLaSSiC2.0/CLaSSiC2.0/*.cpp",
				"-o",
				"'Your location to CLaSSiC'/CLaSSiC2.0/CLaSSiC2.0/model.out",
				"-std=c++17",
				"-mavx2",
				"-march=native"
			],
			"problemMatcher": ["$gcc"],
			"group": {
				"kind": "build",
				"isDefault": true
			},
			"detail": "compiler: C:/msys64/ucrt64/bin/g++.exe"
		},
}
```

## Guide to run CLaSSiC
With the 'tasks.json' set-up and the environment running, then open the 'run.bash' and press "ctrl/command"+"shift"+"b" to run the 'task.json' file and this should make the 'model.out' file.
With this file created, CLaSSiC is ready to run.


To run CLaSSiC, you specify the structure and variables in the 'run.bash' file and then to compute the simulation you "cd 'Location of run.bash'" in the terminal and run the command "bash run.bash".

### Variable documentation
* Time Steps: **dt**
	* Step width [sec]
* Steps: **steps**
	* Number of steps [#]
* Exchange Constant: **J**
	* [Joule]
* Relaxation Constant: **Lambda**
	* [Dimension less]
* Magnetic field: **B**
	* External magnetic field, [Tesla]
* Axial Anisotropy Strength: **AnisotropyAxis**
	* Axial anisotropy field strength [Tesla]
* Planar Anisotropy Strength: **AnisotropyPlane**
	* Planar anisotropy field strength [Tesla]
* Temperature: **T**
	* Temperature of sample/system [Kelvin]
* Initial Spin: **Init**
	* Case determination for spin calculation 0,1,…,7,8 [#]
		- Case 0: Aligned with z-axis
		- Case 1: Angle with z-axis
		- Case 2: Small z angle for spin waves
		- Case 3: Rotor mode for two spins
		- Case 4: 2D spin waves
		- Case 5: Triangle
		- Case 6: Kagome afm sqrt(3) x sqrt(3)
		- Case 7: Kagome Q0
		- Case 8: Uses spin orientations in the Saved_spin.dat file to initialize 
* Angle: **angle**
	* Angle between spins [degree]
	* Only affects the system in Case 1, 3, 6 & 7 
* Manual zero mode: **Mode**
	* Only used for case 2 [#]
* Structure: **structure**
	* Possible structures: [text]
		* single
		* chain
		* square
		* triangle
		* cubic
		* hexagonal
		* kagome
		* hyperkagome	
* Total unit cells: **nCellsX**
	* Number of unit cells calculated [#]
		* Need to correlate with the structure for AntiFerroMagnetic spin	
* Periodic Boundary: **PeriodicBoundary**
	* [True/False]
* Field stabilization: **Stabalize**
	* Only effects if the structure is kagome [True/False]
		* Simplifies spin calculation to only take iteration angle and not specified angle into account

