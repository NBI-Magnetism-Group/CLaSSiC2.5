# CLaSSiC2.5

## Installation Guide
### Github installation
Easiest way to get CLaSSiC is to git clone from this repository.
*  Install GIT from: https://git-scm.com/downloads 
*  Go to the folder you want to clone to
*  Git clone https://github.com/NBI-Magnetism-Group/CLaSSiC2.5.git

### C++ Setup
To run CLaSSiC a C++ enviroment is needed.
* Install MSYS2 from: https://www.msys2.org/
* Add the path to system variables 'C:\msys64\mingw64\bin' 

### Code editor
To run CLaSSiC a code editor is needed that both can run python and C++.
Personal preference is Visual Studio Code, which then needs the following extensions:
* C/C++
* python
* pylance
* github: pull requests and issues

### Python Enviroment
To run the plotting tools in the simulator a pyhton enviroment is needed.

Packages needed:
* pyfftw
* matplotlib
* numpy
* scipy


### Task.json structure
When everything is download and located on your computer and all the above steps are complete. Then open the 'run.bash' and press "ctrl/command"+"shift"+"b" to create the task.json file.

If this returns an error, then in the root of CLaSSiC folder, create ".vscode" and in that folder create a "tasks.json" file.

Else the "tasks.json" file needs to contain:

```
{
  "version":"2.0.0",
  "tasks": [
    {
			"type": "cppbuild",
			"label": "release",
			"command": "C:/msys64/mingw64/bin/g++.exe",
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
			"detail": "compiler: C:/msys64/mingw64/bin/g++.exe"
		},
}
```





