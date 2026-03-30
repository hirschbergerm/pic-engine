# pic-engine
A simple electrostatic particle-in-cell (ESPIC) simulation program. This is my first ESPIC and I used _Plasma Simulations by Example_ by Lubos Brieda as a guide for the overall architecture and logic. 

# Build Steps
1. Clone the project with git from the command line ```git clone https://github.com/hirschbergerm/pic-engine.git```
2. Change into root directory and create a build directory with ```mkdir build```
3. Change into the build directory ```cd build``` and run the CMake configuration step ```cmake .. -G [your build tool]```
4. Now build the exectuable by invoking your build tool ```ninja```, ```make```, or another tool

**Note:** I used Ninja for this project and so far I haven't tested any other build tools, so Ninja is highly recommended for now.

# Running and Viewing Results
1. You can run the executable by typing its path in the console ```./build/src/pic-sim``` As the program runs, it should write .vti files into the "results" folder.
2. When the program is finished, open Paraview
3. Within Paraview, select File > Open then navigate to the "results" folder in the root directory 
