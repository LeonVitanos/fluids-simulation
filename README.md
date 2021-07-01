# SiCG-Part-1

Simulation in Computer Graphics Project Part 1

## Running

### Makefile (Mac, VS Code)

- Make sure to `brew install libpng` and install C/C++ and MakeFile Tools packages
- Run using the makefile extension

### CMake (Linux, VS Code)

- Make sure to install CMake, C/C++, libpng and freeglut
- Then run `cmake -B build`
- Next `cmake --build build`
- And lastly, to run, `./build/2imv15_project_1`

## Next steps:

- [ ] Add boundaries to the objects
  - [ ] Fix crashing when moving object to top
  - [ ] Moving object e.g to the right, makes it go to the left (max dis=distance to center of object)
- [ ] (maybe) have a second look at the fluid movement of the objects
- [ ] Rotations
- [ ] Two-way coupling
- [ ] Particles and boundaries
