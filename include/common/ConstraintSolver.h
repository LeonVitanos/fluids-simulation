#pragma once

#include "Force.h"
#include "Particle.h"

#include <vector>

class ConstraintSolver
{

private:
  std::vector<Particle *> particles;
  std::vector<Force *> constraints;

public:
    ConstraintSolver(std::vector<Particle *> particles, std::vector<Force *> constraints);

    void calculate();
};