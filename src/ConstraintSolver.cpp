#include "ConstraintSolver.h"

#include <vector>
#include "Force.h"

ConstraintSolver::ConstraintSolver(std::vector<Particle *> particles, std::vector<Force *> constraints): particles(particles), constraints(constraints) {}

void ConstraintSolver::calculate()
{
	for (auto &constraint : constraints)
	{
		constraint->calculate();
	}

    // First get all the J, J_dot, C, C_dot
}