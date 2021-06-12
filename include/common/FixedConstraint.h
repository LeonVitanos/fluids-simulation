#pragma once

#include "Particle.h"
#include "Force.h"

class FixedConstraint : Force
{
public:
    FixedConstraint(Particle *p);

    void draw(bool draw[]) override;
    void calculate() override;
    Particle *p;
};