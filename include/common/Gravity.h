#pragma once
#include "Particle.h"
#include "Force.h"

class Gravity : Force
{
public:
    Gravity(Particle *particle);
    ~Gravity();
    void draw(bool draw[]) override;
    void calculate() override;
    Particle *p;
};