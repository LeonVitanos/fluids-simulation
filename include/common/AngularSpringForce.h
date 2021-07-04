#pragma once

#include "Particle.h"
#include "Force.h"

class AngularSpringForce : Force
{
public:
    AngularSpringForce(Particle *p1, Particle *p2, Particle *p3, double dist, double k_s, double k_d);

    void draw(bool draw[]) override;
    void calculate() override;

private:
    Particle *const m_p1;
    Particle *const m_p2;
    Particle *const m_p3;
    double const m_dist;
    double const m_ks, m_kd;
};