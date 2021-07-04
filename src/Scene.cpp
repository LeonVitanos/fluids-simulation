#include "Scene.h"
#define PI 3.1415926535897932384626433832795

void Scene::setup(std::vector<Particle *> &particles, std::vector<Force *> &forces, std::vector<Force *> &constraints, std::vector<Wall *> &walls, int scene)
{

    switch (scene)
    {
    case 1: // Simple Spring
    {        
        const double dist = 0.6;
        const Vec2f center(0.0, 0.0);
        const Vec2f offset(0.9, 0.0);

        particles.push_back(new Particle(center - offset));
        particles.push_back(new Particle(center + offset));

        forces.push_back((Force *)new SpringForce(particles[0], particles[1], dist, 1, 1));
        break;
    }
    case 2: // Spring with FixedConstraint (Pendulum)
    {
        const double dist = 0.3;
        const Vec2f center(0.0, 0.0);
        const Vec2f offset(dist, 0.0);

        particles.push_back(new Particle(center - 2 * offset));
        particles.push_back(new Particle(center));

        forces.push_back((Force *)new Gravity(particles[0]));
        forces.push_back((Force *)new SpringForce(particles[0], particles[1], dist, 1, 1));

        constraints.push_back((Force *)new FixedConstraint(particles[1]));
        break;
    }
    case 3:
    {
        const double dist = 0.2;
        const Vec2f center(0.0, 0.0);
        const Vec2f offset(dist, 0.0);

        // Create three particles, attach them to each other, then add a
        // circular wire constraint to the first.

        particles.push_back(new Particle(center + offset));
        particles.push_back(new Particle(center + offset + offset + offset));
        particles.push_back(new Particle(center + offset + offset + offset + offset));

        forces.push_back((Force *)new Gravity(particles[0]));
        forces.push_back((Force *)new Gravity(particles[1]));
        forces.push_back((Force *)new Gravity(particles[2]));

        // Set the spring force
        forces.push_back((Force *)new SpringForce(particles[0], particles[1], dist, 1, 1));

        constraints.push_back((Force *)new CircularWireConstraint(particles[0], center, dist));
        constraints.push_back((Force *)new RodConstraint(particles[1], particles[2], dist * 2));
        break;
    }
    case 4:
    { // Cloth
        const double dist = 0.05;
        Vec2f center(0.0, 0.0);
        const Vec2f offset(dist, 0.0);

        for (int i = 0; i < 11; i++)
        {
            center = Vec2f(0.0, 0.0 + 5 * dist - dist * i);
            for (int j = -4; j <= 4; j++)
            {
                particles.push_back(new Particle(center + j * offset));
                forces.push_back((Force *)new Gravity(particles.back()));
                if (i == 0)
                {
                    //Apply fixed constraint on the first row particles
                    constraints.push_back((Force *)new SlidingConstraint(particles.back(), 0.25));
                    if (j != -4)
                    {
                    }
                    // Rod constraints here
                }
                if (j != -4)
                    //spring connecting horizontally
                    forces.push_back((Force *)new SpringForce(particles[i * 9 + j + 3], particles[i * 9 + j + 4], dist, 0.7, 0.8));
                if (i != 0)
                {
                    //spring connecting vertically
                    forces.push_back((Force *)new SpringForce(particles[i * 9 - 9 + j + 4], particles[i * 9 + j + 4], dist, 0.7, 0.8));
                    if (j != 4)
                        //spring connecting diagonally to the right
                        forces.push_back((Force *)new SpringForce(particles[i * 9 + j + 4], particles[i * 9 + j + 4 - 8], dist, 0.05, 0.8));
                    if (j != -4)
                        //spring connecting diagonally to the left
                        forces.push_back((Force *)new SpringForce(particles[i * 9 + j + 4], particles[i * 9 + j + 4 - 10], dist, 0.05, 0.8));
                }
            }
        }

        Vec2f p1(-0.5, 0.6);
        Vec2f p2(-0.5, -0.6);
        walls.push_back(new Wall(p1, p2));

        Vec2f p3(0.5, 0.6);
        Vec2f p4(0.5, -0.6);
        walls.push_back(new Wall(p3, p4));
        break;
    }
    case 5:
    {
        // Angular spring

        // Set up particle locations
        const Vec2f p1(0.0, 0.0);
        const Vec2f p2(0.1, 0.0);
        const Vec2f p3(0.1, 0.1);
        const Vec2f p4(0.2, 0.1);
        const Vec2f p5(0.3, 0.1);

        // Create particles
        particles.push_back(new Particle(p1));
        particles.push_back(new Particle(p2));
        particles.push_back(new Particle(p3));
        particles.push_back(new Particle(p4));
        particles.push_back(new Particle(p5));

        // Apply gravity
        for (int i = 0; i < particles.size(); i++)
        {
            forces.push_back((Force *)new Gravity(particles[i]));
        }

        // Set up Angular Springs
        for (int i = 0; i < particles.size() - 2; i++)
        {
            forces.push_back((Force *)new AngularSpringForce(particles[i], particles[i + 1], particles[i + 2], 0.1, 0.005, 1));
        }

        // Set up rod constraints to prevent lenghtening of hair
        for (int i = 0; i < particles.size() - 1; i++)
        {
            constraints.push_back((Force *)new RodConstraint(particles[i], particles[i + 1], 0.1));
        }

        // Fix the first particle
        constraints.push_back((Force *)new FixedConstraint(particles[0]));
    }
    }
};
