#include <gfx/vec2.h>
#include "Particle.h"
#include "Force.h"
#include "SpringForce.h"
#include "Gravity.h"
#include "AngularSpringForce.h"
#include "CircularWireConstraint.h"
#include "SlidingConstraint.h"
#include "RodConstraint.h"
#include "FixedConstraint.h"
#include "Wall.h"

#include <vector>

class Scene
{
public:
    Scene();
    ~Scene();
    static void setup(std::vector<Particle *> &particles, std::vector<Force *> &forces, std::vector<Force *> &constraints, std::vector<Wall *> &walls, int scene);
};