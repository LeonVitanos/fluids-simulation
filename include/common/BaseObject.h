#include <vector>
#include <math.h>
#include "Boundary.h"

class BaseObject
{
public:
    BaseObject();

    virtual void draw() = 0;
    virtual void update(BoundaryCell *boundaries, float dt) = 0;
    virtual void setPosition(float x, float y) = 0;
    virtual void setVelocity(float u, float v) = 0;
    virtual bool isOnCell(float x, float y) = 0;
    virtual std::vector<float> getPosition() = 0;
    virtual std::vector<float> getVelocity() = 0;
};