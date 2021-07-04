#pragma once

#include <gfx/vec2.h>
#include <gfx/mat2.h>
#include <gfx/vec4.h>
#include <gfx/mat4.h>
#include <vector>
#include <Particle.h>

class Force
{
public:
    Force();
    virtual ~Force(void);

    virtual void draw(bool draw[]);
    virtual void calculate();
    Vec2f force;
};