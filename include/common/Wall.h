#pragma once

#include <gfx/vec2.h>

class Wall
{
public:
  Wall(const Vec2f p1, const Vec2f p2);

  void draw();

  Vec2f const m_p1;
  Vec2f const m_p2;
};