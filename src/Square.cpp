#include "Square.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define IX(i, j) ((i) + (grid_size + 2) * (j))

Square::Square(float x, float y, int height, int width, int N) : o_x(x), o_y(y), o_h(height), o_w(width), grid_size(N)
{
}

void Square::draw()
{
    glLineWidth(1.0f);
    glBegin(GL_LINES);

    float h = 1.0f / grid_size;
    float x1 = (o_x - o_w / 2) * h;
    float x2 = (o_x + o_w / 2 ) * h;
    float y1 = (o_y - o_h / 2 ) * h;
    float y2 = (o_y + o_h / 2) * h;
    glColor3f(0.6, 0.2, 0.4);
    glVertex2f(x1, y1);
    glVertex2f(x2, y1);
    glVertex2f(x1, y1);
    glVertex2f(x1, y2);
    glVertex2f(x2, y2);
    glVertex2f(x1, y2);
    glVertex2f(x2, y1);
    glVertex2f(x2, y2);
    glEnd();
}

void Square::update(BoundaryCell *boundaries, float dt)
{
    float h = 1.0f / grid_size;
    float x1 = (o_x - o_w / 2) * h;
    float x2 = (o_x + o_w / 2) * h;
    float y1 = (o_y - o_h / 2) * h;
    float y2 = (o_y + o_h / 2) * h;

    std::vector<float> pos;
    pos.push_back(o_x - o_w / 2);
    pos.push_back(o_y + o_h / 2);
    pos.push_back(o_x + o_w / 2);
    pos.push_back(o_y - o_h / 2);

    for (int j = pos[3]; j < pos[1]; j++)
    {
        boundaries[IX((int)pos[0], j)].b_left = true;
        boundaries[IX((int)pos[0],j)].b_x = pos[0];
        boundaries[IX((int)pos[0],j)].b_y = j;
        boundaries[IX((int)pos[2]-1, j)].b_right = true;
        boundaries[IX((int)pos[2],j)].b_x = pos[2];
        boundaries[IX((int)pos[2],j)].b_y = j;
    }

    for (int i = pos[0]; i < pos[2]; i++) {
        boundaries[IX(i, (int)pos[1])].b_bottom = true;
        boundaries[IX(i, (int)pos[1])].b_x = i;
        boundaries[IX(i, (int)pos[1])].b_y = pos[1];
        boundaries[IX(i, (int)pos[3] - 1)].b_top = true;
        boundaries[IX(i, (int)pos[3])].b_x = i;
        boundaries[IX(i, (int)pos[3])].b_y = pos[3];
    }

    // Put square back in boundaries as its not allowed to leave
    if (x2 >= 1)
        o_x = grid_size - o_w / 2 - 0.5f;
    if (x1 <= 0)
        o_x = o_w / 2 + 0.5f;
    if (y2 >= 1)
        o_y = grid_size - o_h / 2 - 0.5f;
    if (y1 <= 0)
        o_y = o_h / 2 + 0.5f;
}

void Square::setPosition(float x, float y)
{
    o_x = x;
    o_y = y;
}

void Square::setVelocity(float u, float v)
{
    o_u = u;
    o_v = v;
}

bool Square::isOnCell(float x, float y)
{
    float x1 = (o_x - o_w / 2);
    float x2 = (o_x + o_w / 2);
    float y1 = (o_y - o_h / 2);
    float y2 = (o_y + o_h / 2);

    return (x <= x2 && y <= y2 && y >= y1 && x >= x1);
}

/**
 *
 * @return Array[4]: 0: left x, 1: top y, 2: right x, 3: bottom y
 */
std::vector<float> Square::getPosition()
{
    std::vector<float> vec;
    vec.push_back(o_x - o_w / 2);
    vec.push_back(o_y + o_h / 2);
    vec.push_back(o_x + o_w / 2);
    vec.push_back(o_y - o_h / 2);
    return vec;
}

std::vector<float> Square::getVelocity()
{
    std::vector<float> vec;
    vec.push_back(o_u);
    vec.push_back(o_v);
    return vec;
}