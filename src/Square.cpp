#include "Square.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

Square::Square(float x, float y, int height, int width, int N) : o_x(x), o_y(y), o_h(height), o_w(width), grid_size(N) {

}

void Square::draw()
{
    glLineWidth(1.0f);
    glBegin(GL_QUADS);

    float h = 1.0f / grid_size;
    float x1 = (o_x - o_w / 2 - 0.5f) * h;
    float x2 = (o_x + o_w / 2 + 0.5f) * h;
    float y1 = (o_y - o_h / 2 - 0.5f) * h;
    float y2 = (o_y + o_h / 2 + 0.5f) * h;
    glColor3f(0.6, 0.2, 0.4);
    glVertex2f(x1, y1);
    glVertex2f(x1, y2);
    glVertex2f(x2, y2);
    glVertex2f(x2, y1);
    glEnd();
}

void Square::update()
{
    float h = 1.0f / grid_size;
    float x1 = (o_x - o_w / 2 - 0.5f) * h;
    float x2 = (o_x + o_w / 2 + 0.5f) * h;
    float y1 = (o_y - o_h / 2 - 0.5f) * h;
    float y2 = (o_y + o_h / 2 + 0.5f) * h;

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
    float h = 1.0f / grid_size;
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