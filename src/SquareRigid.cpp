#include "SquareRigid.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

SquareRigid::SquareRigid(float x, float y, int height, int width, int N) : o_x(x), o_y(y), o_h(height), o_w(width), grid_size(N)
{
    /*rotation.push_back(0);
    rotation.push_back(0);
    rotation.push_back(0);
    rotation.push_back(0);*/
    rotation = 0.0f;
    float h = 1.0f / grid_size;
    x1 = (-o_w / 2) * h;
    x2 = (o_w / 2) * h;
    y1 = (-o_h / 2) * h;
    y2 = (o_h / 2) * h;
    coordinates.push_back(Vec2f(x1, y1));
    coordinates.push_back(Vec2f(x2, y1));
    coordinates.push_back(Vec2f(x1, y2));
    coordinates.push_back(Vec2f(x2, y2));
}

void SquareRigid::draw()
{
//    glPushMatrix();
    float h = 1.0f / grid_size;
//    glTranslatef(o_x * h, o_y * h, 0.0);
//    glRotatef(rotation, 0.0, 0.0, 1.0);
//
//    glLineWidth(1.0f);
//    glBegin(GL_LINES);
//    glColor3f(0.6, 0.2, 0.4);
//    glVertex2f(x1, y1);
//    glVertex2f(x2, y1);
//    glVertex2f(x1, y1);
//    glVertex2f(x1, y2);
//    glVertex2f(x2, y2);
//    glVertex2f(x1, y2);
//    glVertex2f(x2, y1);
//    glVertex2f(x2, y2);
//    glEnd();
//    glPopMatrix();
    glBegin(GL_LINES);

//    float h = 1.0f / grid_size;

    glColor3f(0.2, 0.6, 0.8);
    glVertex2f(coordinates[0][0] * h, coordinates[0][1] * h);
    glVertex2f(coordinates[1][0] * h, coordinates[1][1] * h);
    glVertex2f(coordinates[0][0] * h, coordinates[0][1] * h);
    glVertex2f(coordinates[2][0] * h, coordinates[2][1] * h);
    glVertex2f(coordinates[1][0] * h, coordinates[1][1] * h);
    glVertex2f(coordinates[3][0] * h, coordinates[3][1] * h);
    glVertex2f(coordinates[2][0] * h, coordinates[2][1] * h);
    glVertex2f(coordinates[3][0] * h, coordinates[3][1] * h);
    glEnd();
}

void SquareRigid::update(BoundaryCell *boundaries, float dt)
{
    float h = 1.0f / grid_size;

    rotation += 0.5f;
    for (int i = 0; i < coordinates.size(); i++)
    {
        coordinates[i][0] = o_x + ( - o_w / 2) * (float)std::cos(rotation) + (o_h / 2) * std::sin(rotation);
        coordinates[i][1] = o_y + (o_w / 2) * (float)std::sin(rotation) + (- o_h / 2) * std::cos(rotation);
        fprintf(stderr, "Coordinates: X=%g Y=%g rotation=%g i=%d\n",
                coordinates[i][0], coordinates[i][1], rotation, i);
    }
    //
    //    float h = 1.0f / grid_size;
    //    float x1 = (o_x - o_w / 2 + rotation[0]) * h;
    //    float x2 = (o_x + o_w / 2 + rotation[1]) * h;
    //    float y1 = (o_y - o_h / 2 + rotation[2]) * h;
    //    float y2 = (o_y + o_h / 2 + rotation[3]) * h;
}

void SquareRigid::setPosition(float x, float y)
{
    o_x = x;
    o_y = y;
}

void SquareRigid::setVelocity(float u, float v)
{
    o_u = u;
    o_v = v;
}

bool SquareRigid::isOnCell(float x, float y)
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
std::vector<float> SquareRigid::getPosition()
{
    std::vector<float> vec;
    vec.push_back(o_x);
    vec.push_back(o_y);
    /*vec.push_back(x1);
    vec.push_back(x2);
    vec.push_back(y1);
    vec.push_back(y2);*/
    return vec;
}

std::vector<float> SquareRigid::getVelocity()
{
    std::vector<float> vec;
    vec.push_back(o_u);
    vec.push_back(o_v);
    return vec;
}

std::vector<float> SquareRigid::getRotationMatrix(float angle)
{
    std::vector<float> vec;
    vec.push_back((float)std::cos(angle));
    vec.push_back((float)-std::sin(angle));
    vec.push_back((float)std::sin(angle));
    vec.push_back((float)std::cos(angle));

    return vec;
}