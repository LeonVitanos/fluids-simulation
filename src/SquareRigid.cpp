#include "SquareRigid.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define IX(i, j) ((i) + (grid_size + 2) * (j))

#define FOR_EACH_CELL            \
	for (i = 1; i <= grid_size; i++)     \
	{                            \
		for (j = 1; j <= grid_size; j++) \
		{
#define END_FOR \
	}           \
	}

SquareRigid::SquareRigid(float x, float y, int height, int width, int N) : o_x(x), o_y(y), o_h(height), o_w(width), grid_size(N)
{
    //o_x = 0.0f;
    //o_y = 0.0f;
    rotation = 0.0f;
    float h = 1.0f / grid_size;
    x1 = (o_x - o_w / 2) * h;
    x2 = (o_x + o_w / 2) * h;
    y1 = (o_y - o_h / 2) * h;
    y2 = (o_y + o_h / 2) * h;
    coordinates.push_back(Vec2f(x1, y1));
    coordinates.push_back(Vec2f(x2, y1));
    coordinates.push_back(Vec2f(x2, y2));
    coordinates.push_back(Vec2f(x1, y2));
    cells = (bool *)malloc((N+2) * (N+2) * sizeof(bool));
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
    // Bottom left to bottom right
    glVertex2f(coordinates[0][0], coordinates[0][1]);
    glVertex2f(coordinates[1][0], coordinates[1][1]);
    // Bottom right to top right
    glVertex2f(coordinates[1][0], coordinates[1][1]);
    glVertex2f(coordinates[2][0], coordinates[2][1]);
    // Top right to top left
    glVertex2f(coordinates[2][0], coordinates[2][1]);
    glVertex2f(coordinates[3][0], coordinates[3][1]);
    // Top left to bottom left
    glVertex2f(coordinates[3][0], coordinates[3][1]);
    glVertex2f(coordinates[0][0], coordinates[0][1]);
    glEnd();
}

void SquareRigid::update(BoundaryCell *boundaries, float dt)
{
    // Works fine the first time, core dumped second time.
    /*int size = (grid_size + 2) * (grid_size + 2);
    for (int i; i < size; i++) {
        cells[i] = false;
    }*/

    float h = 1.0f / grid_size;
    float tempX;
    float tempY;

    rotation = 0.05f;
    /*if (rotation > 360.0f) {
        rotation -= 360;
    }*/
    for (int i = 0; i < coordinates.size(); i++)
    {
        tempX = coordinates[i][0] - (o_x * h);
        tempY = coordinates[i][1] - (o_y * h);

        coordinates[i][0] = (tempX * std::cos(rotation) - tempY * std::sin(rotation)) + (o_x * h);
        coordinates[i][1] = (tempX * std::sin(rotation) + tempY * std::cos(rotation)) + (o_y * h);
        //fprintf(stderr, "Coordinates: X=%g Y=%g rotation=%g i=%d\n",
        //        coordinates[i][0], coordinates[i][1], rotation, i);
    }

    /*int i, j;
    FOR_EACH_CELL
        if (isOnCell(i, j)) cells[IX(i, j)] = true;
    END_FOR

    FOR_EACH_CELL
        if (!cells[IX(i + 1, j)]) {
            boundaries[IX(i, j)].b_right = true;
            boundaries[IX(i, j)].b_x = i;
            boundaries[IX(i, j)].b_y = j;
        }
        if (!cells[IX(i - 1, j)]) {
            boundaries[IX(i, j)].b_left = true;
            boundaries[IX(i, j)].b_x = i;
            boundaries[IX(i, j)].b_y = j;
        }
        if (!cells[IX(i, j + 1)]) {
            boundaries[IX(i, j)].b_top = true;
            boundaries[IX(i, j)].b_x = i;
            boundaries[IX(i, j)].b_y = j;
        }
        if (!cells[IX(i, j - 1)]) {
            boundaries[IX(i, j)].b_bottom = true;
            boundaries[IX(i, j)].b_x = i;
            boundaries[IX(i, j)].b_y = j;
        }
    END_FOR*/


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

    float c_x = (x + 1 / 2) * h;
    float c_y = (y + 1 / 2) * h;

    float sign;

    for (int i = 0; i < coordinates.size(); i++)
    {
        sign = ((c_x - coordinates[i][0]) * (coordinates[(i+1)%4][1] - coordinates[i][1])) - ((c_y - coordinates[i][1]) * (coordinates[(i+1)%4][0] - coordinates[i][0]));

        if (sign > 0) return false;
    }

    glBegin(GL_QUADS);
    glColor3f(0.1, 0.3, 0.4);
    // Bottom left to bottom right
    glVertex2f(x * h, y * h);
    glVertex2f((x + 1) * h, y * h);
    // Bottom right to top right
    glVertex2f((x + 1) * h, y * h);
    glVertex2f((x + 1) * h, (y + 1) * h);
    // Top right to top left
    glVertex2f((x + 1) * h, (y + 1) * h);
    glVertex2f(x * h, (y + 1) * h);
    // Top left to bottom left
    glVertex2f(x * h, (y + 1) * h);
    glVertex2f(x * h, y * h);
    glEnd();

    return true;
    /*float x1 = (o_x - o_w / 2);
    float x2 = (o_x + o_w / 2);
    float y1 = (o_y - o_h / 2);
    float y2 = (o_y + o_h / 2);

    return (x <= x2 && y <= y2 && y >= y1 && x >= x1);*/
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