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

#define FOR_EACH_CELL            \
	for (i = 1; i <= grid_size; i++)     \
	{                            \
		for (j = 1; j <= grid_size; j++) \
		{
#define END_FOR \
	}           \
	}

Square::Square(float x, float y, int height, int width, int N) : o_x(x), o_y(y), o_h(height), o_w(width), grid_size(N)
{
    float x1 = (o_x - o_w / 2) / grid_size;
    float x2 = (o_x + o_w / 2 ) / grid_size;
    float y1 = (o_y - o_h / 2 ) / grid_size;
    float y2 = (o_y + o_h / 2) / grid_size;
    coordinates.push_back(Vec2f(x1, y1));
    coordinates.push_back(Vec2f(x2, y1));
    coordinates.push_back(Vec2f(x2, y2));
    coordinates.push_back(Vec2f(x1, y2));
}

void Square::draw()
{
    glLineWidth(1.0f);
    float h = 1.0f / grid_size;
    int i, j;
    FOR_EACH_CELL
        if (isOnCell(i, j)){
            glBegin(GL_QUADS);
            glColor3f(0.3, 0.1, 0.2);
            // Bottom left to bottom right
            glVertex2f(i * h, j * h);
            glVertex2f((i + 1) * h, j * h);
            // Bottom right to top right
            glVertex2f((i + 1) * h, j * h);
            glVertex2f((i + 1) * h, (j + 1) * h);
            // Top right to top left
            glVertex2f((i + 1) * h, (j + 1) * h);
            glVertex2f(i * h, (j + 1) * h);
            // Top left to bottom left
            glVertex2f(i * h, (j + 1) * h);
            glVertex2f(i * h, j * h);
            glEnd();
        }
    END_FOR

    glBegin(GL_LINES);
    glColor3f(0.6, 0.2, 0.4);
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

    int i, j;
    FOR_EACH_CELL
        if (isOnCell(i, j)) {
            if (!isOnCell(i + 1, j)) {
                boundaries[IX(i, j)].b_right = true;
                boundaries[IX(i, j)].b_x = i;
                boundaries[IX(i, j)].b_y = j;
            }
            if (!isOnCell(i - 1, j)) {
                boundaries[IX(i, j)].b_left = true;
                boundaries[IX(i, j)].b_x = i;
                boundaries[IX(i, j)].b_y = j;
            }
            if (!isOnCell(i, j + 1)) {
                boundaries[IX(i, j)].b_top = true;
                boundaries[IX(i, j)].b_x = i;
                boundaries[IX(i, j)].b_y = j;
            }
            if (!isOnCell(i, j - 1)) {
                boundaries[IX(i, j)].b_bottom = true;
                boundaries[IX(i, j)].b_x = i;
                boundaries[IX(i, j)].b_y = j;
            }
        }
    END_FOR

}

void Square::setPosition(float x, float y)
{
    float diffX = (o_x - x) / grid_size;
    float diffY = (o_y - y) / grid_size;

    coordinates[0][0] -= diffX;
    coordinates[0][1] -= diffY;
    coordinates[1][0] -= diffX;
    coordinates[1][1] -= diffY;
    coordinates[2][0] -= diffX;
    coordinates[2][1] -= diffY;
    coordinates[3][0] -= diffX;
    coordinates[3][1] -= diffY;

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
    float c_x = (x - 1/2) / grid_size;
    float c_y = (y - 1/2) / grid_size;

    float sign;

    for (int i = 0; i < coordinates.size(); i++)
    {
        sign = ((c_x - coordinates[i][0]) * (coordinates[(i+1)%4][1] - coordinates[i][1])) - ((c_y - coordinates[i][1]) * (coordinates[(i+1)%4][0] - coordinates[i][0]));

        if (sign > 0) return false;
    }

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