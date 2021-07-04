#include "Boundary.h"

#include <stdlib.h>
#include <stdio.h>

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

BoundaryCell::BoundaryCell(float x, float y, bool left, bool top, bool right, bool bottom) {

}

void BoundaryCell::draw(){
    glLineWidth(1.0f);
    glBegin(GL_LINES);
    float h = 1.0f / 64;
    float x1 = (b_x - 1 ) * h;
    float x2 = (b_x) * h;
    float y1 = (b_y - 1) * h;
    float y2 = (b_y) * h;
    glColor3f(0.6, 0.2, 0.4);
    glVertex2f(x1, y1);
    glVertex2f(x2, y1);
    glVertex2f(x1, y1);
    glVertex2f(x1, y2);
    glEnd();
}

void BoundaryCell::clear(){
    b_x = 0;
    b_y = 0;
    b_left = false;
    b_right = false;
    b_bottom = false;
    b_top = false;
}
