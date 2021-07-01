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
    float x1 = (b_x - 1 - 0.5f) * h;
    float x2 = (b_x - 1+ 0.5f) * h;
    float y1 = (b_y - 1- 0.5f) * h;
    float y2 = (b_y - 1+ 0.5f) * h;
    glColor3f(0.6, 0.2, 0.4);
    glVertex2f(x1, y1);
    glVertex2f(x2, y1);
    glVertex2f(x1, y1);
    glVertex2f(x1, y2);
    glEnd();
}
