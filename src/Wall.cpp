#include "Wall.h"

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

Wall::Wall(const Vec2f p1, const Vec2f p2) : m_p1(p1), m_p2(p2)
{
}

void Wall::draw()
{
    glBegin(GL_LINES);
    glColor3f(0.8, 0.7, 0.6);
    glVertex2f(m_p1[0], m_p1[1]);
    glColor3f(0.8, 0.7, 0.6);
    glVertex2f(m_p2[0], m_p2[1]);
    glEnd();
}