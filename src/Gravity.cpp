#include "Gravity.h"

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

Gravity::Gravity(Particle *particle)
{
    p = particle;
}

Gravity::~Gravity() {}

void Gravity::draw(bool draw[])
{
    if (draw[0])
    {
        glBegin(GL_LINES);
        glColor3f(0.0, 1.0, 0.0);
        glVertex2f(p->m_Position[0], p->m_Position[1]);
        glColor3f(0.0, 1.0, 0.0);
        glVertex2f(p->m_Position[0], p->m_Position[1] + 50 * this->force[1]);
        glEnd();

        const double h = 0.03;
        glColor3f(0.0, 1.0, 0.0);
        glBegin(GL_TRIANGLES);
        glVertex2f(p->m_Position[0] - h / 2.0, p->m_Position[1] + 50 * this->force[1] + h / 2.0);
        glVertex2f(p->m_Position[0], p->m_Position[1] + 50 * this->force[1] - h);
        glVertex2f(p->m_Position[0] + h / 2.0, p->m_Position[1] + 50 * this->force[1] + h / 2.0);
        glEnd();
    }
}

void Gravity::calculate()
{
    this->force = Vec2f(0, -(p->m_Mass * 9.81f * 0.0001f));
    p->m_Force += this->force - p->m_Velocity * 0.01;
}
