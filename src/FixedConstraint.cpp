#include "FixedConstraint.h"

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

FixedConstraint::FixedConstraint(Particle *particle)
{
    p = particle;
}

void FixedConstraint::draw(bool draw[]) {
    if (draw[4]) {
        const double h = 0.03;

		glBegin(GL_LINES);
		glColor3f(1.0, 0.0, 0.6);
		glVertex2f(p->m_Position[0], p->m_Position[1]);
		glColor3f(1.0, 0.0, 0.6);
		glVertex2f(p->m_Position[0] + 50 * this->force[0], p->m_Position[1] + 50 * this->force[1]);
		glEnd();

		glPushMatrix();
		glTranslatef(p->m_Position[0] + 50 * this->force[0], p->m_Position[1] + 50 * this->force[1], 0.0);
		glRotatef((atan2(this->force[1], this->force[0]) * 180) / 3.14159265, 0.0, 0.0, 1.0);
		glColor3f(1.0, 0.0, 0.6);
		glBegin(GL_TRIANGLES);
		glVertex2f(-h / 2.0, -h / 2.0);
		glVertex2f(h / 1.0, 0.0);
		glVertex2f(-h / 2.0, h / 2.0);
		glEnd();
		glPopMatrix();
    }
}

void FixedConstraint::calculate()
{
    this->force = -p->m_Force;
    p->m_Force = Vec2f(0.0, 0.0);
}