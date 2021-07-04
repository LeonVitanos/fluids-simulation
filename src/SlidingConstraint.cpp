#include "SlidingConstraint.h"

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

SlidingConstraint::SlidingConstraint(Particle *p, double height) : m_p(p), m_height(height)
{
}

void SlidingConstraint::draw(bool draw[])
{
    if (draw[5]) {
        const double h = 0.03;

		glBegin(GL_LINES);
		glColor3f(0.6, 0.0, 1.0);
		glVertex2f(m_p->m_Position[0], m_p->m_Position[1]);
		glColor3f(0.6, 0.0, 1.0);
		glVertex2f(m_p->m_Position[0] + 50 * this->force[0], m_p->m_Position[1] + 50 * this->force[1]);
		glEnd();

		glPushMatrix();
		glTranslatef(m_p->m_Position[0] + 50 * this->force[0], m_p->m_Position[1] + 50 * this->force[1], 0.0);
		glRotatef((atan2(this->force[1], this->force[0]) * 180) / 3.14159265, 0.0, 0.0, 1.0);
		glColor3f(0.6, 0.0, 1.0);
		glBegin(GL_TRIANGLES);
		glVertex2f(-h / 2.0, -h / 2.0);
		glVertex2f(h / 1.0, 0.0);
		glVertex2f(-h / 2.0, h / 2.0);
		glEnd();
		glPopMatrix();
    }
}

void SlidingConstraint::calculate()
{

    Vec2f position = m_p->m_Position;
    Vec2f velocity = m_p->m_Velocity;
    Vec2f Q = m_p->m_Force;
    //Mat2 W = Mat2(Vec2f(1 / m_p->m_Mass, 0), Vec2f(0, 1 / m_p->m_Mass));

    Vec2f J = Vec2f(0, 1);
    Vec2f J_dot = Vec2f(0, 0);

    double C = position[1] - m_height;
    double C_dot = velocity[1];

    float J_dot_q_dot = J_dot * velocity;
    // Vector-Matrix-Vector multiplication so get ready for ugly code
    //float JWQ = J[0] * W[0][0] * Q[0] + J[1] * W[1][1] * Q[1];
    //float JWJ = J[0] * W[0][0] * J[0] + J[1] * W[1][1] * J[1];
    float JWQ = J * (1 / m_p->m_Mass) * Q;
    float JWJ = J * (1 / m_p->m_Mass) * J;
    float lambda = (-J_dot_q_dot - JWQ - C - C_dot) / JWJ;

    this->force[0] = lambda * J[0];
    this->force[1] = lambda * J[1];

    m_p->m_Force[0] += lambda * J[0];
    m_p->m_Force[1] += lambda * J[1];
}
