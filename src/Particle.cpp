#include "Particle.h"

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

Particle::Particle(const Vec2f &ConstructPos) : m_ConstructPos(ConstructPos), m_Position(Vec2f(0.0, 0.0)), m_Velocity(Vec2f(0.0, 0.0)), m_Force(Vec2f(0.0, 0.0))
{
	m_Mass = 2.5;
}

Particle::~Particle(void)
{
}

void Particle::reset()
{
	m_Position = m_ConstructPos;
	m_Velocity = Vec2f(0.0, 0.0);
	m_Force = Vec2f(0.0, 0.0);
}
void Particle::draw(bool force, bool vel)
{
	const double h = 0.015;
	glColor3f(1.f, 1.f, 1.f);
	glBegin(GL_QUADS);
	glVertex2f(m_Position[0] - h / 2.0, m_Position[1] - h / 2.0);
	glVertex2f(m_Position[0] + h / 2.0, m_Position[1] - h / 2.0);
	glVertex2f(m_Position[0] + h / 2.0, m_Position[1] + h / 2.0);
	glVertex2f(m_Position[0] - h / 2.0, m_Position[1] + h / 2.0);
	glEnd();

	if (force) {
		const double h = 0.03;

		glBegin(GL_LINES);
		glColor3f(1.0, 1.0, 1.0);
		glVertex2f(m_Position[0], m_Position[1]);
		glColor3f(1.0, 1.0, 1.0);
		glVertex2f(m_Position[0] + 50 * m_Force[0], m_Position[1] + 50 * m_Force[1]);
		glEnd();

		glPushMatrix();
		glTranslatef(m_Position[0] + 50 * m_Force[0], m_Position[1] + 50 * m_Force[1], 0.0);
		glRotatef((atan2(m_Force[1], m_Force[0]) * 180) / 3.14159265, 0.0, 0.0, 1.0);
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_TRIANGLES);
		glVertex2f(-h / 2.0, -h / 2.0);
		glVertex2f(h / 1.0, 0.0);
		glVertex2f(-h / 2.0, h / 2.0);
		glEnd();
		glPopMatrix();
	}

	if (vel) {
		const double h = 0.03;

		glBegin(GL_LINES);
		glColor3f(0.6, 0.6, 0.6);
		glVertex2f(m_Position[0], m_Position[1]);
		glColor3f(0.6, 0.6, 0.6);
		glVertex2f(m_Position[0] + 50 * m_Velocity[0], m_Position[1] + 50 * m_Velocity[1]);
		glEnd();

		glPushMatrix();
		glTranslatef(m_Position[0] + 50 * m_Velocity[0], m_Position[1] + 50 * m_Velocity[1], 0.0);
		glRotatef((atan2(m_Velocity[1], m_Velocity[0]) * 180) / 3.14159265, 0.0, 0.0, 1.0);
		glColor3f(0.6, 0.6, 0.6);
		glBegin(GL_TRIANGLES);
		glVertex2f(-h / 2.0, -h / 2.0);
		glVertex2f(h / 1.0, 0.0);
		glVertex2f(-h / 2.0, h / 2.0);
		glEnd();
		glPopMatrix();
	}
}
