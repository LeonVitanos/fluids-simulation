#include "SpringForce.h"
#include <chrono>

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

SpringForce::SpringForce(Particle *p1, Particle *p2, double dist, double ks, double kd) : m_p1(p1), m_p2(p2), m_dist(dist), m_ks(ks), m_kd(kd) {}

void SpringForce::draw(bool draw[])
{
  glBegin(GL_LINES);
  glColor3f(0.0, 0.7, 0.8);
  glVertex2f(m_p1->m_Position[0], m_p1->m_Position[1]);
  glColor3f(1.0, 0.7, 0.2);
  glVertex2f(m_p2->m_Position[0], m_p2->m_Position[1]);
  glEnd();

  if (draw[1])
  {
    const double h = 0.03;

    glBegin(GL_LINES);
    glColor3f(0.0, 0.0, 1.0);
    glVertex2f(m_p1->m_Position[0], m_p1->m_Position[1]);
    glColor3f(0.0, 0.0, 1.0);
    glVertex2f(m_p1->m_Position[0] + 50 * -this->force[0], m_p1->m_Position[1] + 50 * -this->force[1]);
    glEnd();

    glPushMatrix();
    glTranslatef(m_p1->m_Position[0] + 50 * -this->force[0], m_p1->m_Position[1] + 50 * -this->force[1], 0.0);
    glRotatef((atan2(-this->force[1], -this->force[0]) * 180) / 3.14159265, 0.0, 0.0, 1.0);
    glColor3f(0.0, 0.0, 1.0);
    glBegin(GL_TRIANGLES);
    glVertex2f(-h / 2.0, -h / 2.0);
    glVertex2f(h / 1.0, 0.0);
    glVertex2f(-h / 2.0, h / 2.0);
    glEnd();
    glPopMatrix();

    glBegin(GL_LINES);
    glColor3f(1.0, 0.0, 0.0);
    glVertex2f(m_p2->m_Position[0], m_p2->m_Position[1]);
    glColor3f(1.0, 0.0, 0.0);
    glVertex2f(m_p2->m_Position[0] + 50 * this->force[0], m_p2->m_Position[1] + 50 * this->force[1]);
    glEnd();

    glPushMatrix();
    glTranslatef(m_p2->m_Position[0] + 50 * this->force[0], m_p2->m_Position[1] + 50 * this->force[1], 0.0);
    glRotatef((atan2(this->force[1], this->force[0]) * 180) / 3.14159265, 0.0, 0.0, 1.0);
    glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_TRIANGLES);
    glVertex2f(-h / 2.0, -h / 2.0);
    glVertex2f(h / 1.0, 0.0);
    glVertex2f(-h / 2.0, h / 2.0);
    glEnd();
    glPopMatrix();
  }
}

void SpringForce::calculate()
/**
 * f_a: forces on a
 * f_b: forces on b
 * l: a - b 
 * r: rest length
 * k_s: spring constant
 * k_d: dampening constant
 * l_dot: time derivative of l, difference between velocities of a and b
 * f_a = - [ k_s (l_abs - r) + k_d * ( (l_dot * l) / l_abs)] * l/l_abs.
 */
{
  Vec2f l = m_p1->m_Position - m_p2->m_Position;
  float dist = norm(l);
  Vec2f l_dot = m_p1->m_Velocity - m_p2->m_Velocity;
  float test = (m_ks * (dist - this->m_dist) + m_kd * (l_dot * l) / dist);
  this->force = test * (l / dist);
  m_p1->m_Force -= this->force;// - 0.01 * m_p1->m_Velocity;
  m_p2->m_Force += this->force;// - 0.01 * m_p2->m_Velocity;
}