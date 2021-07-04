#include "RodConstraint.h"
#include "vector"

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

RodConstraint::RodConstraint(Particle *p1, Particle *p2, double dist) : m_p1(p1), m_p2(p2), m_dist(dist) {
  Vec2f force(0, 0);
  f1 = force;
  f2 = force;
}

void RodConstraint::calculate()
{
  // If you thought the rod constraint was ugly
  Vec2f position1 = m_p1->m_Position;
  Vec2f position2 = m_p2->m_Position;
  Vec2f velocity1 = m_p1->m_Velocity;
  Vec2f velocity2 = m_p2->m_Velocity;
  Vec4f position = Vec4f(position1[0], position1[1], position2[0], position2[1]);
  Vec4f Q = Vec4f(m_p1->m_Force[0], m_p1->m_Force[1], m_p2->m_Force[0], m_p2->m_Force[1]);
  //Mat4 W = Mat4(Vec4f(1 / m_p1->m_Mass, 0, 0, 0), Vec4f(0, 1 / m_p1->m_Mass, 0, 0),
  //Vec4f(0, 0, 1 / m_p2->m_Mass, 0), Vec4f(0, 0, 0, 1 / m_p2->m_Mass));

  Vec4f J = Vec4f(2 * (position1[0] - position2[0]), 2 * (position1[1] - position2[1]),
                  2 * (position2[0] - position1[0]), 2 * (position2[1] - position1[1]));
  Vec4f Jdot = Vec4f(2 * (velocity1[0] - velocity2[0]), 2 * (velocity1[1] - velocity2[1]),
                     2 * (velocity2[0] - velocity1[0]), 2 * (velocity2[1] - velocity1[1]));

  double C = (position1 - position2) * (position1 - position2) - m_dist * m_dist;
  double Cdot = (position1 - position2) * (velocity1 - velocity2);
  //float JWJ = J[0] * W[0][0] * J[0] + J[1] * W[1][1] * J[1] + J[2] * W[2][2] * J[2] + J[3] * W[3][3] * J[3];
  //float JWQ = J[0] * W[0][0] * Q[0] + J[1] * W[1][1] * Q[1] + J[2] * W[2][2] * Q[2] + J[3] * W[3][3] * Q[3];
  //Since both particles have the same mass
  float JWJ = J * (1 / m_p1->m_Mass) * J;
  float JWQ = J * (1 / m_p1->m_Mass) * Q;
  float Jdotqdot = Jdot[0] * velocity1[0] + Jdot[1] * velocity1[1] + Jdot[2] * velocity2[0] + Jdot[3] * velocity2[1];
  float lambda = (0 - Jdotqdot - JWQ - C - Cdot) / JWJ;

  f1[0] = lambda * J[0];
  f1[1] = lambda * J[1];
  f2[0] = lambda * J[2];
  f2[1] = lambda * J[3];

  m_p1->m_Force[0] += lambda * J[0];
  m_p1->m_Force[1] += lambda * J[1];
  m_p2->m_Force[0] += lambda * J[2];
  m_p2->m_Force[1] += lambda * J[3];
}

void RodConstraint::draw(bool draw[])
{
  glBegin(GL_LINES);
  glColor3f(0.8, 0.7, 0.6);
  glVertex2f(m_p1->m_Position[0], m_p1->m_Position[1]);
  glColor3f(0.8, 0.7, 0.6);
  glVertex2f(m_p2->m_Position[0], m_p2->m_Position[1]);
  glEnd();

  if (draw[2])
  {
    const double h = 0.03;

    glBegin(GL_LINES);
    glColor3f(0.0, 1.0, 1.0);
    glVertex2f(m_p1->m_Position[0], m_p1->m_Position[1]);
    glColor3f(0.0, 1.0, 1.0);
    glVertex2f(m_p1->m_Position[0] + 50 * f1[0], m_p1->m_Position[1] + 50 * f1[1]);
    glEnd();

    glPushMatrix();
    glTranslatef(m_p1->m_Position[0] + 50 * f1[0], m_p1->m_Position[1] + 50 * f1[1], 0.0);
    glRotatef((atan2(f1[1], f1[0]) * 180) / 3.14159265, 0.0, 0.0, 1.0);
    glColor3f(0.0, 1.0, 1.0);
    glBegin(GL_TRIANGLES);
    glVertex2f(-h / 2.0, -h / 2.0);
    glVertex2f(h / 1.0, 0.0);
    glVertex2f(-h / 2.0, h / 2.0);
    glEnd();
    glPopMatrix();

    glBegin(GL_LINES);
    glColor3f(1.0, 1.0, 0.0);
    glVertex2f(m_p2->m_Position[0], m_p2->m_Position[1]);
    glColor3f(1.0, 1.0, 0.0);
    glVertex2f(m_p2->m_Position[0] + 50 * f2[0], m_p2->m_Position[1] + 50 * f2[1]);
    glEnd();

    glPushMatrix();
    glTranslatef(m_p2->m_Position[0] + 50 * f2[0], m_p2->m_Position[1] + 50 * f2[1], 0.0);
    glRotatef((atan2(f2[1], f2[0]) * 180) / 3.14159265, 0.0, 0.0, 1.0);
    glColor3f(1.0, 1.0, 0.0);
    glBegin(GL_TRIANGLES);
    glVertex2f(-h / 2.0, -h / 2.0);
    glVertex2f(h / 1.0, 0.0);
    glVertex2f(-h / 2.0, h / 2.0);
    glEnd();
    glPopMatrix();
  }
}
