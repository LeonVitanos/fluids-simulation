#include "AngularSpringForce.h"

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define PI 3.1415926535897932384626433832795

AngularSpringForce::AngularSpringForce(Particle *p1, Particle *p2, Particle *p3, double dist, double k_s, double k_d) : m_p1(p1), m_p2(p2), m_p3(p3), m_dist(dist), m_ks(k_s), m_kd(k_d) {}

void AngularSpringForce::calculate()
{
    /**
     * f_a = - [k_s (l_abs - r) + k_d * ((l_dot * l) / l_abs)] * l/l_abs
     * But now the r is the distance according to the cosine rule
     */
    // Basically just the spring force, except we now use the cosine rule for distance
    Vec2f p1p2 = m_p1->m_Position - m_p2->m_Position;
    Vec2f p2p3 = m_p2->m_Position - m_p3->m_Position;
    double angle = (p1p2 * p2p3) / (norm(p1p2) * norm(p2p3));

    // Numerical errors, since cosine allows for at most 1 or -1
    angle = angle > 1 ? 1 : angle;
    angle = angle < -1 ? -1 : angle;
    angle = acos(angle);
    Vec2f p1p3 = m_p1->m_Position - m_p3->m_Position;
    Vec2f v1v3 = m_p1->m_Velocity - m_p3->m_Velocity;

    // Cosine rule distance computation
    auto cos_dist = sqrt(norm(p1p2) * norm(p1p2) + norm(p2p3) * norm(p2p3) - 2 * norm(p1p2) * norm(p2p3) * cos(m_dist));

    // Final spring computation
    Vec2f force = (m_ks * (norm(p1p3) - cos_dist) + m_kd * (p1p3 * v1v3) / norm(p1p3)) * (p1p3 / norm(p1p3));

    this->force = force;

    m_p1->m_Force += force - m_p1->m_Velocity * 0.1;
    m_p3->m_Force += force - m_p1->m_Velocity * 0.1;
}

void AngularSpringForce::draw(bool draw[])
{
    glBegin(GL_LINES);
    glColor3f(1.0f, 1.0f, 0.0f);
    glVertex2f(m_p1->m_Position[0], m_p1->m_Position[1]);
    glVertex2f(m_p2->m_Position[0], m_p2->m_Position[1]);
    glVertex2f(m_p3->m_Position[0], m_p3->m_Position[1]);
    glEnd();

    if (draw[6])
    {
        const double h = 0.03;

        glBegin(GL_LINES);
        glColor3f(0.0, 0.6, 1.0);
        glVertex2f(m_p1->m_Position[0], m_p1->m_Position[1]);
        glColor3f(0.0, 0.6, 1.0);
        glVertex2f(m_p1->m_Position[0] + 50 * this->force[0], m_p1->m_Position[1] + 50 * this->force[1]);
        glEnd();

        glPushMatrix();
        glTranslatef(m_p1->m_Position[0] + 50 * this->force[0], m_p1->m_Position[1] + 50 * this->force[1], 0.0);
        glRotatef((atan2(this->force[1], this->force[0]) * 180) / 3.14159265, 0.0, 0.0, 1.0);
        glColor3f(0.0, 0.6, 1.0);
        glBegin(GL_TRIANGLES);
        glVertex2f(-h / 2.0, -h / 2.0);
        glVertex2f(h / 1.0, 0.0);
        glVertex2f(-h / 2.0, h / 2.0);
        glEnd();
        glPopMatrix();

        glBegin(GL_LINES);
        glColor3f(1.0, 0.6, 0.0);
        glVertex2f(m_p2->m_Position[0], m_p2->m_Position[1]);
        glColor3f(1.0, 0.6, 0.0);
        glVertex2f(m_p2->m_Position[0] + 50 * this->force[0], m_p2->m_Position[1] + 50 * this->force[1]);
        glEnd();

        glPushMatrix();
        glTranslatef(m_p2->m_Position[0] + 50 * this->force[0], m_p2->m_Position[1] + 50 * this->force[1], 0.0);
        glRotatef((atan2(this->force[1], this->force[0]) * 180) / 3.14159265, 0.0, 0.0, 1.0);
        glColor3f(1.0, 0.6, 0.0);
        glBegin(GL_TRIANGLES);
        glVertex2f(-h / 2.0, -h / 2.0);
        glVertex2f(h / 1.0, 0.0);
        glVertex2f(-h / 2.0, h / 2.0);
        glEnd();
        glPopMatrix();
    }
}