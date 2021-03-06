/*
  ======================================================================
   demo.c --- protoype to show off the simple solver
  ----------------------------------------------------------------------
   Author : Jos Stam (jstam@aw.sgi.com)
   Creation Date : Jan 9 2003

   Description:

	This code is a simple prototype that demonstrates how to use the
	code provided in my GDC2003 paper entitles "Real-Time Fluid Dynamics
	for Games". This code uses OpenGL and GLUT for graphics and interface

  =======================================================================
*/

#include "SquareRigid.h"
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <vector>
#include "Particle.h"
#include "Force.h"
#include "SpringForce.h"
#include "Wall.h"

/* macros */

extern void simulation_step(int N, float *u, float *v, float *dens, std::vector<Particle *> pVector, std::vector<Force *> forces, std::vector<Force *> constraints, std::vector<Wall *> walls, float dt, int solver);

#define IX(i, j) ((i) + (N + 2) * (j))

#define FOR_EACH_CELL            \
	for (i = 1; i <= N; i++)     \
	{                            \
		for (j = 1; j <= N; j++) \
		{
#define END_FOR \
	}           \
	}

/* external definitions (from solver.c) */

extern void dens_step(int N, float *x, float *x0, float *u, float *v, float diff, float dt, BoundaryCell *boundaries, std::vector<BaseObject *> objects);
extern void vel_step(int N, float *u, float *v, float *u0, float *v0, float visc, float dt, float *d0, BoundaryCell *boundaries, std::vector<BaseObject *> objects);
extern void update_velocities(std::vector<BaseObject *> objects, float *u, float *v, float *d, int N);

/* global variables */

static int N;
static float dt, diff, visc;
static float force, source;

BoundaryCell *boundaries;
static int dvel;

static float *u, *v, *u_prev, *v_prev;
static float *dens, *dens_prev;
static float *curl;
std::vector<BaseObject *> objects;

static int win_id;
static int win_x, win_y;
static int mouse_down[3];
static int omx, omy, mx, my;

BaseObject *selectedObject;
std::vector<Particle *> pVector;
std::vector<Force *> forces;
std::vector<Force *> constraints;
std::vector<Wall *> walls;
bool drawForces[7] = {false, false, false, false, false, false, false};

/*
  ----------------------------------------------------------------------
   free/clear/allocate simulation data
  ----------------------------------------------------------------------
*/

static void free_data(void)
{
	if (u)
		free(u);
	if (v)
		free(v);
	if (u_prev)
		free(u_prev);
	if (v_prev)
		free(v_prev);
	if (dens)
		free(dens);
	if (dens_prev)
		free(dens_prev);
	if (curl)
		free(curl);
	if (boundaries)
		free(boundaries);
}

static void clear_data(void)
{
	int i, size = (N + 2) * (N + 2);

	for (i = 0; i < size; i++)
	{
		u[i] = v[i] = u_prev[i] = v_prev[i] = dens[i] = dens_prev[i] = curl[i] /*= boundaries[i]*/ = 0.0f;
	}

	int ii;
	size = pVector.size();

	for (ii = 0; ii < size; ii++)
	{
		pVector[ii]->reset();
	}
}

static int allocate_data(void)
{
	int size = (N + 2) * (N + 2);

	u = (float *)malloc(size * sizeof(float));
	v = (float *)malloc(size * sizeof(float));
	u_prev = (float *)malloc(size * sizeof(float));
	v_prev = (float *)malloc(size * sizeof(float));
	dens = (float *)malloc(size * sizeof(float));
	dens_prev = (float *)malloc(size * sizeof(float));
	curl = (float *)malloc(size * sizeof(float));
	boundaries = (BoundaryCell *)malloc(size * sizeof(BoundaryCell));

	if (!u || !v || !u_prev || !v_prev || !dens || !dens_prev || !curl)
	{
		fprintf(stderr, "cannot allocate data\n");
		return (0);
	}

	return (1);
}

/*
  ----------------------------------------------------------------------
   OpenGL specific drawing routines
  ----------------------------------------------------------------------
*/

static void pre_display(void)
{
	glViewport(0, 0, win_x, win_y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, 1.0, 0.0, 1.0);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
}

static void post_display(void)
{
	glutSwapBuffers();
}

static void draw_velocity(void)
{
	int i, j;
	float x, y, h;

	h = 1.0f / N;

	glColor3f(1.0f, 1.0f, 1.0f);
	glLineWidth(1.0f);

	glBegin(GL_LINES);

	for (i = 1; i <= N; i++)
	{
		x = (i - 0.5f) * h;
		for (j = 1; j <= N; j++)
		{
			y = (j - 0.5f) * h;

			glVertex2f(x, y);
			glVertex2f(x + u[IX(i, j)], y + v[IX(i, j)]);
		}
	}

	glEnd();
}

static void draw_density(void)
{
	int i, j;
	float x, y, h, d00, d01, d10, d11;

	h = 1.0f / N;

	glBegin(GL_QUADS);

	for (i = 0; i <= N; i++)
	{
		x = (i - 0.5f) * h;
		for (j = 0; j <= N; j++)
		{
			y = (j - 0.5f) * h;

			d00 = dens[IX(i, j)];
			d01 = dens[IX(i, j + 1)];
			d10 = dens[IX(i + 1, j)];
			d11 = dens[IX(i + 1, j + 1)];

			glColor3f(d00, d00, d00);
			glVertex2f(x, y);
			glColor3f(d10, d10, d10);
			glVertex2f(x + h, y);
			glColor3f(d11, d11, d11);
			glVertex2f(x + h, y + h);
			glColor3f(d01, d01, d01);
			glVertex2f(x, y + h);
		}
	}

	glEnd();
}

static void draw_particles(void)
{
	int size = pVector.size();

	for (int ii = 0; ii < size; ii++)
	{
		pVector[ii]->draw(false, false);
	}
}

static void draw_forces(void)
{
	int size = forces.size();

	for (int ii = 0; ii < size; ii++)
	{
		forces[ii]->draw(drawForces);
	}
}

/*
  ----------------------------------------------------------------------
   relates mouse movements to forces sources
  ----------------------------------------------------------------------
*/

static void get_from_UI(float *d, float *u, float *v)
{
	int i, j, size = (N + 2) * (N + 2);

	for (i = 0; i < size; i++)
	{
		u[i] = v[i] = d[i] = 0.0f;
	}

	if (!mouse_down[0] && !mouse_down[2])
	{
		selectedObject = NULL;
		return;
	}

	i = (int)((mx / (float)win_x) * N + 1);
	j = (int)(((win_y - my) / (float)win_y) * N + 1);

	if (i < 1 || i > N || j < 1 || j > N)
		return;

	if (mouse_down[0])
	{
		for (int k = 0; k < objects.size(); k++)
		{
			if (objects[k]->isOnCell(i, j))
			{
				objects[k]->setPosition(i, j);
				objects[k]->setVelocity(0, 0);
				selectedObject = objects[k];
			}
		}
		if (selectedObject != NULL)
		{
			selectedObject->setPosition(i, j);
			selectedObject->setVelocity((mx - omx), (omy - my));
		}
		else
		{
			u[IX(i, j)] = force * (mx - omx);
			v[IX(i, j)] = force * (omy - my);
		}
	}
	if (mouse_down[2])
	{
		d[IX(i, j)] = source;
	}

	omx = mx;
	omy = my;

	return;
}

/*
  ----------------------------------------------------------------------
   GLUT callback routines
  ----------------------------------------------------------------------
*/

static void key_func(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'c':
	case 'C':
		clear_data();
		break;

	case 'q':
	case 'Q':
		free_data();
		exit(0);
		break;

	case 'v':
	case 'V':
		dvel = !dvel;
		break;
	}
}

static void mouse_func(int button, int state, int x, int y)
{
	omx = mx = x;
	omx = my = y;

	mouse_down[button] = state == GLUT_DOWN;
}

static void motion_func(int x, int y)
{
	mx = x;
	my = y;
}

static void reshape_func(int width, int height)
{
	glutSetWindow(win_id);
	glutReshapeWindow(width, height);

	win_x = width;
	win_y = height;
}

static void idle_func(void)
{
	get_from_UI(dens_prev, u_prev, v_prev);
	int i, j;

	//Clear boundaries
	FOR_EACH_CELL
	boundaries[IX(i, j)].clear();
	END_FOR

	// Update object position and its boundaries
	for (int i = 0; i < objects.size(); i++)
	{
		objects[i]->update(boundaries, dt);
	}

	update_velocities(objects, u_prev, v_prev, dens, N);
	vel_step(N, u, v, u_prev, v_prev, visc, dt, curl, boundaries, objects);
	dens_step(N, dens, dens_prev, u, v, diff, dt, boundaries, objects);
	simulation_step(N, u, v, dens, pVector, forces, constraints, walls, dt, 2);

	glutSetWindow(win_id);
	glutPostRedisplay();
}

static void display_func(void)
{
	int i, j;
	float h = 1.0f / N;
	pre_display();

	if (dvel)
		draw_velocity();
	else
	{
		draw_density();
		/*FOR_EACH_CELL
		boundaries[IX(i, j)].draw();
		END_FOR*/
		for (int k = 0; k < objects.size(); k++)
		{
			FOR_EACH_CELL
			objects[k]->isOnCell(i, j);
			END_FOR
			objects[k]->draw();
		}
		draw_particles();
		draw_forces();
	}

	post_display();
}

/*
  ----------------------------------------------------------------------
   open_glut_window --- open a glut compatible window and set callbacks
  ----------------------------------------------------------------------
*/

static void open_glut_window(void)
{
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

	glutInitWindowPosition(0, 0);
	glutInitWindowSize(win_x, win_y);
	win_id = glutCreateWindow("Alias | wavefront");

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();

	pre_display();

	glutKeyboardFunc(key_func);
	glutMouseFunc(mouse_func);
	glutMotionFunc(motion_func);
	glutReshapeFunc(reshape_func);
	glutIdleFunc(idle_func);
	glutDisplayFunc(display_func);
}

/*
  ----------------------------------------------------------------------
   main --- main routine
  ----------------------------------------------------------------------
*/

int main(int argc, char **argv)
{
	glutInit(&argc, argv);

	if (argc != 1 && argc != 6)
	{
		fprintf(stderr, "usage : %s N dt diff visc force source\n", argv[0]);
		fprintf(stderr, "where:\n");
		fprintf(stderr, "\t N      : grid resolution\n");
		fprintf(stderr, "\t dt     : time step\n");
		fprintf(stderr, "\t diff   : diffusion rate of the density\n");
		fprintf(stderr, "\t visc   : viscosity of the fluid\n");
		fprintf(stderr, "\t force  : scales the mouse movement that generate a force\n");
		fprintf(stderr, "\t source : amount of density that will be deposited\n");
		exit(1);
	}

	if (argc == 1)
	{
		N = 64;
		dt = 0.1f;
		diff = 0.0f;
		visc = 0.0f;
		force = 5.0f;
		source = 100.0f;
		fprintf(stderr, "Using defaults : N=%d dt=%g diff=%g visc=%g force = %g source=%g\n",
				N, dt, diff, visc, force, source);
	}
	else
	{
		N = atoi(argv[1]);
		dt = atof(argv[2]);
		diff = atof(argv[3]);
		visc = atof(argv[4]);
		force = atof(argv[5]);
		source = atof(argv[6]);
	}

	printf("\n\nHow to use this demo:\n\n");
	printf("\t Add densities with the right mouse button\n");
	printf("\t Add velocities with the left mouse button and dragging the mouse\n");
	printf("\t Toggle density/velocity display with the 'v' key\n");
	printf("\t Clear the simulation by pressing the 'c' key\n");
	printf("\t Quit by pressing the 'q' key\n");

	dvel = 0;

	if (!allocate_data())
		exit(1);
	clear_data();

	// Square
//	objects.push_back((BaseObject *)new Square(20, 20, 10, 10, N));

	// Rigid Square
//	objects.push_back((BaseObject *)new SquareRigid(20, 20, 20, 20, N));
//
//	// Cloth
	const double dist = 0.05;
	Vec2f center(0.0, 0.0);
	const Vec2f offset(dist, 0.0);

	for (int i = 0; i < 6; i++)
	{
		center = Vec2f(0.5, 0.5 + 5 * dist - dist * i);
		for (int j = -2; j <= 2; j++)
		{
			pVector.push_back(new Particle(center + j * offset));
			//forces.push_back((Force *)new Gravity(particles.back()));
			if (j != -2)
				//spring connecting horizontally
				forces.push_back((Force *)new SpringForce(pVector[i * 5 + j + 1], pVector[i * 5 + j + 2], dist, 0.7, 0.8));
			if (i != 0)
			{
				//spring connecting vertically
				forces.push_back((Force *)new SpringForce(pVector[i * 5 - 5 + j + 2], pVector[i * 5 + j + 2], dist, 0.7, 0.8));
				if (j != 2)
					//spring connecting diagonally to the right
					forces.push_back((Force *)new SpringForce(pVector[i * 5 + j + 2], pVector[i * 5 + j + 2 - 4], dist, 0.05, 0.8));
				if (j != -2)
					//spring connecting diagonally to the left
					forces.push_back((Force *)new SpringForce(pVector[i * 5 + j + 2], pVector[i * 5 + j + 2 - 6], dist, 0.05, 0.8));
			}
		}
	}

	int size = pVector.size();

	for (int ii = 0; ii < size; ii++)
	{
		pVector[ii]->reset();
	}
	// End Cloth

	win_x = 512;
	win_y = 512;
	open_glut_window();

	glutMainLoop();

	exit(0);
}