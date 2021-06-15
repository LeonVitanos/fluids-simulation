#define IX(i, j) ((i) + (N + 2) * (j))
#define SWAP(x0, x)      \
	{                    \
		float *tmp = x0; \
		x0 = x;          \
		x = tmp;         \
	}
#define FOR_EACH_CELL            \
	for (i = 1; i <= N; i++)     \
	{                            \
		for (j = 1; j <= N; j++) \
		{
#define END_FOR \
	}           \
	}

#include <math.h>
#include "Particle.h"
#include "Force.h"
#include "Wall.h"
#include <chrono>
#include "ConstraintSolver.h"

void add_source(int N, float *x, float *s, float dt)
{
	int i, size = (N + 2) * (N + 2);
	for (i = 0; i < size; i++)
		x[i] += dt * s[i];
}

/**
 * @param (int) N: number of blocks in the grid
 * @param (int) b: in {1,2} => (1 -> boundary on the right, 2-> boundary on the left)
 * @param (float) x: place on the grid (counting from top left to bottom right)
 * @param (float) boundaries: place on the grid where boundaries should be
 */
void set_bnd(int N, int b, float *x, float *boundaries)
{
	int i, j;

	for (i = 1; i <= N; i++)
	{
		x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
		x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
		x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
		x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
	}
	// Create a boundary each block where this is defined
	// Here we most likely will want to recompute after forces are exerted on boundaries
	FOR_EACH_CELL
	if (boundaries[IX(i, j)] == 1)
	{
		x[IX(i, j)] = 0; // The center has no impact on any fluids
						 // TODO: The below comments should make sense but cause weird behaviour
						 // // On the left and right, invert the x as happens above
						 // x[IX(i - 1, j)] = b == 1 ? -x[IX(i - 2, j)] : x[IX(i - 2, j)];
						 // x[IX(i + 1, j)] = b == 1 ? -x[IX(i + 2, j)] : x[IX(i + 2, j)];
						 // // Also invert y vertically
						 // x[IX(i, j - 1)] = b == 2 ? -x[IX(i, j - 2)] : x[IX(i, j - 2)];
						 // x[IX(i, j + 1)] = b == 2 ? -x[IX(i, j + 2)] : x[IX(i, j + 2)];
	}
	END_FOR

	x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
	x[IX(0, N + 1)] = 0.5f * (x[IX(1, N + 1)] + x[IX(0, N)]);
	x[IX(N + 1, 0)] = 0.5f * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
	x[IX(N + 1, N + 1)] = 0.5f * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}

void lin_solve(int N, int b, float *x, float *x0, float a, float c, float *boundaries)
{
	int i, j, k;

	for (k = 0; k < 20; k++)
	{
		FOR_EACH_CELL
		x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
		END_FOR
		set_bnd(N, b, x, boundaries);
	}
}

void diffuse(int N, int b, float *x, float *x0, float diff, float dt, float *boundaries)
{
	float a = dt * diff * N * N;
	lin_solve(N, b, x, x0, a, 1 + 4 * a, boundaries);
}

void advect(int N, int b, float *d, float *d0, float *u, float *v, float dt, float *boundaries)
{
	int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, dt0;

	dt0 = dt * N;
	FOR_EACH_CELL
	x = i - dt0 * u[IX(i, j)];
	y = j - dt0 * v[IX(i, j)];
	if (x < 0.5f)
		x = 0.5f;
	if (x > N + 0.5f)
		x = N + 0.5f;
	i0 = (int)x;
	i1 = i0 + 1;
	if (y < 0.5f)
		y = 0.5f;
	if (y > N + 0.5f)
		y = N + 0.5f;
	j0 = (int)y;
	j1 = j0 + 1;
	s1 = x - i0;
	s0 = 1 - s1;
	t1 = y - j0;
	t0 = 1 - t1;
	d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
				  s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
	END_FOR
	set_bnd(N, b, d, boundaries);
}

void project(int N, float *u, float *v, float *p, float *div, float *boundaries)
{
	int i, j;

	FOR_EACH_CELL
	div[IX(i, j)] = -0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
	p[IX(i, j)] = 0;
	END_FOR
	set_bnd(N, 0, div, boundaries);
	set_bnd(N, 0, p, boundaries);

	lin_solve(N, 0, p, div, 1, 4, boundaries);

	FOR_EACH_CELL
	u[IX(i, j)] -= 0.5f * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
	v[IX(i, j)] -= 0.5f * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
	END_FOR
	set_bnd(N, 1, u, boundaries);
	set_bnd(N, 2, v, boundaries);
}

void vorticity_confinement(int N, float *u, float *v, float *u0, float *v0, float *d0, float dt)
{
	int size = (N + 2) * (N + 2);

	int i, j, ij;
	float *curlx = (float *)malloc(size * sizeof(float));
	float *curly = (float *)malloc(size * sizeof(float));
	float *curl = d0; //(float *)malloc(size * sizeof(float));
	float dt0 = dt * 5.0f;

	FOR_EACH_CELL
	ij = IX(i, j);
	curlx[ij] = (v[IX(i + 1, j)] - v[IX(i - 1, j)]) * 0.5f;
	curly[ij] = (u[IX(i, j + 1)] - u[IX(i, j - 1)]) * 0.5f;
	curl[ij] = abs(curly[ij] - curlx[ij]);
	END_FOR

	FOR_EACH_CELL
	ij = IX(i, j);
	float Nx = (curl[IX(i + 1, j)] - curl[IX(i - 1, j)]) * 0.5;
	float Ny = (curl[IX(i, j + 1)] - curl[IX(i, j - 1)]) * 0.5;
	float len = 1.0f / (sqrtf(Nx * Nx + Ny * Ny) + 0.0000001f);
	Nx *= len;
	Ny *= len;
	u0[ij] += Nx * (curlx[ij] - curly[ij]) * dt0;
	v0[ij] += Ny * (curly[ij] - curlx[ij]) * dt0;
	END_FOR
}

// Delp/Delt = -(u * D)p + kD^2p + S
void dens_step(int N, float *x, float *x0, float *u, float *v, float diff, float dt, float *boundaries)
{
	add_source(N, x, x0, dt); // S
	SWAP(x0, x);
	diffuse(N, 0, x, x0, diff, dt, boundaries); // kD^2p
	SWAP(x0, x);
	advect(N, 0, x, x0, u, v, dt, boundaries); // -(u * D)p
}

void vel_step(int N, float *u, float *v, float *u0, float *v0, float visc, float dt, float *d0, float *boundaries)
{
	add_source(N, u, u0, dt);
	add_source(N, v, v0, dt);
	vorticity_confinement(N, u, v, u0, v0, d0, dt);
	SWAP(u0, u);
	diffuse(N, 1, u, u0, visc, dt, boundaries);
	SWAP(v0, v);
	diffuse(N, 2, v, v0, visc, dt, boundaries);
	project(N, u, v, u0, v0, boundaries);
	SWAP(u0, u);
	SWAP(v0, v);
	advect(N, 1, u, u0, u0, v0, dt, boundaries);
	advect(N, 2, v, v0, u0, v0, dt, boundaries);
	project(N, u, v, u0, v0, boundaries);
}

void Clear_Forces(std::vector<Particle *> pVector)
{
	//Clear forces: zero the force accumulators
	for (auto &particle : pVector)
	{
		particle->m_Force = Vec2f(0.0, 0.0);
	}
}

void Compute_Forces(std::vector<Force *> forces)
{
	for (auto &force : forces)
	{
		force->calculate();
	}
}

void Compute_Collision(Particle *particle, Vec2f position, Vec2f velocity, std::vector<Wall *> walls, float dt)
{
	// Determine if particle crosses wall
	// Current particle coordinates (position vector)
	float p0_x = particle->m_Position[0];
	float p0_y = particle->m_Position[1];

	// Expected particle coordinates
	Vec2f newPos = position + dt * velocity;

	float p1_x = newPos[0];
	float p1_y = newPos[1];

	// Difference (direction vector)
	float s1_x = p1_x - p0_x;
	float s1_y = p1_y - p0_y;

	float mindt = -1;
	int iWall = -1;
	int wSize = walls.size();
	for (int i = 0; i < wSize; i++)
	{
		//Wall endpoint coordinates (p2 as position vector)
		float p2_x = walls[i]->m_p1[0];
		float p3_x = walls[i]->m_p2[0];
		float p2_y = walls[i]->m_p1[1];
		float p3_y = walls[i]->m_p2[1];

		// Difference (direction vector)
		float s2_x = p3_x - p2_x;
		float s2_y = p3_y - p2_y;

		float d = (-s2_x * s1_y + s1_x * s2_y);
		if (!(d == 0))
		{ // if not (near) parallel
			float s = (s1_x * (p0_y - p2_y) - s1_y * (p0_x - p2_x)) / d;
			float t = (s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / d;

			if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
			{ // if intersecting
				if (mindt < 0 || mindt < dt * t)
				{
					mindt = dt * t;
					iWall = i;
				}
			}
		}
	}

	if (mindt < 0)
	{													 // if no intersection
		particle->m_Position = position + dt * velocity; // xdot=v
	}
	else
	{ // if intersection
		// Position where particle hits wall
		newPos = position + (mindt - 0.05f) * velocity;

		// Calculate remaining velocity parallel to the wall
		Vec2f partVeloc = (dt - (mindt - 0.05f)) * velocity;
		Vec2f lineDir = walls[iWall]->m_p1 - walls[iWall]->m_p2;
		float length = sqrt(pow((walls[iWall]->m_p1[0] - walls[iWall]->m_p2[0]), 2) + pow((walls[iWall]->m_p1[1] - walls[iWall]->m_p2[1]), 2));
		lineDir = lineDir / length;
		particle->m_Velocity = ((partVeloc * lineDir) / (lineDir * lineDir)) * lineDir;

		// Final new position
		particle->m_Position = newPos + particle->m_Velocity;
	}
}

void simulation_step(std::vector<Particle *> pVector, std::vector<Force *> forces, std::vector<Force *> constraints, std::vector<Wall *> walls, float dt, int solver)
{
	int ii, size = pVector.size();

	std::vector<Particle *> initial;
	auto start = std::chrono::system_clock::now();

	switch (solver)
	{
	case 0: //Euler
		Clear_Forces(pVector);
		Compute_Forces(forces);
		ConstraintSolver(pVector, constraints).calculate();

		for (ii = 0; ii < size; ii++)
		{
			pVector[ii]->m_Velocity += dt * pVector[ii]->m_Force / pVector[ii]->m_Mass; // vdot = f/m

			Compute_Collision(pVector[ii], pVector[ii]->m_Position, pVector[ii]->m_Velocity, walls, dt);
		}

		break;
	case 1: //MidPoint
		Clear_Forces(pVector);
		Compute_Forces(forces);
		ConstraintSolver(pVector, constraints).calculate();

		for (ii = 0; ii < size; ii++)
		{ //Half euler step
			initial.push_back(pVector[ii]);
			pVector[ii]->m_Velocity += (dt / 2) * pVector[ii]->m_Force / pVector[ii]->m_Mass; // vdot = f/m
			//pVector[ii]->m_Position += (dt / 2) * pVector[ii]->m_Velocity;					  // xdot=v

			Compute_Collision(pVector[ii], pVector[ii]->m_Position, pVector[ii]->m_Velocity, walls, dt / 2);
		}

		Clear_Forces(pVector);
		Compute_Forces(forces);
		ConstraintSolver(pVector, constraints).calculate();

		for (ii = 0; ii < size; ii++)
		{
			pVector[ii]->m_Velocity = initial[ii]->m_Velocity + dt * pVector[ii]->m_Force / pVector[ii]->m_Mass; // vdot = f/m
			//pVector[ii]->m_Position = initial[ii]->m_Position + dt * pVector[ii]->m_Velocity;					 // xdot=v

			Compute_Collision(pVector[ii], initial[ii]->m_Position, pVector[ii]->m_Velocity, walls, dt);
		}

		break;
	case 2: //RungeKutta
	{
		std::vector<Vec2f> p1, p2, p3, p4, v1, v2, v3, v4;

		Clear_Forces(pVector);

		for (ii = 0; ii < size; ii++)
		{
			initial.push_back(pVector[ii]);
		}

		Compute_Forces(forces);
		ConstraintSolver(pVector, constraints).calculate();

		for (ii = 0; ii < size; ii++)
		{
			p1.push_back(pVector[ii]->m_Velocity);
			v1.push_back(pVector[ii]->m_Force / pVector[ii]->m_Mass);

			pVector[ii]->m_Velocity += dt / 2 * v1[ii]; // vdot = f/m
			//pVector[ii]->m_Position += dt / 2 * p1[ii]; // xdot=v

			Compute_Collision(pVector[ii], pVector[ii]->m_Position, p1[ii], walls, dt / 2);
		}

		Clear_Forces(pVector);
		Compute_Forces(forces);
		ConstraintSolver(pVector, constraints).calculate();

		for (ii = 0; ii < size; ii++)
		{
			p2.push_back(pVector[ii]->m_Velocity);
			v2.push_back(pVector[ii]->m_Force / pVector[ii]->m_Mass);

			pVector[ii]->m_Velocity = initial[ii]->m_Velocity + dt / 2 * v2[ii]; // vdot = f/m
			//pVector[ii]->m_Position = initial[ii]->m_Position + dt / 2 * p2[ii]; // xdot=v

			Compute_Collision(pVector[ii], initial[ii]->m_Position, p2[ii], walls, dt / 2);
		}

		Clear_Forces(pVector);
		Compute_Forces(forces);
		ConstraintSolver(pVector, constraints).calculate();

		for (ii = 0; ii < size; ii++)
		{
			p3.push_back(pVector[ii]->m_Velocity);
			v3.push_back(pVector[ii]->m_Force / pVector[ii]->m_Mass);

			pVector[ii]->m_Velocity = initial[ii]->m_Velocity + dt * v3[ii]; // vdot = f/m
			//pVector[ii]->m_Position = initial[ii]->m_Position + dt * p3[ii]; // xdot=v

			Compute_Collision(pVector[ii], initial[ii]->m_Position, p3[ii], walls, dt);
		}

		Clear_Forces(pVector);
		Compute_Forces(forces);
		ConstraintSolver(pVector, constraints).calculate();

		for (ii = 0; ii < size; ii++)
		{
			p4.push_back(pVector[ii]->m_Velocity);
			v4.push_back(pVector[ii]->m_Force / pVector[ii]->m_Mass);

			pVector[ii]->m_Velocity = initial[ii]->m_Velocity + 1 / 6 * v1[ii] + 1 / 3 * v2[ii] + 1 / 3 * v3[ii] + 1 / 6 * v4[ii]; // vdot = f/m
			//pVector[ii]->m_Position = initial[ii]->m_Position + 1 / 6 * p1[ii] + 1 / 3 * p2[ii] + 1 / 3 * p3[ii] + 1 / 6 * p4[ii]; // xdot=v

			Compute_Collision(pVector[ii], initial[ii]->m_Position, 1 / 6 * p1[ii] + 1 / 3 * p2[ii] + 1 / 3 * p3[ii] + 1 / 6 * p4[ii], walls, 1);
		}

		break;
	}
	case 3: // Implicit
		break;
	}

	//CPU timing
	/*auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	t.push_back(elapsed_seconds.count() * 1000);*/
	//std::cout << "Average " << accumulate( t.begin(), t.end(), 0.0) / t.size() << "ms\n";
}
