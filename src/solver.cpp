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
