#pragma once

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>


class Geometry {

public:
	friend class Particle;
	int pointsCounter;
	std::vector<std::vector<double>> CoordLeft;
	std::vector<std::vector<double>> CoordRight;
	std::vector<std::vector<double>> JFaceVector;
	__device__ __host__ Geometry();
	__device__ __host__ ~Geometry();
};

class Particle {
private:
	double XCoord, YCoord;
	double VX, VY;
	float r, g, b;
	double Lifetime;
public:
	__device__ __host__ Particle();
	__device__ __host__ Particle(double XCoord_inp, double YCoord_inp, double VX_inp, double VY_inp, float r_inp, float g_inp, float b_inp, double Life_inp);
	__device__ __host__ ~Particle();
	void UpdateParticle(double x_m1, double y_m1, Geometry Geom_inp, double timestep);
	void UpdateLifeStatus(double XCoord_inp, double YCoord_inp, double VX_inp, double VY_inp, float r_inp, float g_inp, float b_inp, double Life_inp);
	std::vector<double> GetCoords();
	std::vector<double> GetVelocity();
	double GetLifetime();
	void SetVelocity(double VX_inp, double VY_inp);
	void SetCoord(double XCoord_inp, double YCoord_inp);
};


