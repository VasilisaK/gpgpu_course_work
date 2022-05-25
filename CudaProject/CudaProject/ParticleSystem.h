#pragma once

#include <stdio.h>

#include <../opengl/glew.h>
#include <../opengl/freeglut.h>

//#include <../opengl/glew.h>
//#include <../opengl/freeglut.h>

#include "Particle.h"

#include "kernel.cuh"

#define MAX_PARTICLES 20

class ParticleSystem {
private:
	Particle particles[MAX_PARTICLES];
	Particle particles_type2[MAX_PARTICLES];
	double ParticlesDist;
	double SourceCoordX_1, SourceCoordY_1;
	double SourceCoordX_2, SourceCoordY_2;
	float Red_1, Green_1, Blue_1;
	float Red_2, Green_2, Blue_2;
	double dt;	
	double Life;
	double BasketLevel;
	int ParticlesInBasketNeeded;
	int BasketCounter;
public:
	ParticleSystem();
	~ParticleSystem();
	ParticleSystem(double SourceCoordX_1_inp, double SourceCoordY_1_inp, double SourceCoordX_2_inp, double SourceCoordY_2_inp, double dt_inp, double Life_inp, double BasketLevel_inp, int ParticlesInBasket_inp);
	void InitSystem(double Vx, double Vy, float r_1_inp, float g_1_inp, float b_1_inp, double Life_1_inp, float r_2_inp, float g_2_inp, float b_2_inp, double Life_2_inp);
	void UpdateSystem(Geometry Geom);
	void GetSystem(GLfloat* vertices);
	std::vector<double> GetColor();
	double GetLifetime() { return Life; };
	int GetParticlesInBasketNeeded(){ return ParticlesInBasketNeeded; };
	int GetBasketCounter() { return BasketCounter; };
	void GetParticles(Particle* p);
	void SetParticles(Particle* p);
//	void Calc(Particle* h_a, Geometry g, int n);
};

