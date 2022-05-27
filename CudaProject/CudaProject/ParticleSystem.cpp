#include "ParticleSystem.h"
#include "Functions.h"
// #include "kernel.cuh"


ParticleSystem::ParticleSystem() {
	ParticlesDist = 100;
	SourceCoordX_1 = 0;
	SourceCoordY_1 = 0;
	SourceCoordX_2 = 0;
	SourceCoordY_2 = 0;
	Red_1 = 1;
	Green_1 = 1;
	Blue_1 = 1;
	Red_2 = 1;
	Green_2 = 1;
	Blue_2 = 1;
	dt = 0;
	Life = 0;
	BasketLevel = 0;
	ParticlesInBasketNeeded = 0;
	BasketCounter = 0;
}

ParticleSystem::ParticleSystem(double SourceCoordX_1_inp, double SourceCoordY_1_inp, double SourceCoordX_2_inp, double SourceCoordY_2_inp, double dt_inp, double Life_inp, double BasketLevel_inp, int ParticlesInBasket_inp) {
	ParticlesDist = 100;
	SourceCoordX_1 = SourceCoordX_1_inp;
	SourceCoordY_1 = SourceCoordY_1_inp;
	SourceCoordX_2 = SourceCoordX_2_inp;
	SourceCoordY_2 = SourceCoordY_2_inp;
	Red_1 = 1;
	Green_1 = 1;
	Blue_1 = 1;
	Red_2 = 1;
	Green_2 = 1;
	Blue_2 = 1;
	dt = dt_inp;
	Life = Life_inp;
	BasketLevel = BasketLevel_inp;
	ParticlesInBasketNeeded = ParticlesInBasket_inp;
	BasketCounter = 0;
}

void ParticleSystem::InitSystem(double Vx_mag, double Vy_mag, float r_1_inp, float g_1_inp, float b_1_inp, double Life_1_inp, float r_2_inp, float g_2_inp, float b_2_inp, double Life_2_inp) {
	
	Red_1 = r_1_inp;
	Green_1 = g_1_inp;
	Blue_1 = b_1_inp;
	Red_2 = r_2_inp;
	Green_2 = g_2_inp;
	Blue_2 = b_2_inp;
	
	// 1st-type particles
	for (int i = 0; i < MAX_PARTICLES; ++i)
		particles[i] = Particle(SourceCoordX_1, SourceCoordY_1, Vx_mag*0.01 * (rand() % 101), Vy_mag*0.01 * (rand() % 101), Red_1, Green_1, Blue_1, Life_1_inp);

	// 2nd-type particles
	for (int i = 0; i < MAX_PARTICLES; ++i)
		particles_type2[i] = Particle(SourceCoordX_2, SourceCoordY_2, Vx_mag * 0.01 * (rand() % 101), Vy_mag * 0.01 * (rand() % 101), Red_2, Green_2, Blue_2, Life_2_inp);
}

void ParticleSystem::UpdateSystem(Geometry Geom) {

	int i, j;
	double NewXCoord, NewYCoord, NewVX, NewVY;
	double NewXCoord_2, NewYCoord_2, NewVX_2, NewVY_2;
	double tempVx, tempVy;

/*
	// Checking wall interaction for 1st-type particles
	for (i = 0; i < MAX_PARTICLES; ++i) {

		NewXCoord = particles[i].GetParams().x + particles[i].GetParams().Vx * dt;
		NewYCoord = particles[i].GetParams().y + particles[i].GetParams().Vy * dt;
		NewVX = particles[i].GetParams().Vx;
		NewVY = particles[i].GetParams().Vy;

		particles[i].UpdateParticle(NewXCoord, NewYCoord, Geom, dt);

		if (particles[i].GetParams().Life <= 0)
			particles[i].UpdateLifeStatus(SourceCoordX_1, SourceCoordY_1, 0.01 * (rand() % 101), 0.01 * (rand() % 101), 0, 0.01 * (rand() % 101), 0, Life);
			
		if (particles[i].GetParams().y > BasketLevel) {
			BasketCounter += 1;
			particles[i].UpdateLifeStatus(SourceCoordX_1, SourceCoordY_1, 0.01 * (rand() % 101), 0.01 * (rand() % 101), 0, 0.01 * (rand() % 101), 0, Life);
		}

	}
*/

	Calc(particles, Geom, MAX_PARTICLES);
	
	// Checking wall interaction for 2nd-type particles
	for (i = 0; i < MAX_PARTICLES; ++i) {

		NewXCoord_2 = particles_type2[i].GetParams().x + particles_type2[i].GetParams().Vx * dt;
		NewYCoord_2 = particles_type2[i].GetParams().y + particles_type2[i].GetParams().Vy * dt;
		NewVX_2 = particles_type2[i].GetParams().Vx;
		NewVY_2 = particles_type2[i].GetParams().Vy;

		particles_type2[i].UpdateParticle(NewXCoord_2, NewYCoord_2, Geom, dt);

		if (particles_type2[i].GetParams().Life <= 0)
			particles_type2[i].UpdateLifeStatus(SourceCoordX_2, SourceCoordY_2, 0.01 * (rand() % 101), 0.01 * (rand() % 101), 0, 0, 0.01 * (rand() % 101), Life);

		if (particles_type2[i].GetParams().y > BasketLevel) {
			BasketCounter += 1;
			particles_type2[i].UpdateLifeStatus(SourceCoordX_2, SourceCoordY_2, 0.01 * (rand() % 101), 0.01 * (rand() % 101), 0, 0, 0.01 * (rand() % 101), Life);
		}
			
	}

	// Checking particle-particle interaction
	for (i = 0; i < MAX_PARTICLES; ++i) {

		for (j = 0; j < MAX_PARTICLES; ++j) {

			ParticlesDist = RDistance(particles[i].GetParams().x, particles_type2[j].GetParams().x, particles[i].GetParams().y, particles_type2[j].GetParams().y);

			if (ParticlesDist <= 14.0) {

				tempVx = particles[i].GetParams().Vx;
				tempVy = particles[i].GetParams().Vy;

				particles[i].SetVelocity(particles_type2[j].GetParams().Vx, particles_type2[j].GetParams().Vy);
				particles_type2[j].SetVelocity(tempVx, tempVy);

			}
		}
	}


}

void ParticleSystem::GetSystem(GLfloat* vertices) {
	for (int i = 0; i < MAX_PARTICLES; ++i) {
		vertices[i * 2] = float(particles[i].GetParams().x);
		vertices[i * 2 + 1] = float(particles[i].GetParams().y);
		vertices[2 * MAX_PARTICLES + i * 2] = float(particles_type2[i].GetParams().x);
		vertices[2 * MAX_PARTICLES + i * 2 + 1] = float(particles_type2[i].GetParams().y);
	}
}

/*
std::vector<double> ParticleSystem::GetColor() {
	std::vector<double> Colors(6);
	Colors[0] = Red_1;
	Colors[1] = Green_1;
	Colors[2] = Blue_1;
	Colors[3] = Red_2;
	Colors[4] = Green_2;
	Colors[5] = Blue_2;
	return Colors;
}



void ParticleSystem::GetParticles(Particle* p) {
	for (int i = 0; i < MAX_PARTICLES; ++i)
		p[i] = particles[i];
}

void ParticleSystem::SetParticles(Particle* p) {
	for (int i = 0; i < MAX_PARTICLES; ++i)
		particles[i] = p[i];
}
*/

ParticleSystem::~ParticleSystem() {
}
