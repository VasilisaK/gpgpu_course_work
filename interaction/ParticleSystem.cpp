#include "ParticleSystem.h"
#include "Functions.h"


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


	// Checking wall interaction for 1st-type particles
	for (i = 0; i < MAX_PARTICLES; ++i) {

		NewXCoord = particles[i].GetCoords()[0] + particles[i].GetVelocity()[0] * dt;
		NewYCoord = particles[i].GetCoords()[1] + particles[i].GetVelocity()[1] * dt;
		NewVX = particles[i].GetVelocity()[0];
		NewVY = particles[i].GetVelocity()[1];

		particles[i].UpdateParticle(NewXCoord, NewYCoord, Geom, dt);

		if (particles[i].GetLifetime() <= 0)
			particles[i].UpdateLifeStatus(SourceCoordX_1, SourceCoordY_1, 0.01 * (rand() % 101), 0.01 * (rand() % 101), 0, 0.01 * (rand() % 101), 0, Life);
			
		if (particles[i].GetCoords()[1] > BasketLevel) {
			BasketCounter += 1;
			particles[i].UpdateLifeStatus(SourceCoordX_1, SourceCoordY_1, 0.01 * (rand() % 101), 0.01 * (rand() % 101), 0, 0.01 * (rand() % 101), 0, Life);
		}

	}

	// Checking wall interaction for 2nd-type particles
	for (i = 0; i < MAX_PARTICLES; ++i) {

		NewXCoord_2 = particles_type2[i].GetCoords()[0] + particles_type2[i].GetVelocity()[0] * dt;
		NewYCoord_2 = particles_type2[i].GetCoords()[1] + particles_type2[i].GetVelocity()[1] * dt;
		NewVX_2 = particles_type2[i].GetVelocity()[0];
		NewVY_2 = particles_type2[i].GetVelocity()[1];

		particles_type2[i].UpdateParticle(NewXCoord_2, NewYCoord_2, Geom, dt);

		if (particles_type2[i].GetLifetime() <= 0)
			particles_type2[i].UpdateLifeStatus(SourceCoordX_2, SourceCoordY_2, 0.01 * (rand() % 101), 0.01 * (rand() % 101), 0, 0, 0.01 * (rand() % 101), Life);

		if (particles_type2[i].GetCoords()[1] > BasketLevel) {
			BasketCounter += 1;
			particles_type2[i].UpdateLifeStatus(SourceCoordX_2, SourceCoordY_2, 0.01 * (rand() % 101), 0.01 * (rand() % 101), 0, 0, 0.01 * (rand() % 101), Life);
		}
			
	}

	// Checking particle-particle interaction
	for (i = 0; i < MAX_PARTICLES; ++i) {

		for (j = 0; j < MAX_PARTICLES; ++j) {

			ParticlesDist = RDistance(particles[i].GetCoords()[0], particles_type2[j].GetCoords()[0], particles[i].GetCoords()[1], particles_type2[j].GetCoords()[1]);

			if (ParticlesDist <= 1.0) {

				tempVx = particles[i].GetVelocity()[0];
				tempVy = particles[i].GetVelocity()[1];

				particles[i].SetVelocity(particles_type2[j].GetVelocity()[0], particles_type2[j].GetVelocity()[1]);
				particles_type2[j].SetVelocity(tempVx, tempVy);

			}
		}
	}


}

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

void ParticleSystem::GetSystem(GLfloat* vertices) {
	for (int i = 0; i < MAX_PARTICLES; ++i) {
		vertices[i * 2] = float(particles[i].GetCoords()[0]);
		vertices[i * 2 + 1] = float(particles[i].GetCoords()[1]);
		vertices[2 * MAX_PARTICLES + i * 2] = float(particles_type2[i].GetCoords()[0]);
		vertices[2 * MAX_PARTICLES + i * 2 + 1] = float(particles_type2[i].GetCoords()[1]);
	}
}

ParticleSystem::~ParticleSystem() {
}
