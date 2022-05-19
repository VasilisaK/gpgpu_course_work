#include "Particle.h"
#include "Functions.h"


int main() {
	

#define MAX_PARTICLES 10

	int i, j;

	Particle particles[MAX_PARTICLES];
	Particle particles_type2[MAX_PARTICLES];
	double ParticlesDist;

	double SourceCoordX_1 = 200;
	double SourceCoordY_1 = 100;
	double SourceCoordX_2 = 600;
	double SourceCoordY_2 = 100;
	double Life = 1000;;
	double Time = 2000;
	double dt = 100;

	double CurrTime = 0;
	double NewXCoord, NewYCoord, NewVX, NewVY;
	double NewXCoord_2, NewYCoord_2, NewVX_2, NewVY_2;
	double tempVx, tempVy;
	int BasketCounter = 0;
	int ParticlesInBasket = 10;
	double BasketLevel = 400;

	std::cout << "Start" << std::endl;

	Geometry Geom;

	// Initializing particle system with properties
	// Particle(Coordinate_X, Coordinate_Y, Velocity_X, Velocity_Y, Color_Red, Color_Green, Color_Blue, Lifetime)

	// 1st-type particles, green
	for (i = 0; i < MAX_PARTICLES; ++i)
		particles[i] = Particle(SourceCoordX_1, SourceCoordY_1, 0.01 * (rand() % 101), 0.01 * (rand() % 101), 0, 0.01 * (rand() % 101), 0, Life);

	// 2nd-type particles, blue
	for (i = 0; i < MAX_PARTICLES; ++i)
		particles_type2[i] = Particle(SourceCoordX_2, SourceCoordY_2, 0.01 * (rand() % 101), 0.01 * (rand() % 101), 0, 0, 0.01 * (rand() % 101), Life);



	while (CurrTime < Time) {

		std::cout << "Time = " << CurrTime << std::endl;

		// Checking wall interaction for 1st-type particles
		for (i = 0; i < MAX_PARTICLES; ++i) {

			NewXCoord = particles[i].GetCoords()[0] + particles[i].GetVelocity()[0] * dt;
			NewYCoord = particles[i].GetCoords()[1] + particles[i].GetVelocity()[1] * dt;
			NewVX = particles[i].GetVelocity()[0];
			NewVY = particles[i].GetVelocity()[1];

			particles[i].UpdateParticle(NewXCoord, NewYCoord, Geom, dt);

			if (particles[i].GetLifetime() <= 0) 
				particles[i].UpdateLifeStatus(SourceCoordX_1, SourceCoordY_1, 0.01 * (rand() % 101), 0.01 * (rand() % 101), 0, 0.01 * (rand() % 101), 0, Life);

			if (particles[i].GetCoords()[1] > BasketLevel)
				BasketCounter += 1;

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
			
			if (particles_type2[i].GetCoords()[1] > BasketLevel)
				BasketCounter += 1;

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

		std::cout << BasketCounter << " particles in bascket\n";
		if (BasketCounter >= ParticlesInBasket)
			std::cout << "Victory!\n";

		CurrTime += dt;		

		std::cout << std::endl;

	}
	
	std::cout << "End calculating" << std::endl;

	return 0;
}
