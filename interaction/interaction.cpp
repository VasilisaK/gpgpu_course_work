#include "ParticleSystem.h"
//#include "Functions.h"

#include <../../../opengl/glew.h>
#include <../../../opengl/freeglut.h>

int main() {


#define MAX_PARTICLES 10

	double Time = 2000;
	double CurrTime = 0;
	double dt = 100;

	std::cout << "Start" << std::endl;

	// Initializing geometry
	Geometry Geom;

	// Initializing particle system properties
	// Particle(SourceCoordX_1, SourceCoordY_1, SourceCoordX_2, SourceCoordY_2, dt, LifeTime, BasketLevel, ParticlesInBasketNeeded)
	ParticleSystem PlayingSystem;
	PlayingSystem = ParticleSystem(200, 100, 600, 100, dt, 1000, 400, 20);

	// Initialize particle system
	// InitSystem(Vx_magnitude, Vy_magnitude, Red_1, Green_1, Blue_1, LifeTime_1, Red_2, Green_2, Blue_2, LifeTime_2)
	PlayingSystem.InitSystem(1.0, 1.0, 0, 1, 0, 1000, 0, 0, 1, 1000);


	while (CurrTime < Time) {

		std::cout << "Time = " << CurrTime << std::endl;

		PlayingSystem.UpdateSystem(Geom);

		std::cout << PlayingSystem.GetBasketCounter() << " particles in basket\n";
		if (PlayingSystem.GetBasketCounter() >= PlayingSystem.GetParticlesInBasketNeeded())
			std::cout << "Victory!\n";		
		
		GLfloat* vertices = PlayingSystem.GetSystem();
		std::cout << "( " << vertices[0] << ", " << vertices[1] << " )\n";

		std::cout << "Color indexes type1: ( " << PlayingSystem.GetColor()[0] << ", " << PlayingSystem.GetColor()[1] << ", " << PlayingSystem.GetColor()[2] << " )\n";

		CurrTime += dt;

		std::cout << std::endl;

	}

	std::cout << "End calculating" << std::endl;

	return 0;
}
