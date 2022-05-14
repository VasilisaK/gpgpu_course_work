#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

#define MAX_PARTICLES 1
int pointsCounter = 8;
int i;
int particles_counter = 1;

std::vector<std::vector<double>> CoordLeft(pointsCounter);
std::vector<std::vector<double>> CoordRight(pointsCounter);
std::vector<std::vector<double>> JFaceVector(pointsCounter);

std::vector<std::vector<double>> ParticlesDist(particles_counter);

std::vector<double> r(2);
void IsOpen(std::ifstream& oi);
void IsOpen(std::ofstream& of);

// double Icr;

typedef struct
{
//		float life;                     // Продолжительность жизни частицы
//		float fade;                     // Быстрота гибели частицы
//		float x, y;                     // Позиция частицы по X и Y
//		float r, g, b;                  // Цвет частицы в RGB
//		float xi, yi;                   // Направление частицы по X и Y

	double XCoord;
	double YCoord;
	double VX;
	double VY;

}
particle;

particle particles[MAX_PARTICLES];
particle particles_type2[MAX_PARTICLES];

double Time = 5000;
double dt = 100;

std::vector<double> Cross(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
double RDistance(double x1, double x2, double y1, double y2);
std::vector<double> Reflect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double Xcross, double Ycross);
particle C_Boundary(particle p, double x_m1, double y_m1, std::vector<std::vector<double>> JFaceVector, std::vector<std::vector<double>> CoordLeft, std::vector<std::vector<double>> CoordRight);

int main()
{
	std::cout << "Start" << std::endl;

	// Initializing geometry

	for (i = 0; i < CoordLeft.size(); ++i)
		CoordLeft[i].resize(2);
	for (i = 0; i < CoordRight.size(); ++i)
		CoordRight[i].resize(2);
	for (i = 0; i < JFaceVector.size(); ++i)
		JFaceVector[i].resize(2);

	CoordLeft[0][0] = 0;
	CoordLeft[1][0] = 200;
	CoordLeft[2][0] = 350;
	CoordLeft[3][0] = 450;
	CoordLeft[4][0] = 450;
	CoordLeft[5][0] = 550;
	CoordLeft[6][0] = 600;
	CoordLeft[7][0] = 650;

	CoordLeft[0][1] = 200;
	CoordLeft[1][1] = 400;
	CoordLeft[2][1] = 400;
	CoordLeft[3][1] = 600;
	CoordLeft[4][1] = 400;
	CoordLeft[5][1] = 400;
	CoordLeft[6][1] = 250;
	CoordLeft[7][1] = 400;

	CoordRight[0][0] = 200;
	CoordRight[1][0] = 350;
	CoordRight[2][0] = 350;
	CoordRight[3][0] = 450;
	CoordRight[4][0] = 550;
	CoordRight[5][0] = 600;
	CoordRight[6][0] = 650;
	CoordRight[7][0] = 800;

	CoordRight[0][1] = 400;
	CoordRight[1][1] = 400;
	CoordRight[2][1] = 600;
	CoordRight[3][1] = 400;
	CoordRight[4][1] = 400;
	CoordRight[5][1] = 250;
	CoordRight[6][1] = 400;
	CoordRight[7][1] = 400;


	for (i = 0; i < CoordRight.size(); ++i) {

		r[0] = CoordRight[i][0] - CoordLeft[i][0]; // r = vector from one node to another
		r[1] = CoordRight[i][1] - CoordLeft[i][1]; 

		JFaceVector[i][0] = - r[1]; // JFaceVector = r rotated on -90 degree
		JFaceVector[i][1] = r[0]; // JFaceVector directed to increasing J - index
		
	}


	// Initializing particle properties

	particles[0].XCoord = 150.0;
	particles[0].YCoord = 130.0;
	particles[0].VX = 0.2;
	particles[0].VY = 0.1;
	
	particles_type2[0].XCoord = 230.0;
	particles_type2[0].YCoord = 130.0;
	particles_type2[0].VX = -0.2;
	particles_type2[0].VY = 0.1;

/*
	particles[0].XCoord = 150.0;
	particles[0].YCoord = 130.0;
	particles[0].VX = 0.2;
	particles[0].VY = 0.1;

	particles_type2[0].XCoord = 650.0;
	particles_type2[0].YCoord = 100.0;
	particles_type2[0].VX = -0.1;
	particles_type2[0].VY = 0.3;
*/
	/*
	particles[0].XCoord = 150.0;
	particles[0].YCoord = 130.0;
	particles[0].VX = 0.2;
	particles[0].VY = 0.1;

	particles_type2[0].XCoord = 700.0;
	particles_type2[0].YCoord = 300.0;
	particles_type2[0].VX = 0.0;
	particles_type2[0].VY = 0.1;
	*/
	for (i = 0; i < ParticlesDist.size(); ++i)
		ParticlesDist[i].resize(particles_counter);


	double CurrTime = 0;

	double NewXCoord, NewYCoord, NewVX, NewVY;
	double NewXCoord_2, NewYCoord_2, NewVX_2, NewVY_2;
	particle ParticleCorrected;
	particle ParticleCorrected_2;

	// Open file for patticle tracking

	std::ofstream OutputFile; //output unit
	char OutputFileName[] = "track.plt";

	OutputFile.open(OutputFileName);
	IsOpen(OutputFile);
	OutputFile << "VARIABLES = \"P1_x\", \"P1_y\"\n";


	std::ofstream OutputFile_2; //output unit
	char OutputFileName_2[] = "track_2.plt";

	OutputFile_2.open(OutputFileName_2);
	IsOpen(OutputFile_2);
	OutputFile_2 << "VARIABLES = \"P1_x\", \"P1_y\"\n";


	OutputFile << std::setprecision(8) << std::fixed << particles[0].XCoord << "\t" << std::setprecision(8) << std::fixed << particles[0].YCoord << "\n";
	OutputFile_2 << std::setprecision(8) << std::fixed << particles_type2[0].XCoord << "\t" << std::setprecision(8) << std::fixed << particles_type2[0].YCoord << "\n";


	while (CurrTime < Time) {

		NewXCoord = particles[0].XCoord + particles[0].VX * dt;
		NewYCoord = particles[0].YCoord + particles[0].VY * dt;
		NewVX = particles[0].VX;
		NewVY = particles[0].VY;

		NewXCoord_2 = particles_type2[0].XCoord + particles_type2[0].VX * dt;
		NewYCoord_2 = particles_type2[0].YCoord + particles_type2[0].VY * dt;
		NewVX_2 = particles_type2[0].VX;
		NewVY_2 = particles_type2[0].VY;

		std::cout << "Trying particle at time " << CurrTime << std::endl;

		double tempVx, tempVy;

		// Checking particle-particle interaction
		for (i = 0; i < ParticlesDist.size(); ++i) {
			ParticlesDist[i][0] = RDistance(NewXCoord, NewXCoord_2, NewYCoord, NewYCoord_2);
//			std::cout << "particle distance = " << ParticlesDist[i][0] << "\n";
			if (ParticlesDist[i][0] < 1e-6) {

				tempVx = NewVX;
				tempVy = NewVY;
				NewVX = NewVX_2;
				NewVY = NewVY_2;
				NewVX_2 = tempVx;
				NewVY_2 = tempVy;

				particles[0].VX = NewVX;
				particles[0].VY = NewVY;
				particles_type2[0].VX = NewVX_2;
				particles_type2[0].VY = NewVY_2;

				NewXCoord += NewVX*dt;
				NewYCoord += NewVY*dt;
				NewXCoord_2 += NewVX_2*dt;
				NewYCoord_2 += NewVY_2*dt;
				
			}
		}


		// Checking bounbary interaction

//		std::cout << "\nParticle 1\n";
		ParticleCorrected = C_Boundary(particles[0], NewXCoord, NewYCoord, JFaceVector, CoordLeft, CoordRight);
//		std::cout << "\nParticle 2\n";
		ParticleCorrected_2 = C_Boundary(particles_type2[0], NewXCoord_2, NewYCoord_2, JFaceVector, CoordLeft, CoordRight);

/*
		ParticleCorrected.XCoord = NewXCoord;
		ParticleCorrected.YCoord = NewYCoord;
		ParticleCorrected.VX = NewVX;
		ParticleCorrected.VY = NewVY;
*/
		OutputFile << std::setprecision(8) << std::fixed << ParticleCorrected.XCoord << "\t" << std::setprecision(8) << std::fixed << ParticleCorrected.YCoord << "\n";
		OutputFile_2 << std::setprecision(8) << std::fixed << ParticleCorrected_2.XCoord << "\t" << std::setprecision(8) << std::fixed << ParticleCorrected_2.YCoord << "\n";

		particles[0] = ParticleCorrected;
		particles_type2[0] = ParticleCorrected_2;

		CurrTime += dt;		

		std::cout << std::endl;

	}
	
	OutputFile.close();
	OutputFile_2.close();
	
	std::cout << "End calculating" << std::endl;

	return 0;
}

std::vector<double> Cross(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {

	double a1, a2, b1, b2;
	std::vector<double> CrossCoords(3);
	double eps = 1.0e-5;
	double x_cross = -900000.0;
	double y_cross = -900000.0;
	int pr = 0;
	double Icr;

	if (abs(x1 - x2) >= eps && abs(x3 - x4) >= eps) {
		a1 = (y2 - y1) / (x2 - x1);
		b1 = y1 - a1 * x1;
		a2 = (y4 - y3) / (x4 - x3);
		b2 = y3 - a2 * x3;
		pr = 1;
			if (abs(a1 - a2) > eps) {
				x_cross = (b2 - b1) / (a1 - a2);
				y_cross = a1 * x_cross + b1;
				pr = 2;
//				std::cout << "Cross1" << std::endl;
			}
	}

	if (abs(x3 - x4) < eps && abs(x1 - x2) > eps && pr == 0) {
		a1 = (y2 - y1) / (x2 - x1);
		b1 = y1 - a1 * x1;

		x_cross = x3;
		y_cross = a1 * x_cross + b1;
		pr = 2;
//		std::cout << "Cross2" << std::endl;
	}

	if (abs(x3 - x4) > eps && abs(x1 - x2) < eps && pr == 0) {
		a2 = (y4 - y3) / (x4 - x3);
		b2 = y3 - a2 * x3;

		x_cross = x1;
		y_cross = a2 * x_cross + b2;
		pr = 2;
//		std::cout << "Cross3" << std::endl;
	}

	
	Icr = 0;

	if (pr == 2) {
//		std::cout << "pr = 2 hooray" << std::endl;

		if (RDistance(x3, x_cross, y3, y_cross) + RDistance(x4, x_cross, y4, y_cross) - RDistance(x3, x4, y3, y4) < eps) {
//			std::cout << "Dist1" << std::endl;
//			std::cout << RDistance(x1, x_cross, y1, y_cross) + RDistance(x2, x_cross, y2, y_cross) - RDistance(x1, x2, y1, y2) << std::endl;
			
			if (RDistance(x1, x_cross, y1, y_cross) + RDistance(x2, x_cross, y2, y_cross) - RDistance(x1, x2, y1, y2) < eps) {
//			std::cout << "Dist2" << std::endl;
//			std::cout << RDistance(x1, x_cross, y1, y_cross) + RDistance(x2, x_cross, y2, y_cross) - RDistance(x1, x2, y1, y2) << std::endl;
				Icr = 1;
			} 
		}
	}
	


	CrossCoords[0] = x_cross;
	CrossCoords[1] = y_cross;
	CrossCoords[2] = Icr;

	return CrossCoords;

}



std::vector<double> Reflect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double Xcross, double Ycross) {
	
	double x2_new, y2_new;
	double a, b, ap, bp, xp, yp, eps;

	std::vector<double> NewCoords(2);

	eps = 1e-8;

//	std::cout << "x2 = " << x2 << ", y2 = " << y2 << " from Reflect\n";
/*
	if (x1 == Xcross && y1 == Ycross) {
		x2 = 10;
		y2 = 240;
	}
*/

	/*
	if ((x1 - Xcross < eps) && (y1 - Ycross < eps))
	{
		x2_new = x2;
		y2_new = y2;
	}
	else {
	*/


		if (abs(x3 - x4) >= eps) {
			a = (y4 - y3) / (x4 - x3);
			b = y3 - a * x3;
			if (abs(a) >= eps) {
				ap = -1. / a;
				bp = y2 - ap * x2;
				xp = (b - bp) / (ap - a);
				yp = ap * xp + bp;
			}
			else {
				xp = x2;
				yp = Ycross;
			}
		}
		else {
			xp = Xcross;
			yp = y2;
		}

		// Case xp = x2 and yp = y2 !

		x2_new = 2.0 * xp - x2;
		y2_new = 2.0 * yp - y2;
//	}

	NewCoords[0] = x2_new;
	NewCoords[1] = y2_new;

//	std::cout << "In Reflect function new coordinates\t" << NewCoords[0] << "\t" << NewCoords[1] << std::endl;

	return NewCoords;
}

double RDistance(double x1, double x2, double y1, double y2) {
	double Distance;
	Distance = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	return Distance;
}

particle C_Boundary(particle p, double x_m1, double y_m1, std::vector<std::vector<double>> JFaceVector, std::vector<std::vector<double>> CoordLeft, std::vector<std::vector<double>> CoordRight) {
	
	double x_m, y_m, u_m, v_m, d, Vn, Vtau;
	std::vector<double> CrossCoords(2);
	std::vector<double> NewCoords(2);
	std::vector<double> nv(2);
	std::vector<double> tv(2);
	particle particleTrue;
	x_m = p.XCoord;
	y_m = p.YCoord;
	u_m = p.VX;
	v_m = p.VY;
	int counter = 0;
	double Icr1 = 0;
	i = 0;
	double eps = 1e-6;

	particleTrue.XCoord = x_m1;
	particleTrue.YCoord = y_m1;
	particleTrue.VX = p.VX;
	particleTrue.VY = p.VY;

//	std::cout << particleTrue.XCoord << "\n";

//	std::cout << "\nfunction C_Boundary called" << std::endl;

		while ((counter < pointsCounter)&&(Icr1 == 0)) {

//			std::cout << "\nCheck crossection\n";

//			std::cout << "x_m = " << x_m << "\n";
//			std::cout << "y_m = " << y_m << "\n";
//			std::cout << "x_m1 = " << x_m1 << "\n";
//			std::cout << "y_m1 = " << y_m1 << "\n";
//			std::cout << "Previous step coordinates\t" << p.XCoord << "\t" << p.YCoord << std::endl;
//			std::cout << "Old coordinates 2\t" << particleTrue.XCoord << "\t" << particleTrue.YCoord << std::endl;
//			std::cout << "Old velocity\t" << particleTrue.VX << "\t" << particleTrue.VY << std::endl;
			CrossCoords = Cross(x_m, y_m, x_m1, y_m1, CoordLeft[counter][0], CoordLeft[counter][1], CoordRight[counter][0], CoordRight[counter][1]);
//			std::cout << "Bounbary # " << counter << std::endl;
//			std::cout << "Crossection coordinates\t" << CrossCoords[0] << "\t" << CrossCoords[1] << std::endl;
//			std::cout << "Particle velocity\t" << particleTrue.VX << "\t" << particleTrue.VY << std::endl;

			Icr1 = CrossCoords[2];

			if (Icr1 == 1) {


				if (abs(x_m1- CrossCoords[0]) < eps && abs(y_m1 - CrossCoords[1]) < eps) {
					x_m1 += p.VX * 0.01;
					y_m1 += p.VY * 0.01;
				}

//				std::cout << "Here is crossection\n";
				NewCoords = Reflect(x_m, y_m, x_m1, y_m1, CoordLeft[counter][0], CoordLeft[counter][1], CoordRight[counter][0], CoordRight[counter][1], CrossCoords[0], CrossCoords[1]);
				particleTrue.XCoord = NewCoords[0];
				particleTrue.YCoord = NewCoords[1];
//				std::cout << "New coordinates in C_Boundary\t" << NewCoords[0] << "\t" << NewCoords[1] << std::endl;

				d = RDistance(JFaceVector[counter][0], 0.0, JFaceVector[counter][1], 0.0);
				nv[0] = -JFaceVector[counter][0] / d;
				nv[1] = -JFaceVector[counter][1] / d;
				tv[0] = -nv[1];
				tv[1] = nv[0];

				Vn = u_m * nv[0] + v_m * nv[1];
				Vtau = u_m * tv[0] + v_m * tv[1];

				particleTrue.VX = Vtau * tv[0] - Vn * nv[0];
				particleTrue.VY = Vtau * tv[1] - Vn * nv[1];			

//				std::cout << "New particle velocity\t" << particleTrue.VX << "\t" << particleTrue.VY << std::endl;
			}
			else {
//				std::cout << "No crossection\n";
			}

			counter += 1;
		}

		return particleTrue;
}

void IsOpen(std::ifstream& oi)
{
	if (!oi.is_open())
	{
		exit(EXIT_FAILURE);
	}
}

void IsOpen(std::ofstream& of)
{
	if (!of.is_open())
	{
		exit(EXIT_FAILURE);
	}
}