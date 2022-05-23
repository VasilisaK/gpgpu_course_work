#include "Particle.h"
#include "Functions.h"

Geometry::Geometry() {

	pointsCounter = 8;
	int i;
	std::vector<double> r(2);

	CoordLeft.resize(pointsCounter);
	for (i = 0; i < CoordLeft.size(); ++i)
		CoordLeft[i].resize(2);
	CoordRight.resize(pointsCounter);
	for (i = 0; i < CoordRight.size(); ++i)
		CoordRight[i].resize(2);
	JFaceVector.resize(pointsCounter);
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

		JFaceVector[i][0] = -r[1]; // JFaceVector = r rotated on -90 degree
		JFaceVector[i][1] = r[0]; // JFaceVector directed to increasing J - index

	}

}

Geometry::~Geometry() {
}

Particle::Particle() {
	XCoord = 0.0;
	YCoord = 0.0;
	VX = 0.0;
	VY = 0.0;
	r = 0;
	g = 0;
	b = 0;
	Lifetime = 0;
}

Particle::Particle(double XCoord_inp, double YCoord_inp, double VX_inp, double VY_inp, float r_inp, float g_inp, float b_inp, double Life_inp) {
	XCoord = XCoord_inp;
	YCoord = YCoord_inp;
	VX = VX_inp;
	VY = VY_inp;
	r = r_inp;
	g = g_inp;
	b = b_inp;
	Lifetime = Life_inp;
}

void Particle::UpdateLifeStatus(double XCoord_inp, double YCoord_inp, double VX_inp, double VY_inp, float r_inp, float g_inp, float b_inp, double Life_inp) {
	XCoord = XCoord_inp;
	YCoord = YCoord_inp;
	VX = VX_inp;
	VY = VY_inp;
	r = r_inp;
	g = g_inp;
	b = b_inp;
	Lifetime = Life_inp;
}

void Particle::UpdateParticle(double x_m1, double y_m1, Geometry Geom_inp, double timestep) {

	double x_m, y_m, u_m, v_m, d, Vn, Vtau;
	std::vector<double> CrossCoords(2);
	std::vector<double> NewCoords(2);
	std::vector<double> nv(2);
	std::vector<double> tv(2);

	x_m = XCoord;
	y_m = YCoord;
	u_m = VX;
	v_m = VY;
	int counter = 0;
	double Icr1 = 0;
	int i = 0;
	double eps = 1e-6;


	while ((counter < Geom_inp.pointsCounter) && (Icr1 == 0)) {

		CrossCoords = Cross(x_m, y_m, x_m1, y_m1, Geom_inp.CoordLeft[counter][0], Geom_inp.CoordLeft[counter][1], Geom_inp.CoordRight[counter][0], Geom_inp.CoordRight[counter][1]);

		Icr1 = CrossCoords[2];

		if (Icr1 == 1) {


			if (abs(x_m1 - CrossCoords[0]) < eps && abs(y_m1 - CrossCoords[1]) < eps) {
				x_m1 += VX * 0.01;
				y_m1 += VY * 0.01;
			}

			NewCoords = Reflect(x_m, y_m, x_m1, y_m1, Geom_inp.CoordLeft[counter][0], Geom_inp.CoordLeft[counter][1], Geom_inp.CoordRight[counter][0], Geom_inp.CoordRight[counter][1], CrossCoords[0], CrossCoords[1]);
			XCoord = NewCoords[0];
			YCoord = NewCoords[1];

			d = RDistance(Geom_inp.JFaceVector[counter][0], 0.0, Geom_inp.JFaceVector[counter][1], 0.0);
			nv[0] = -Geom_inp.JFaceVector[counter][0] / d;
			nv[1] = -Geom_inp.JFaceVector[counter][1] / d;
			tv[0] = -nv[1];
			tv[1] = nv[0];

			Vn = u_m * nv[0] + v_m * nv[1];
			Vtau = u_m * tv[0] + v_m * tv[1];

			VX = Vtau * tv[0] - Vn * nv[0];
			VY = Vtau * tv[1] - Vn * nv[1];

		}
		else {
			XCoord = x_m1;
			YCoord = y_m1;
		}

		counter += 1;
	}

	Lifetime -= timestep;

}

std::vector<double> Particle::GetCoords() {
	std::vector<double> Coords(2);
	Coords[0] = XCoord;
	Coords[1] = YCoord;
	return Coords;
}

std::vector<double> Particle::GetVelocity() {
	std::vector<double> Vel(2);
	Vel[0] = VX;
	Vel[1] = VY;
	return Vel;
}

double Particle::GetLifetime() {
	return Lifetime;
}

void Particle::SetVelocity(double VX_inp, double VY_inp) {
	VX = VX_inp;
	VY = VY_inp;
}

void Particle::SetCoord(double XCoord_inp, double YCoord_inp) {
	XCoord = XCoord_inp;
	YCoord = YCoord_inp;
}

Particle::~Particle() {
}
