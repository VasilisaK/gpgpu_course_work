#include "Functions.h"

CrossPoint::CrossPoint() {
	CrossPointX = 0;
	CrossPointY = 0;
	Icr = 0;
}

CrossPoint::CrossPoint(CrossPoint& point) {
	CrossPointX = point.CrossPointX;
	CrossPointY = point.CrossPointY;
	Icr = point.Icr;
}

CrossPoint::~CrossPoint() {
}

TrueCoords::TrueCoords() {
	TrueCoordX = 0;
	TrueCoordY = 0;
}

TrueCoords::TrueCoords(TrueCoords& coords) {
	TrueCoordX = coords.TrueCoordX;
	TrueCoordY = coords.TrueCoordY;
}

TrueCoords::~TrueCoords() {
}

CrossPoint Cross(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {

	double a1, a2, b1, b2;
	CrossPoint point;
	double eps = 1.0e-5;
	double x_cross = -900000.0;
	double y_cross = -900000.0;
	int pr = 0;
	int Icr;

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
		}
	}

	if (abs(x3 - x4) < eps && abs(x1 - x2) > eps && pr == 0) {
		a1 = (y2 - y1) / (x2 - x1);
		b1 = y1 - a1 * x1;

		x_cross = x3;
		y_cross = a1 * x_cross + b1;
		pr = 2;
	}

	if (abs(x3 - x4) > eps && abs(x1 - x2) < eps && pr == 0) {
		a2 = (y4 - y3) / (x4 - x3);
		b2 = y3 - a2 * x3;

		x_cross = x1;
		y_cross = a2 * x_cross + b2;
		pr = 2;
	}


	Icr = 0;

	if (pr == 2) {

		if (RDistance(x3, x_cross, y3, y_cross) + RDistance(x4, x_cross, y4, y_cross) - RDistance(x3, x4, y3, y4) < eps) {

			if (RDistance(x1, x_cross, y1, y_cross) + RDistance(x2, x_cross, y2, y_cross) - RDistance(x1, x2, y1, y2) < eps) {
				Icr = 1;
			}
		}
	}

	point.CrossPointX = x_cross;
	point.CrossPointY = y_cross;
	point.Icr = Icr;

	return point;

}

TrueCoords Reflect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double Xcross, double Ycross) {

	double x2_new, y2_new;
	double a, b, ap, bp, xp, yp, eps;

	TrueCoords NewCoords;

	eps = 1e-8;

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

	x2_new = 2.0 * xp - x2;
	y2_new = 2.0 * yp - y2;

	NewCoords.TrueCoordX = x2_new;
	NewCoords.TrueCoordY = y2_new;

	return NewCoords;
}

double RDistance(double x1, double x2, double y1, double y2) {
	double Distance;
	Distance = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	return Distance;
}
