#pragma once
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>


std::vector<double> Cross(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);

std::vector<double> Reflect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double Xcross, double Ycross);

double RDistance(double x1, double x2, double y1, double y2);
