#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

#include<iostream>
#include <cuda.h>


class CrossPoint {
public:
	double CrossPointX;
	double CrossPointY;
	int Icr;
	__device__ __host__ CrossPoint();
	__device__ __host__ CrossPoint(CrossPoint& point);
	__device__ __host__ ~CrossPoint();
};

class TrueCoords {
public:
	double TrueCoordX;
	double TrueCoordY;
	__device__ __host__ TrueCoords();
	__device__ __host__ TrueCoords(TrueCoords& point);
	__device__ __host__ ~TrueCoords();
};

__device__ __host__ CrossPoint Cross(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
__device__ __host__ TrueCoords Reflect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double Xcross, double Ycross);
__device__ __host__ double RDistance(double x1, double x2, double y1, double y2);

