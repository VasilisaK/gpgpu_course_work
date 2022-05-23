
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

#include<iostream>
#include <cuda.h>

#include "kernel.cuh"

#include "ParticleSystem.h"

#define N 256

// Kernel Definition
__global__ void iter(Particle* p, Geometry g, int n)
{
	double NewXCoord, NewYCoord, NewVX, NewVY;
	double SourceCoordX_1 = 100;
	double SourceCoordY_1 = 200;
	double BasketLevel = 400;
	double Life = 1000;
	double dt = 1;
	int BasketCounter = 0;

	int i = threadIdx.x;

	if (i < n) {
		

		NewXCoord = p[i].GetCoords()[0] + p[i].GetVelocity()[0] * dt;
		NewYCoord = p[i].GetCoords()[1] + p[i].GetVelocity()[1] * dt;
		NewVX = p[i].GetVelocity()[0];
		NewVY = p[i].GetVelocity()[1];

		p[i].UpdateParticle(NewXCoord, NewYCoord, g, dt);

		if (p[i].GetLifetime() <= 0)
			p[i].UpdateLifeStatus(SourceCoordX_1, SourceCoordY_1, 0.01 * (rand() % 101), 0.01 * (rand() % 101), 0, 0.01 * (rand() % 101), 0, Life);

		if (p[i].GetCoords()[1] > BasketLevel) {
			BasketCounter += 1;
			p[i].UpdateLifeStatus(SourceCoordX_1, SourceCoordY_1, 0.01 * (rand() % 101), 0.01 * (rand() % 101), 0, 0.01 * (rand() % 101), 0, Life);
		}


	}
}

void Calc(Particle* h_a) {

	int i, j;
	double NewXCoord, NewYCoord, NewVX, NewVY;

	Geometry geom;

	// Allocate host memory

	// Initialize host array

	// Allocate arrays in Device memory
	Particle* d_a;
	Geometry g;
	
	cudaMalloc((void**)&d_a, MAX_PARTICLES * sizeof(Particle));
//	cudaMalloc((void)geom, sizeof(Geometry));

	// Copy memory from Host to Device
	cudaMemcpy(d_a, h_a, MAX_PARTICLES * sizeof(Particle), cudaMemcpyHostToDevice);
//	cudaMemcpy(g, geom, sizeof(Geometry), cudaMemcpyHostToDevice);

	// Block and Grid dimentions
	dim3 grid_size(1); dim3 block_size(N);

	// Launch Kernel
	iter << <grid_size, block_size >> > (d_a, g, N);

	// Some kind of synchronization
	cudaDeviceSynchronize();

	cudaMemcpy(h_a, d_a, N * sizeof(int), cudaMemcpyDeviceToHost);


	for (int i = 0; i < 10; ++i) {
		//		printf("c[%d] = %d\n", i, h_a[i]);
		//		printf("c[%d] = %d\n", i, h_b[i]);
		std::cout << "h_a: c[" << i << "] = " << h_a[i] << "\n";

	}



//	free(h_a);
	cudaFree(d_a);

	//	return 0;
}

/*
// Kernel Definition
__global__ void iter(int* a, int* b, int n)
{
	int i = threadIdx.x;
	if (i < n) {
		a[i] = a[i] * 2;
		b[i] = a[i] + 1;
	}
}

//void CalcFunction();

// int main() {
void Calc() {

	int* h_a;
	int* h_b;
	// Allocate host memory
	h_a = (int*)malloc(sizeof(int) * N);
	h_b = (int*)malloc(sizeof(int) * N);

	// Initialize host array
	for (int i = 0; i < N; i++) {
		h_a[i] = i;
		h_b[i] = i;
	}

	// Allocate arrays in Device memory
	int* d_a;
	int* d_b;
	cudaMalloc((void**)& d_a, N * sizeof(int));
	cudaMalloc((void**)& d_b, N * sizeof(int));

	// Copy memory from Host to Device
	cudaMemcpy(d_a, h_a, N * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, h_b, N * sizeof(int), cudaMemcpyHostToDevice);

	// Block and Grid dimentions
	dim3 grid_size(1); dim3 block_size(N);

	// Launch Kernel
	iter << <grid_size, block_size >> > (d_a, d_b, N);

	// Some kind of synchronization
	cudaDeviceSynchronize();

	cudaMemcpy(h_a, d_a, N * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_b, d_b, N * sizeof(int), cudaMemcpyDeviceToHost);

	for (int i = 0; i < 10; ++i) {
		//		printf("c[%d] = %d\n", i, h_a[i]);
		//		printf("c[%d] = %d\n", i, h_b[i]);
		std::cout << "h_a: c[" << i << "] = " << h_a[i] << "\n";
		std::cout << "h_b: c[" << i << "] = " << h_b[i] << "\n";
	}

	free(h_a);
	free(h_b);
	cudaFree(d_a);
	cudaFree(d_b);

	//	return 0;
}
*/