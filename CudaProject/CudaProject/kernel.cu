
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

#include<iostream>
#include <cuda.h>

#include "kernel.cuh"

#include "ParticleSystem.h"


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

	int i = blockDim.x * blockIdx.x + threadIdx.x;

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

void Calc(Particle* h_a, Geometry g, int n) {

	int i, j;
	double NewXCoord, NewYCoord, NewVX, NewVY;

	Geometry geom;

	// Allocate host memory

	// Initialize host array

	// Allocate arrays in Device memory
	Particle* d_a;
	size_t size= BLOCKS_NUMBER*THREADS_NUMBER*sizeof(Particle);

	cudaMalloc(&d_a, size);

	// Copy memory from Host to Device
	cudaMemcpy(d_a, h_a, size, cudaMemcpyHostToDevice);
//	cudaMemcpy(g, geom, sizeof(Geometry), cudaMemcpyHostToDevice);

	// Block and Grid dimentions
	
	iter<<<BLOCKS_NUMBER,THREADS_NUMBER>>>(dev_classarray);
	// Launch Kernel
	iter << <grid_size, block_size >> > (d_a, g, n);

	// Some kind of synchronization
	cudaDeviceSynchronize();

	cudaMemcpy(h_a, d_a, size, cudaMemcpyDeviceToHost);


//	free(h_a);
	cudaFree(d_a);

	//	return 0;
}