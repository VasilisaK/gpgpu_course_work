
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

#include<iostream>
#include <cuda.h>

#include <curand.h>
#include <curand_kernel.h>

#include "kernel.cuh"

#include "ParticleSystem.h"
#include "Particle.h"

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

	Coords coords;

	curandState_t state;
	double numb;
	double RANGE = 1;

	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < n) {
		
		coords = p[i].GetCoords();
		NewXCoord = coords.x + coords.Vx * dt;
		NewYCoord = coords.y + coords.Vy * dt;
		NewVX = coords.Vx;
		NewVY = coords.Vy;

		p[i].UpdateParticle(NewXCoord, NewYCoord, g, dt);

		if (coords.y <= 0) {
			numb = (double)((curand_uniform(&state)*(RANGE+1)));
			p[i].UpdateLifeStatus(SourceCoordX_1, SourceCoordY_1, numb, numb, 0, numb, 0, Life);
		}

		if (coords.y > BasketLevel) {
			BasketCounter += 1;
			numb = (double)((curand_uniform(&state) * (RANGE + 1)));
			p[i].UpdateLifeStatus(SourceCoordX_1, SourceCoordY_1, numb, numb, 0, numb, 0, Life);
		}


	}
}

void Particle::Calc(Particle* h_a, Geometry g, int n) {

	int i, j;
	double NewXCoord, NewYCoord, NewVX, NewVY;

//	Geometry geom;

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
	
	iter<<<BLOCKS_NUMBER,THREADS_NUMBER>>>(d_a, g, n);
	// Launch Kernel
//	iter << <grid_size, block_size >> > (d_a, g, n);

	// Some kind of synchronization
	cudaDeviceSynchronize();

	cudaMemcpy(h_a, d_a, size, cudaMemcpyDeviceToHost);


//	free(h_a);
	cudaFree(d_a);

}