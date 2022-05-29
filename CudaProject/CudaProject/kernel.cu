
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cooperative_groups.h"
#include <device_functions.h>
#include <cuda_device_runtime_api.h>

#include <stdio.h>

#include<iostream>
#include <cuda.h>

#include <curand.h>
#include <curand_kernel.h>

#include "kernel.cuh"

//#include "ParticleSystem.h"
#include "Particle.h"
#include "Functions.h"


// Kernel Definition
__global__ void iter(Particle* p1, Particle* p2, Geometry g, int n)
{
	double NewXCoord, NewYCoord, NewVX, NewVY;
	double SourceCoordX_1 = 300;
	double SourceCoordY_1 = 100;
	double SourceCoordX_2 = 500;
	double SourceCoordY_2 = 100;
	double BasketLevel = 400;
	double Life = 10000;
	double dt = 0.3;
	int BasketCounter = 0;

	ParticleParams params1;
	ParticleParams params2;

	curandState_t state;
	curand_init(0, 0, 0, &state);

	double numb;
	int RANGE = 1;
	int MAX = 101;

//	numb = curand_uniform(&state);
//	printf("%d\n", numb);

	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < n) {
		
		params1 = p1[i].GetParams();
		NewXCoord = params1.x + params1.Vx * dt;
		NewYCoord = params1.y + params1.Vy * dt;
		NewVX = params1.Vx;
		NewVY = params1.Vy;

		p1[i].UpdateParticle(NewXCoord, NewYCoord, g, dt);
		if (params1.Life <= 0) {
			numb = (double)((curand_uniform(&state)*(RANGE+1)));		
			p1[i].UpdateLifeStatus(SourceCoordX_1, SourceCoordY_1, 0.1, 0.1, 0, 1, 0, Life);
		}

		if (params1.y > BasketLevel) {
			BasketCounter += 1;
			numb = (double)((curand_uniform(&state) * (RANGE + 1)));
			p1[i].UpdateLifeStatus(SourceCoordX_1, SourceCoordY_1, 0.1, 0.1, 0, 1, 0, Life);
		}

		params2 = p2[i].GetParams();
		NewXCoord = params2.x + params2.Vx * dt;
		NewYCoord = params2.y + params2.Vy * dt;
		NewVX = params2.Vx;
		NewVY = params2.Vy;

		p2[i].UpdateParticle(NewXCoord, NewYCoord, g, dt);

		if (params2.Life <= 0) {
			numb = (double)((curand_uniform(&state) * (RANGE + 1)));
			p2[i].UpdateLifeStatus(SourceCoordX_2, SourceCoordY_2, 0.1, 0.1, 0, 0, 1, Life);
		}

		if (params2.y > BasketLevel) {
			BasketCounter += 1;
			numb = (double)((curand_uniform(&state) * (RANGE + 1)));
			p2[i].UpdateLifeStatus(SourceCoordX_2, SourceCoordY_2, 0.1, 0.1, 0, 0, 1, Life);
		}

	}

	__syncthreads();
}

__global__ void iter2(Particle* p1, Particle* p2, int n) {

	double tempVx, tempVy;
	double ParticlesDist;
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j = blockDim.y * blockIdx.y + threadIdx.y;

	// Checking particle-particle interaction

	if (i >= n || j >= n)
		return;

	ParticlesDist = RDistance(p1[i].GetParams().x, p2[j].GetParams().x, p1[i].GetParams().y, p2[j].GetParams().y);

	if (ParticlesDist <= 14.0) {

		tempVx = p1[i].GetParams().Vx;
		tempVy = p1[i].GetParams().Vy;

		p1[i].SetVelocity(p2[j].GetParams().Vx, p2[j].GetParams().Vy);
		p2[j].SetVelocity(tempVx, tempVy);

	}

		__syncthreads();
	
}

void Calc(Particle* h_a, Particle* h_b, Geometry g, int n) {

	// Allocate host memory
	// Initialize host array

	// Allocate arrays in Device memory
	Particle* d_a;
	Particle* d_b;

	size_t size= BLOCKS_NUMBER*THREADS_NUMBER*sizeof(Particle);

	cudaMalloc(&d_a, size);
	cudaMalloc(&d_b, size);

	// Copy memory from Host to Device
	cudaMemcpy(d_a, h_a, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, h_b, size, cudaMemcpyHostToDevice);

	// Block and Grid dimentions
	
	iter<<<BLOCKS_NUMBER,THREADS_NUMBER>>>(d_a, d_b, g, n);
	// Launch Kernel

	// Some kind of synchronization
	cudaDeviceSynchronize();

	cudaMemcpy(h_a, d_a, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_b, d_b, size, cudaMemcpyDeviceToHost);

	cudaFree(d_a);
	cudaFree(d_b);

}

void Calc2(Particle* h_a, Particle* h_b, int n) {
	// Allocate host memory
// Initialize host array

// Allocate arrays in Device memory
	Particle* d_a;
	Particle* d_b;

	size_t size = BLOCKS_NUMBER * THREADS_NUMBER * sizeof(Particle);

	cudaMalloc(&d_a, size);
	cudaMalloc(&d_b, size);

	// Copy memory from Host to Device
	cudaMemcpy(d_a, h_a, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, h_b, size, cudaMemcpyHostToDevice);

	// Block and Grid dimentions

	dim3 block(BLOCKS_NUMBER, BLOCKS_NUMBER);
	dim3 grid(THREADS_NUMBER, THREADS_NUMBER);
	iter2 << <grid, block >> > (d_a, d_b, n);
	// Launch Kernel

	// Some kind of synchronization
	cudaDeviceSynchronize();

	cudaMemcpy(h_a, d_a, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_b, d_b, size, cudaMemcpyDeviceToHost);

	cudaFree(d_a);
	cudaFree(d_b);
}