#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

#include<iostream>
#include <cuda.h>

#include <curand.h>
#include <curand_kernel.h>

#include "ParticleSystem.h"
#include "Particle.h"

#define THREADS_NUMBER 2
#define BLOCKS_NUMBER 1

//void Particle::Calc(Particle* h_a, Geometry g, int n);
