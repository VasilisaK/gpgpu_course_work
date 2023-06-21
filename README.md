# Particle Game
Repo for gpgpu course work "Particle game"

## Authors
Команда курсовой работы:
1. Кудрявцева Василиса (Vasilisa Kudryavtseva)
2. Пестряков Данил (Danil Pestryakov)
3. Суриков Илья (Ilya Surikov)

## Rules

The movement of 10000 particles is calculated on the GPU during the game.

Two sources, each source generates 5000 particles at the beginning of the game at the same time.
Each source generates its own type of particles.
There are two types of particles (blue and green).
Particles of the same type do not interact with each other.
Particles of different types repel absolutely elastically.

Particles spawn as the lifetime of each of them expires, have a random direction, speed, and lifetime of their certain limits.
The particles are repelled from the walls of geometric objects.

There is a particle counter. The game ends when a certain number of particles, set by the user, falls into the basket and dissappear in it.

## Using
* Visual Studio 2019 (C++)
* CUDA 10.2

