# Parallel N-Body Simulation with Barnes-Hut Approximation

The N-body problem is fundamental to many scientific fields and is concerned with the simulation of a system of particles under the influence of forces over time. Performing such simulations involving a large number of particles, typical in many real world problems, requires enormous computational time. This can however be significantly sped up by parallelizing the problem and using approximation algorithms to reduce the number of computations involved, making such simulations viable for practical purposes. 

This repository contains code that demonstrates the effectiveness of parallelization using [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) and the [Barnes-Hut algorithm](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation) as an approximation technique on the N-body problem. 

## Introduction

The N-body problem is fundamental to many scientific fields of study and is concerned with the simulation of a system of particles, under the influence of forces over time. The nature of these forces depends on the problem being solved with typical examples being gravitational forces experienced by a system of astronomical bodies, or the electrostatic Coulomb forces experienced by electrically charged particles such as protons and electrons in an atom, or atoms in a in molecule. Many real world phenomena are modeled and simulated as N-body problems, with the former being an example of its application in astrophysics for simulating planetary systems or stars in a galaxy, and the latter in molecular dynamics for problems such as protein folding. Other interesting applications include vortex flows in fluid dynamics and radiosity computation in computer graphics.

While solutions for 1-body problems are trivial, 2-body problems can be expressed as differential equations and solved using Kepler’s laws of motion. For problems involving three or more bodies, the solutions are either too complex or no known solution exists. As a result such systems are modeled and studied using computational methods and simulations. Simulations carried out on systems with large number of bodies over long periods of time involve a huge number of repetitive calculations. These simulations must also be completed within reasonable time frames for the results to be useful. This paper is concerned with the achievable speed ups on the N-body problem. It begins with a more formal definition of the N-body gravitational problem followed by an analysis of its sequential complexity. Various techniques for improving this complexity are discussed next, with the main focus being on the Barnes-Hut algorithm. This is followed by a discussion on parallelizing the N-body problem. The various experiments carried out are presented and the results are discussed in length.

## Defining the Problem

The N-body gravitational problem begins with a set of
N particles with initial positions P<sub>i</sub> = (x<sub>i</sub> , y<sub>i</sub> , z<sub>i</sub>) in three dimensional space, where i ∈ {1, 2, 3, ..., N }, moving with initial velocities V<sub>i</sub> = (v<sub>i</sub><sup>x</sup>, v<sub>i</sub><sup>y</sup>, v<sub>i</sub><sup>z</sup>). Each of these particles are under the influence of a gravitational force F<sub>i</sub> = (f<sub>i</sub><sup>x</sup>, f<sub>i</sub><sup>y</sup>, f<sub>i</sub><sup>z</sup>) due to their masses m i , according to Newton’s inverse-square law of gravitation. The component of the force in the x direction is given by:

<p style="text-align: center;">f<sub>i</sub><sup>x</sup> = G**Σ**<sub>j≠i</sub> **(**(m<sub>i</sub>m<sub>j</sub>) / (d<sup>2</sup>(i,j))**)** **(** (x<sub>i</sub>x<sub>j</sub>) / (d(i,j))**)**</p>

where G is the gravitational constant and d is the distance between particles i and j. Similar components of force apply in the y and z directions. The acceleration on a particle, A<sub>i</sub> = (a<sub>i</sub><sup>x</sup>, a<sub>i</sub><sup>y</sup>, a<sub>i</sub><sup>z</sup>), as a result of the force experienced in the x direction is given by:

<p style="text-align: center;">a<sub>i</sub><sup>x</sup> = f<sub>i</sub><sup>x</sup> / m<sub>i</sub></p>

This results in a new velocity v<sub>i</sub><sup>x'</sup> after time δt:

<p style="text-align: center;">v<sub>i</sub><sup>x'</sup> = v<sub>i</sub><sup>x</sup> + a<sub>i</sub><sup>x</sup> δt</p>

The N-body simulation then continues for a predetermined time period T , integrating with δt discrete time steps. In order to achieve accurate results, the time intervals between each integration step must be short.

## A Naive Algorithm

A straight forward algorithm for the N-body problem involves computing the forces experienced by each of the N bodies as a result of the gravitational attraction influenced by the other N-1 bodies, for each time step. The complexity of this algorithm for one iteration is Θ(N<sup>2</sup>). Further, since there are N bodies in the system, the spatial (memory) complexity of the algorithm is Θ(N). For practical applications however, where N is large, the quadratic time complexity of the algorithm quickly becomes intractable and better approaches are required.

```c
for (i ← 0 to N ) do
  for (j ← 0 to N ) do
    if (i == j) then
      continue;
    end if
    compute_force_on_i();
  end for
  compute_velocity_on_i();
  update_position_on_i();
end for
```
