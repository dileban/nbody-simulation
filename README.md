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

## Parallelizing the N-Body Problem

Parallelizing the N-body simulation on a cluster can give significant speedups, depending on the number of processors and interconnection network used. An intuitive approach to parallelization is to partition the problem into P groups, where P is the number of processors, each processor being responsible for N/P particles. Each processor is also responsible for computing the forces experienced by the particles in its group and updating their positions. At the end of each iteration, the processors exchange the new positions of their particles with other processors, so that all processors now have complete updated view of the system. While parallel execution can help speedup the simulation, this does not improve the complexity of the algorithm. For large problem sizes this can still be a problem.

### Improving time complexity

The algorithm described so far is a simple particle-particle approach where the mutual forces between every pair of particles is the system is computed. A number of [other](http://www.amara.com/papers/nbody.html) [interesting](http://www.scholarpedia.org/article/N-body) simulation algorithms for N-body simulations exists that can give an improved time complexity. Some of the popular methods use the divide and conquer approach where the system is spatially decomposed into a hierarchical cluster in the form of a tree, with each node describing the mass and cubic volume of a cluster in the system. These algorithms are also referred to as hierarchical or tree (code) methods. Two widely used hierarchical algorithms in real world applications are the [Barnes-Hut](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation) algorithm and [Greengard’s Fast Multipole Method](https://en.wikipedia.org/wiki/Fast_multipole_method).

Both these algorithms give significant improvements, with Barnes-Hut resulting in O(Nlog<sub>2</sub>N) time complexity and the Fast Multipole Method (FMM) in O(N) complexity for uniform distributions. The efficient implementation of the FMM however is significantly more complex for a three dimensional system and in general requires a larger number of iterations before the algorithm can show comparatively better [speedups](http://www.cs.cmu.edu/~scandal/papers/dimacs-nbody.html). The code presented here is primarily concerned with achievable speedups through
parallelization, and will thus focus on the Barnes-Hut algorithm for its simplicity, effectiveness in reducing the overall running time and its wider usage in real world applications.

### The Barnes-Hut algorithm

The use of the Barnes-Hut algorithm for force calculation in the N-body problem proceeds in two phases. The first phase involves constructing a hierarchical octtree representing the three dimensional space through recursive subdivision of the root cell (containing all particles in
the system) into eight cubic subcells of equal size, and then progressing recursively on each of the subcells until all subcells contain at most one particle. The leaves in the tree represent cells containing one particle, while the parent nodes contain the total mass and center of mass of the subtree below it.

The second phase involves computing the forces experienced by each of the particles by traversing the tree. This phase makes use of the fact that the interaction between a particle and a distant cluster of particles can be approximated to a single combined force between the particle and the cluster as a whole. The calculation is then reduced to using the total mass and center of mass of the cluster.

The algorithm for force calculation on each particle proceeds by starting at the root cell. If the distance between the particle and the center of mass of the cluster is D, and the length of the cell is l, the interaction between the particle and the cell (cluster) as a whole is calculated if l/D < θ, where θ is a constant accuracy parameter, typically ≤ 1.0. If not, the current cell is resolved into its eight subcells, and each of them is examined recursively.

Once the computation of forces experienced by all particles is complete, the N-body algorithm proceeds to compute the new velocities and positions, discarding the tree. A new tree is then regenerated for the next iteration and all subsequent iterations. In the parallel N-body algorithm described in this paper, where each processor is responsible for N/P particles, the tree for the entire system is generated on each processor. The interactions however are only computed for the particles belonging to each of the processors. The construction of the tree requires time O(Nlog<sub>2</sub>N). Similar time is taken for computing all the forces using the tree, depending on the value of θ. The overall running time complexity of using the Barnes-Hut algorithm is thus O(Nlog<sub>2</sub>N).

### Strategies for optimal solutions

This project examines the use of the popular Barnes-Hut algorithm as a feasible alternative to the particle-particle method for large values of N , by generating the octtree for the entire system on each of the processors in parallel, for every iteration. One of the problems with this approach is that an uneven distribution of particles in the system will result in an unbalanced tree, and thus an imbalance in the load distribution with some of the processors performing more calculations that others. Depending on the particle distribution this can have a significant impact on the overall running time. However, it is possible to further optimize the parallel computation time by employing a number of other strategies.

A well known approach described by [J. Salmon](https://thesis.library.caltech.edu/6291/) uses the method of orthogonal recursive bisection (ORB) to recursively decompose the system into sub-volumes similar to the approach described above. However the strategy for decomposition is such that an equal number of particles lie on either side of the division. The decomposition into sub-volumes continues for as many levels as required, producing a well balanced tree. Each of the processors are then assigned one or more sub-volumes, balancing the computational load on the cluster.

A further optimization can be made by using the fact that particles move very small distances in every iteration (for small values of δt). Thus by not discarding the entire tree in every iteration, except regenerating those parts
of the tree whose corresponding particles have moved outside the cell, a significant computation time can be saved when N is large.

A third optimization involves parallelizing the tree generation such that each of the node generates a partial tree contain particles belonging to it. Once generated, the nodes in the cluster exchange their partial trees with each other, resulting in each of the nodes having a complete tree representing the system.
