# AlignedGP

AlignedGP uses interval arithmetic to backpropagate the hits metric through a tree, finding optimal additive constants at every node in O(n log n) time. Lexicase selection over the per-point hits vector drives the search.

The key mechanism: given a target interval per data point, interval inversion maps root-level targets backward through each tree node to produce per-point surrogate intervals for child subtrees. A constant-stabbing routine then finds the additive constant that hits the maximum number of surrogate intervals — exactly and efficiently, without gradient search.

For the theoretical background, see [the technical report](doc/techrep.pdf).

## Installation

Requires Julia. Clone the repo and instantiate the included manifest:

```julia
] instantiate
```

## Quickstart

```julia
include("scripts/run_gp.jl")
```

This runs AlignedGP on the Keijzer-4 benchmark (`x³ eˣ sin x cos x`) and prints a progress report every two seconds showing hits, complexity, and effort.
