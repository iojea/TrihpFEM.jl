# TrihpFEM

[![Build Status](https://github.com/iojea/TrihpFEM.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/iojea/TrihpFEM.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://iojea.github.io/TrihpFEM.jl/dev)
<!-- [![Documentation](https://github.com/iojea/TrihpFEM.jl/actions/workflows/Docs.yml/badge.svg)](https://iojea.github.io/TrihpFEM.jl/) -->
[![Coverage](https://codecov.io/gh/iojea/TrihpFEM.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/iojea/TrihpFEM.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)


## To Do

- [] Review all the assembly process and reorder things when needed. Some examples:
  - [] Meshes are based on Dictionaries. Dictionaries iterate over _values_, but I often need to iterate over _keys_ or _pairs_. Options: change the structure of the Dictionaries, or define a speciallized method for iterating.
  - [] Some parts of code are probably misplaced. E.g.: PolyFields is loaded before Meshes because of CMatrix. However, AffineToRef is overloaded in Assembly. This deserves a review.
  - [] Measures are meant to overshadow Meshes, so that Assembly speaks with Measures instead of Meshes, but this is not completely done, yet.
  - [] Some getter functions are named with the `get` prefix, some are not. Some inner functions are named with the prefix `_`, some are not. This should be standarized.
  - [] Some names of structs or functions are not very descriptive. This should be revised.

- [] FESolution should be expanded and `jump` introduced. This would also imply changes in Measures (maybe introducing SkeletonMeasure or something like that).
- [] Problems should be expanded to h and hp problems.