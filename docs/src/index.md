```@meta
# Information for Documenter
CurrentModule = GCATConductance
```

```@contents
Pages = ["index.md"]
```

# Introduction

## Definitions

Let $B$ be an alphabet and $n$ the alphabet size. Let $l$ denote the lengths of words. Let $S\subset B^l$ a set of tuples. For instance, $\mathcal{B} = \{A, T, C, G\}$ denotes the bases of the genetic code.

## Weight matrix

A symmetric $n \times n$-matrix where each row and each column represents a (unique) letter from the alphabet $B$ and the diagonal elements are 0 is called a transition matrix.
$W_{a,b}$ refers to the element in $W$ in row $a$ and column $b$.

Let us create some results. The following data frame contains four codons. 
Column `t` lists the true codons whereas `e` lists the estimated codons.

```@example rt
using GCATConductance, BioSequences, BioSymbols

W = ones_weights([DNA_A, DNA_T], 2)
W[1]
```
```@example rt
W[2]
```

## Conductance and robustness (`set_conductance`)

The summed weights of all edges in a set $S$ is defined as:

``
E(S) =
\underset{\text{ every tuple }}{\underbrace{\sum_{t\in S}}}
\underset{\text{ every position }}{\underbrace{\sum_{i=1}^l}}
\underset{\text{ other letters }}{\underbrace{\sum_{b \in B} W_{t_i, b}}}
``

```@example rt
es = sum_all_wedges([dna"AA", dna"AT"], W, 2)
```

The summed weights of all __internal__ edges in a set $S$ is defined as:

``I(S) = \sum_{\text{all pairs with edges } (t,u) \in S \times S} \sum_{i=1}^l W_{t_i, u_i}``

```@example rt
is = sum_intern_wedges([dna"AA", dna"AT"], W)
```

The (set) conductance is ratio of the number of outgoing edges to all edges of ``S``.

``\varphi(S)=\frac{E(S)-I(S)}{E(S)}``

```@example rt
sc = set_conductance([dna"AA", dna"AT"], W, 2)
```

See [`GCATConductance.set_conductance`](@ref) for details.

The (set) robustness is ratio of the number of internal edges to all edges of ``S``.

``\rho(S)=\frac{I(S)}{E(S)} = 1 - \varphi(S)``

## Conductance for a partition (`part_conductance`)

```@example rt
p = [1, 1, 1, 2] # a partition for the tuples below
pc = part_conductance([dna"AA", dna"AT", dna"TA", dna"TT"], p, W, 2)
```

# API

```@autodocs
Order   = [:module, :type, :function]
Modules = [GCATConductance]
```
