# Copyright 2025 by the authors.
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

using NamedArrays, BioSequences, BioSymbols, LinearAlgebra

export ones_weights, sum_all_wedges, sum_intern_wedges
export set_conductance, part_conductance

# TODO:
# wedges: make independent of numbers, use any symbol.

"""
    ones_weights(alphabet::Vector{<:BioSymbol}, t_size::Int)::Vector{<:Matrix}
    
Weight matrices with all values set to 1. Returns a list of size `tsize` 
where each entry is a matrix of dimension n × n and all values except 
for the diagonal are set to 1.

Functions like `setconductance` or `partconductance` require matrices
with transition values. If the weights of the conductance graph
should be all set to 1, this function can be used.

Arguments

- `alphabet`: Alphabet as a vector of symbols, e.g. A, T, C, G.
- `tsize`: Tuple size
"""
function ones_weights(alphabet::Vector{<:BioSymbol}, t_size::Int)::Vector{<:NamedArray}
    n = length(alphabet)
    W = ones(n, n) - I(n)
    W = NamedArray(W, (alphabet, alphabet))
    return [W for _ in 1:t_size]
end

"""
    sum_all_wedges(tuples::Vector{<:BioSequence}, W::Vector, n::Int)::Number

 Calculate sum of weighted edges for a set of tuples.

 This is a helper function used in setconductance and partconductance.

 Arguments
 - `tuples`: List of tuples represented by a vector of strings.
 - `W`: List of transition weight matrices. The size of the list
 must be the tuple size. Each list entry must have matrices of dimensions
 n × n (alphabet sizes).
 - `n::Int`: Alphabet size, e.g. |{A, T, C, G}| = 4.
"""
function sum_all_wedges(tuples::Vector{<:BioSequence}, W::Vector, n::Int)::Number
    l = length(W) # Tuple length
    sum([sum([W[i][tuple[i], j] for j in 1:n]) for i in 1:l for tuple in tuples])
end

"""
    sum_intern_wedges(tuples::Vector{<:BioSequence}, W::Vector{<:AbstractMatrix{<:Number}})::Number

Calculate sum of internal weighted edges for a set of tuples.

This is a helper function used in setconductance and partconductance.

Arguments
 - `tuples`: List of tuples represented by a vector of strings.
 - `W`: List of transition weight matrices. The size of the list
 must be the tuple size. Each list entry must have matrices of dimensions
 n × n (alphabet sizes).
"""
function sum_intern_wedges(tuples::Vector{<:BioSequence}, W::Vector{<:AbstractMatrix{<:Number}})::Number
    if length(tuples) <= 1 # one node has no internal edges
        return 0
    else # more than one node
        S = to_matrix(tuples)
        l = length(W) # Tuple size

        ro = eachrow(S) # get rows

        k = [(a, b) for a in ro for b in ro] # all combinations

        # idx contains list of tuple-pairs which differ in one letter:
        idx = [sum([a[i] != b[i] for i in 1:l]) == 1 for (a, b) in k]
        pairs = k[idx]
        if isempty(pairs) # cannot happen
            return 0 # nothing to do.
        end

        return sum([W[i][p1[i], p2[i]] for (p1, p2) in pairs for i in 1:l])
    end
end

"""
    set_conductance(tuples::Vector{<:BioSequence}, W::Vector{<:AbstractMatrix{<:Number}}, n::Int)::Number

 Calculate the conductance for a set of tuples.

 Arguments
 - `tuples`: List of tuples represented by a vector of `BioSymbols`.
 - `W`: List of transition weight matrices. The size of the list
 must be the tuple size. Each list entry must have matrices of dimensions
 n × n (alphabet sizes).
 - `n`: Alphabet size, e.g. |{A, T, C, G}| = 4.
"""
function set_conductance(tuples::Vector{<:BioSequence}, W::Vector{<:AbstractMatrix{<:Number}}, n::Int)::Number
    i = sum_intern_wedges(tuples, W)
    e = sum_all_wedges(tuples, W, n)
    return (e - i) / e
end

"""
    part_conductance(tuples::Vector{<:BioSequence}, p::Vector, W::Vector{<:AbstractMatrix{<:Number}}, n::Int)::Vector{<:Number}

Calculate the set conductance for a partition, i.e. a vector of set partitions.

Arguments
 - `tuples`: List of tuples represented by a vector of strings.
 - `p`: Partitions for tuples. The size of `p` must match the size of tuples.
 - `W`: List of transition weight matrices. The size of the list
 must be the tuple size. Each list entry must have matrices of dimensions
 n × n (alphabet sizes).
 - `n`: Alphabet size, e.g. |{A, T, C, G}| = 4.
"""
function part_conductance(tuples::Vector{<:BioSequence}, p::Vector, W::Vector{<:AbstractMatrix{<:Number}}, n::Int)::Vector{<:Number}
    P = [tuples[p.==i] for i in unique(p)]
    return [set_conductance(S, W, n) for S in P]
end

# to base package?
function to_matrix(S)
    l = [collect(s) for s in S]
    return hcat(l...)
end