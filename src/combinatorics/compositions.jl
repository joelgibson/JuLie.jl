################################################################################
# Compositions.
#
# Copyright (C) 2021 Ulrich Thiel, ulthiel.com/math
################################################################################

export Composition, num_compositions, compositions


"""
    Composition{T} <: AbstractArray{T,1}

A **composition** of an integer n ≥ 0 is a sequence (λ₁,…,λᵣ) of positive integers whose sum is equal to n.

# Examples
```julia-repl
julia> c=Composition([2,1]);
julia> sum(c)
3
julia> c=Composition(Int8[2,1]); #You can use smaller integer types
```

# References
1. Wikipedia, [Composition (combinatorics)](https://en.wikipedia.org/wiki/Composition_(combinatorics))
"""
struct Composition{T} <: AbstractArray{T,1}
    p::Array{T,1}
end

# The following are functions to make the Composition struct array-like.
function Base.show(io::IO, ::MIME"text/plain", P::Composition)
    print(io, P.p)
end

function Base.size(P::Composition)
    return size(P.p)
end

function Base.length(P::Composition)
    return length(P.p)
end

function Base.getindex(P::Composition, i::Int)
    return getindex(P.p,i)
end

function Base.setindex!(P::Composition, x::Integer, i::Int)
    return setindex!(P.p,x,i)
end

function Base.copy(P::Composition{T}) where T<:Integer
  return Composition{T}(copy(P.p))
end


"""
    num_compositions(n::Integer)

The number of compositons of an integer ``n>0`` is equal to ``2^{n-1}``. For ``n=0`` there is 1 composition (the empty one) by convention.

# References
1. The On-Line Encyclopedia of Integer Sequences, [A011782](https://oeis.org/A011782)
"""
function num_compositions(n::Integer)
    if n==0
        return ZZ(1)
    else
        return ZZ(2)^(n-1)
    end
end


"""
    compositions(n::Integer, k::Integer)

Returns an array of all compositions of n into k parts. The algorithm used is Algorithm 72 "Composition Generator" by Hellerman & Ogden (1961), which also refers to Chapter 6 of Riordan (1958). De-gotoed by E. Thiel. I don't know if there are faster algorithms but this one is already very fast.

# References
1. Hellerman, L. & Ogden, S. (1961). Algorithm 72: composition generator. *Communications of the ACM, 4*(11), 498. [https://doi.org/10.1145/366813.366837](https://doi.org/10.1145/366813.366837)
2. Riordan, J. (1958). *An introduction to Combinatorial Analysis*. John Wiley & Sons, Inc., New York.
"""
function compositions(n::Integer, k::Integer)

    # Argument checking
    n >= 0 || throw(ArgumentError("n ≥ 0 required"))
    k >= 0 || throw(ArgumentError("k ≥ 0 required"))

    # Use type of n
    T = typeof(n)

    # This will be the array of compositions
    C = Composition{T}[]

    # Special cases
    if k > n
        return C
    elseif n == 0
        c = Composition{T}([])
        push!(C,copy(c))
        return C
    elseif k == 0
        return C
    end

    # Initialize c
    c = Composition{T}([ [1 for i=1:k-1] ; n-k+1 ])
    push!(C,copy(c))

    if k==1
        return C
    end

    # Initialize d
    d = Array{T,1}([0 for i=1:k])

    while true

        for jk=1:k
            d[jk] = c[jk] -1
        end

        j = k
        Lfound = false

        while j>0 && Lfound == false

            if d[j]>0 #label test successful, then go to set
                d[j]=0
                d[j-1]=d[j-1] + 1
                d[k]=c[j]-2
                for jk=1:k
                    c[jk] = d[jk] + 1
                end
                push!(C,copy(c))
                Lfound=true #while-loop can be immediately left here

            else

                j=j-1 #decrease while loop counter and ...

                # ... now check whether we have already reached the first index.
                if j == 1 # if so, return
                    return C
                end

            end

        end

    end

end


"""
    compositions(n::Integer)

Returns an array of all compositions of an integer n. This iterates over compositions of n into k parts for 1 ≤ k ≤ n.
"""
function compositions(n::Integer)

    # Argument checking
    n >= 0 || throw(ArgumentError("n ≥ 0 required"))

    # Use type of n
    T = typeof(n)

    # This will be the array of compositions
    C = Composition{T}[]

    # Special case
    if n == 0
        c = Composition{T}([])
        push!(C,copy(c))
        return C
    end

    for k=1:n
        append!(C,compositions(n,k))
    end

    return C

end
