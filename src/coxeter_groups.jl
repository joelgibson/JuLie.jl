export CoxGrp, CoxElt, lexword, rank, gens, coxeter_group_mat, right_descent_set, left_descent_set, parent, mult_gen

using LinearAlgebra: I

"Abstract supertype for Coxeter groups"
abstract type CoxGrp end

"Abstract supertype for Coxeter group elements"
abstract type CoxElt end

"The length of a Coxeter group element is the length of any reduced expression in the generators."
function Base.length(x::CoxElt)
    count = 0
    while true
        desc = right_descent_set(x)
        if length(desc) == 0
            break
        end
        x = mult_gen(x, first(desc))
        count += 1
    end
    return count
end

"Returns `true` if `x` is the identity element in the group."
function Base.isone(x::CoxElt)
    return length(right_descent_set(x)) == 0
end

"The lexicographically least reduced expression for `x`"
function lexword(x::CoxElt)
	word = Array{Int, 1}()
	while true
        desc = left_descent_set(x)
        if length(desc) == 0
            break
        end
        s = first(desc)
		x = mult_gen(x, s)
		push!(word, s)
	end
	return word
end

"The product ``xy`` of Coxeter group elements."
function Base.:(*)(x::CoxElt, y::CoxElt)
    while true
        desc = left_descent_set(y)
        if length(desc) == 0
            break
        end
        s = first(desc)
        x = mult_gen(x, s)
        y = mult_gen(s, y)
    end
    return x
end

"The inverse ``x^-1`` of a Coxeter group element."
function Base.inv(x::CoxElt)
	xinv = one(parent(x))
	while true
		desc = right_descent_set(x)
		if length(desc) == 0
			break
		end
		s = first(desc)
		x = mult_gen(x, s)
		xinv = mult_gen(xinv, s)
	end
	return xinv
end

"The power ``x^n`` of a Coxeter group element."
function Base.:(^)(x::CoxElt, n::Integer)
	# Deal with negative powers
	if n < 0
		x = inv(x)
		n = -n
	end

	# Exponentiate by squaring
	acc = one(parent(x))
	pow2 = x
	while n > 0
		if n % 2 == 1
			acc *= pow2
		end
		pow2 *= pow2
		n ÷= 2
	end

	return acc
end




"""
A Coxeter group represented by matrices with entries of type `T`. Each matrix is the action of a group
element in the basis of simple roots.
"""
mutable struct CoxGrpMat{T} <: CoxGrp
	cartanMat::Array{T, 2}
	simpRefls::Array{Array{T, 2}, 1}
end

rank(grp::CoxGrpMat) = size(grp.cartanMat)[1]
gens(grp::CoxGrpMat) = [CoxEltMat(grp, refl, refl) for refl=grp.simpRefls]
function Base.one(grp::CoxGrpMat{T}) where {T}
	eye = Matrix{T}(I, rank(grp), rank(grp))
	return CoxEltMat(grp, eye, eye)
end

raw"""
Create a `CoxGrpMat`, a Coxeter group where elements are represented by matrices acting on the
basis of simple roots. The argument should be a Cartan matrix."""
function coxeter_group_mat(cartanMat::Matrix{T}) where {T}
    # Should check that the matrix is a Cartan matrix here.
	n, m = size(cartanMat)
	@assert n == m

	# The simple reflection with index i acts on the weight space as
	#    s_i(x) = x - <x, alpha_i^> alpha_i
	# When x is written in the root basis, then <x, alpha_i^vee> is the
	# inner product of the entries of x with the ith row of the Cartan matrix.
	# In particular this takes the jth coordinate vector to e_j - a_ij e_i.
	simpRefls = [Matrix{T}(I, n, n) for k=1:n]
	for i in 1:n
		for j in 1:n
			simpRefls[i][i, j] -= cartanMat[i, j]
		end
	end

	return CoxGrpMat(cartanMat, simpRefls)
end

struct CoxEltMat{T} <: CoxElt
	grp::CoxGrpMat{T}
	mat::Array{T, 2}
	inv::Array{T, 2}
end

parent(x::CoxEltMat) = x.grp
Base.isequal(x::CoxEltMat, y::CoxEltMat) = x.grp == y.grp && x.mat == y.mat
Base.hash(x::CoxEltMat) = hash(x.grp, hash(x.mat))
Base.:(==)(x::CoxEltMat, y::CoxEltMat) = Base.isequal(x, y)
mult_gen(s::Integer, x::CoxEltMat) = CoxEltMat(x.grp, x.grp.simpRefls[s] * x.mat, x.inv * x.grp.simpRefls[s])
mult_gen(x::CoxEltMat, s::Integer) = CoxEltMat(x.grp, x.mat * x.grp.simpRefls[s], x.grp.simpRefls[s] * x.inv)

raw"""
`x.mat` represents the action of `x` on the simple root basis, hence the columns of the matrix are
``[w(\alpha_1) \cdots w(\alpha_r)]``. We have that ``ws < w`` if and only if ``w(\alpha_s) < 0``, and
negative roots are simple to distinguish in the simple root basis: they have at least one negative
coordinate (in fact all their coordinates will be negative or zero). Therefore the right descent set
is the set of columns which have at least one negative entry.
"""
function right_descent_set(x::CoxEltMat)
	return BitSet(s for s=1:rank(x.grp) if any(x.mat[t, s] < 0 for t=1:rank(x.grp)))
end

raw"""
The left descents of ``w`` are the right descents of ``w^{-1}``.
"""
function left_descent_set(x::CoxEltMat)
	return BitSet(s for s=1:rank(x.grp) if any(x.inv[t, s] < 0 for t=1:rank(x.grp)))
end
