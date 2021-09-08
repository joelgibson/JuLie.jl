using Test
using LinearAlgebra

function count_elements(grp::CoxGrp)
	frontier = Set([one(grp)])
	count = length(frontier)
	while length(frontier) > 0
		frontier = Set(Base.Iterators.flatten((mult_gen(x, s) for s=setdiff(BitSet(1:rank(grp)), right_descent_set(x))) for x=frontier))
		count += length(frontier)
	end
	return count
end

function A(n::Int)
	mat = 2 * Matrix{Int}(LinearAlgebra.I, n, n)
	for i in 1:n-1
		mat[i, i+1] = -1
		mat[i+1, i] = -1
	end
	return mat
end

@testset "Coxeter Groups: Matrix representation" begin
    A2 = coxeter_group_mat([[2, -1] [-1, 2]])
    id = one(A2)
    s, t = gens(A2)

    @test s == s
    @test s != t
    @test t == t

    @test s*s == id
    @test t*t == id

    @test [length(x) for x=[id, s, t, s*t, t*s, s*t*s]] == [0, 1, 1, 2, 2, 3]
    @test inv(s*t*s) == s*t*s

    # Exponentiation
    @test (s*t)^0 == id
    @test (s*t)^1 == s*t
    @test (s*t)^2 == s*t*s*t
    @test (s*t)^-1 == t*s

    @test lexword(s*t*s) == [1, 2, 1]

    @test count_elements(A2) == 6

    # Non-integral Cartan matrices
    H2 = coxeter_group_mat([[1, 5] [5, 1]])
    @test count_elements(H2) == 10

    H3 = coxeter_group_mat([[1, 5, 2] [5, 1, 3] [2, 3, 1]])
    @test count_elements(H3) == 120

    H4 = coxeter_group_mat([
        1 5 2 2
        5 1 3 2
        2 3 1 3
        2 2 3 1
    ])
    #@test count_elements(H4) == 14400

    A7 = coxeter_group_mat(A(7))
    #@test count_elements(A7) == 40320
end
