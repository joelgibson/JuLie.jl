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
end
