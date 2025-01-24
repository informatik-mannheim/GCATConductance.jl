# some tests

using GCATConductance, BioSequences, BioSymbols, Test

@testset "one_weights" begin
    w = ones_weights([DNA_A, DNA_T], 2)
    sc = set_conductance([dna"AA", dna"AT"], w, 2)
    @test sc == 0.5

    w = ones_weights([DNA_A, DNA_T, DNA_C, DNA_G], 3)
    sc = set_conductance([dna"ATG"], w, 4)
    @test sc == 1
end

@testset "all weights" begin
    w = ones_weights([DNA_A, DNA_T], 2)
    s = sum_all_wedges([dna"AT"], w, 2)
    @test s == 1 * (1 + 1)

    w = ones_weights([DNA_A, DNA_T, DNA_C, DNA_G], 3)
    s = sum_all_wedges([dna"ATG"], w, 4)
    @test s == 1 * (3 + 3 + 3)

    s = sum_all_wedges([dna"ATG", dna"ATT"], w, 4)
    @test s == 2 * (3 + 3 + 3)
end

@testset "set conductance" begin
    w = ones_weights([DNA_A, DNA_T], 2)
    sc = set_conductance([dna"AA", dna"AT"], w, 2)
    @test sc == 2 / 4

    w = ones_weights([DNA_A, DNA_T, DNA_C, DNA_G], 3)
    sc = set_conductance([dna"ATG"], w, 4)
    @test sc == 1

    w = ones_weights([DNA_A, DNA_T, DNA_C, DNA_G], 3)
    sc = set_conductance([dna"ATG", dna"ATT"], w, 4)
    @test sc == (2 * 9 - 2) / (2 * 9)
end