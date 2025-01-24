# some tests

using GCATConductance, BioSequences, BioSymbols, Test

@testset "Conductance" begin
    w = ones_weights([DNA_A, DNA_T], 2)
    sc = set_conductance([dna"AA", dna"AT"], w, 2)
    @test sc == 0.5
    
    w = ones_weights([DNA_A, DNA_T, DNA_C, DNA_G], 3)
    sc = set_conductance([dna"ATG"], w, 4)
    @test sc == 1
end