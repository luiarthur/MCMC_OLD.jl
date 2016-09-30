using Base.Test
@testset "Posterior Means" begin
  srand(1)
  include("../src/mcmc.jl")
  using Distributions, MCMC

  const n = 1000
  const μ = 5
  const σ² = 2
  const y = randn(n) * sqrt(σ²) + μ
  const ybar = mean(y)
  const a = 2
  const b = 1

  rig(shape::Float64,rate::Float64) = 1 / rand(Gamma(shape,1/rate))

  function updateσ²(param)
    mu = param["mu"]
    ss = sum((y - mu).^2)
    rig(a+n/2, b+ss/2)
  end

  function updateμ(param)
    s2 = param["sig2"]
    randn() * sqrt(s2/n) + ybar
  end

  const specs = Specifications(
    # init:
    Dict(
      "mu" => 0,
      "sig2" => 1
    ),
    # full conditionals:
    Dict(
      "mu" => updateμ,
      "sig2" => updateσ²
    )
  )

  B = 10000
  @time out = gibbs(specs,B,1000)

  postMu = [out[i]["mu"] for i in 1:B]
  postSig2 = [out[i]["sig2"] for i in 1:B]

  @test true
  @test abs(mean(postMu) - μ) < .05
  @test abs(mean(postSig2) - σ²) < .05
end

#= Personal tests:
# to test: vulture test.jl julia --color=yes test.jl

include("runtests.jl")
using RCall
R"library(rcommon) # install_github('luairthur/rcommon')"
plotPosts = R"plotPosts"
plotPosts([postMu postSig2])
=#
