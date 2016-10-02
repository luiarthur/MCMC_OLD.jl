using MCMC
using Base.Test

# Define immutable class State
immutable State
  mu::Float64
  sig2::Float64
end

@testset "Posterior Means" begin
  srand(1)
  using Distributions

  # Generate data
  const n = 1000
  const (mu,sig2) = (5,2)
  const y = rand(Normal(mu,sqrt(sig2)),n)

  # precomputes
  const ybar = mean(y)
  const (sig2_a, sig2_b) = (2.0, 1.0)
  rig(shape::Float64,rate::Float64) = 1 / rand(Gamma(shape,1/rate))


  # define an update function for State
  function update(state::State)
    const newMu = rand(Normal(ybar,sqrt(state.sig2/n)))
    const newSig2 = rig(sig2_a+n/2, sig2_b + sum((y-newMu).^2)/2)
    State(newMu,newSig2)
  end

  # run gibbs
  @time out = gibbs(State(0,1),update,10000,1000);

  ## post processing
  const postMu = map(o -> o.mu, out)
  const postSig2 = map(o -> o.sig2, out)
  @test abs(mean(postMu)-mu) < .05
  @test abs(mean(postSig2)-sig2) < .05
end

#= Personal tests:
# to test: vulture test.jl julia --color=yes test.jl

include("runtests.jl")
using RCall
R"library(rcommon) # install_github('luairthur/rcommon')"
plotPosts = R"plotPosts"
plotPosts([postMu postSig2])
=#
