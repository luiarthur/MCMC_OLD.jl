__precompile__()
"""
My MCMC stuff

Currently the only method is gibbs. Eventually, I would
like to add a Metropolis Hastings method.
"""
module MCMC

  export gibbs, GibbsState
  """
  Runs gibbs sampler with specifications (spec), samples (B), 
  and burn-in (burn). 

  # Example:

  ```julia
  using Distributions

  # Generate data
  const n = 1000
  const (mu,sig2) = (5,2)
  const y = rand(Normal(mu,sqrt(sig2)),n)

  # precomputes
  const ybar = mean(y)
  const (sig2_a, sig2_b) = (2.0, 1.0)
  rig(shape::Float64,rate::Float64) = 1 / rand(Gamma(shape,1/rate))

  # Define State
  immutable State
    mu::Float64
    sig2::Float64
  end

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

  ```
  """

  function gibbs{T}(init::T, update, B::Int, burn::Int)
    const out = Array{T,1}( (B+burn) )
    out[1] = init
    for i in 2:(B+burn)
      out[i] = update(out[i-1])
    end
    out[ (burn+1):end ]
  end

end # module
