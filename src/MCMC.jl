"""
My MCMC stuff

Currently the only method is gibbs. Eventually, I would
like to add a Metropolis Hastings method.
"""
module MCMC

  export Specifications, gibbs, State

  typealias State Dict{String,Any}
  typealias FullConditionals Dict{String,Any}

  immutable Specifications
    init::State
    fcs::FullConditionals
  end

  """
  Runs gibbs sampler with specifications (spec), samplers (B), 
  and burn-in (burn). 


  # Example:

  ```julia
  using Distributions

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

  ```
  """
  function gibbs(specs::Specifications,B::Int,burn::Int)
    const out = Array{State,1}( (B+burn) )

    const init = specs.init
    const fcs = specs.fcs
    const params = keys(init)
    const d = length(params)

    function update(currState::State)
      newState = Dict( collect(currState) ) # important step!
      for p in params
        newState[p] = fcs[p](newState)
      end
      newState
    end

    out[1] = init
    for i in 2:(B+burn)
      out[i] = update(out[i-1])
    end

    out[(burn+1):end]
  end


end # module
