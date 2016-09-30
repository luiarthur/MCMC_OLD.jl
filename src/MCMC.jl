module MCMC

  export Specifications, gibbs, State

  typealias State Dict{String,Any}
  typealias FullConditionals Dict{String,Any}

  immutable Specifications
    init::State
    fcs::FullConditionals
  end

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
