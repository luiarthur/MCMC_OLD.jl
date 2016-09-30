module MCMC
export Specifications, gibbs

#=
using Lazy
function logfact(n)
 @bounce function loop(N, acc) 
   N == 0? acc : loop(N-1, log(N) + acc)
 end
 x = loop(n,0.0)
 println(x)
 x
end

=#

typealias State Dict{String,Any}

type Specifications
  init::Dict{String,Any}
  fcs::Dict{String,Any}
end

function gibbs(specs::Specifications,B::Int,burn::Int)
  out = [ State for i in 1:(B+burn) ]

  init = specs.init
  fcs = specs.fcs
  params = 

  function update(currState::State) 
    function loop
    end

  end

  out
end

end

#= Tests
f(x) = x + 1
typeof(f)

init = Dict(:mu => 1)
fcs = Dict(:mu => f)
specs = Specifications(init,fcs)
=#

