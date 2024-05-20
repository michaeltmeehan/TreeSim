module TreeSim

using DataFrames
using DataFramesMeta


mutable struct Host
	id::Int64
    type::Int64
    infector::Union{Int64, Nothing}
    infector_type::Union{Int64, Nothing}
	root::Union{Float64, Nothing}
	leaf_times::Vector{Float64}
	leaf_ids::Vector{Int64}
end


include("utils.jl")
include("tree_operations.jl")
include("probabilities.jl")
include("constraints.jl")
include("backward.jl")
include("forward.jl")
include("sampling.jl")
include("tree_construction.jl")

export simulate


end # module TreeSim
