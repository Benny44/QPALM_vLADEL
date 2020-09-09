module QPALM

const depsfile = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
if isfile(depsfile)
    include(depsfile)
else
    error("QPALM not properly installed. Please run Pkg.build(\"QPALM\").")
end

const Maybe{T} = Union{T, Nothing} where T

include("const.jl")
include("types.jl")
include("wrappers.jl")

end # module
