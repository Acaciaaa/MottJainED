#!/usr/bin/env julia
import Pkg
const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT; io=devnull)
include(joinpath(PROJECT_ROOT, "experimental", "FastED.jl"))
include(joinpath(PROJECT_ROOT, "experimental", "FastMuSearch.jl"))

isempty(ARGS) && error("Usage: fast_mu_search.jl plan|run --config=PROFILE")
command = first(ARGS)
options = Dict(split(arg[3:end], "="; limit=2) for arg in ARGS[2:end] if startswith(arg, "--"))
haskey(options, "config") || error("Missing --config=PROFILE")
all(key == "config" for key in keys(options)) || error("Only --config=PROFILE is supported; worker processes have been removed")
spec = FastED.load_spec(FastED.project_path(options["config"]))
search = FastMuSearch.options(spec)
if command == "plan"
    println(FastED.plan(spec))
    println(search)
    println("One Hamiltonian; sectors run sequentially in one process, using disk cache and complete CSVs.")
    println("Slurm requests 8 CPUs; no cache release or next-point submission.")
elseif command == "run"
    FastMuSearch.run(spec)
else
    error("Unknown command: $command")
end
