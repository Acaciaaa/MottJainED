#!/usr/bin/env julia
import Pkg
const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT; io=devnull)
include(joinpath(PROJECT_ROOT, "experimental", "FastED.jl"))
include(joinpath(PROJECT_ROOT, "experimental", "FastMuSearch.jl"))

isempty(ARGS) && error("Usage: fast_mu_search.jl plan|run --config=PROFILE [--threads-per-worker=8]")
command = first(ARGS)
options = Dict(split(arg[3:end], "="; limit=2) for arg in ARGS[2:end] if startswith(arg, "--"))
haskey(options, "config") || error("Missing --config=PROFILE")
spec = FastED.load_spec(FastED.project_path(options["config"]))
search = FastMuSearch.options(spec)
if command == "plan"
    println(FastED.plan(spec))
    println(search)
    println("One Hamiltonian, four resident sector workers, no cache release or next-point submission.")
elseif command == "run"
    FastMuSearch.run(spec; threads_per_worker=parse(Int, get(options, "threads-per-worker", "8")))
else
    error("Unknown command: $command")
end
