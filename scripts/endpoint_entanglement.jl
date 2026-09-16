#!/usr/bin/env julia

println("Endpoint entanglement starting: ", join(ARGS, " "))
flush(stdout)

import Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT; io=devnull)

include(joinpath(PROJECT_ROOT, "experimental", "EndpointEntanglement.jl"))
using .EndpointEntanglement

function parse_options(args)
    options = Dict{String,String}()
    for argument in args
        startswith(argument, "--") || throw(ArgumentError(
            "unknown positional argument '$argument'; use --key=value options",
        ))
        parts = split(argument[3:end], "="; limit=2)
        options[parts[1]] = length(parts) == 2 ? parts[2] : "true"
    end
    return options
end

option_bool(options, key, default=false) = haskey(options, key) ?
    lowercase(options[key]) in ("1", "true", "yes", "on") : default

function usage()
    println("""
Usage:
  julia --project=. scripts/endpoint_entanglement.jl plan --config=PROFILE.toml
  julia --project=. scripts/endpoint_entanglement.jl run --config=PROFILE.toml [--force]

The run command solves the two endpoint ground states, then computes only the
configured fixed-QA/F3A/F8A edge sectors.  Each dense SVD is checkpointed.
""")
end

isempty(ARGS) && (usage(); exit(1))
command = lowercase(ARGS[1])
options = parse_options(ARGS[2:end])
haskey(options, "config") || throw(ArgumentError("missing --config=..."))
config_path = isabspath(options["config"]) ? options["config"] :
    normpath(joinpath(PROJECT_ROOT, options["config"]))
spec = load_spec(config_path)

if command == "plan"
    println(plan(spec))
elseif command == "run"
    println(plan(spec))
    flush(stdout)
    result = run_endpoint_entanglement(spec; force=option_bool(options, "force"))
    println("completed: ", result.output)
else
    usage()
    throw(ArgumentError("unknown command '$command'"))
end
