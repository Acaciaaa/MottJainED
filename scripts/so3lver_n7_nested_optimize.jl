#!/usr/bin/env julia

# Keep the production N=6 driver byte-for-byte stable while reusing the exact
# same search implementation for the requested N=7 finite-size comparison.
# The core driver's historical override is injected only after this wrapper
# has verified that the selected profile is exactly N=7.

import TOML

const WRAPPER_PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))

function wrapper_config_path(args)
    for argument in args
        startswith(argument, "--config=") || continue
        return abspath(split(argument, "="; limit=2)[2])
    end
    return joinpath(
        WRAPPER_PROJECT_ROOT,
        "config", "so3lver", "n7_nested_six_term_search.toml",
    )
end

any(startswith(argument, "--allow-test-size") for argument in ARGS) && error(
    "the N=7 production wrapper manages the core size override internally",
)
config_path = wrapper_config_path(ARGS)
isfile(config_path) || error("configuration not found: $config_path")
config = TOML.parsefile(config_path)
Int(config["run"]["nm"]) == 7 || error(
    "the N=7 production wrapper only accepts run.nm = 7",
)

push!(ARGS, "--allow-test-size=true")
include(joinpath(@__DIR__, "so3lver_n6_nested_optimize.jl"))
