#!/usr/bin/env julia

import Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT; io=devnull)

include(joinpath(PROJECT_ROOT, "experimental", "FastED.jl"))
using .FastED

function parse_options(args)
    options = Dict{String,String}()
    for argument in args
        startswith(argument, "--") || throw(ArgumentError(
            "Unknown positional argument '$argument'; use --key=value options",
        ))
        parts = split(argument[3:end], "="; limit=2)
        options[parts[1]] = length(parts) == 2 ? parts[2] : "true"
    end
    return options
end

option_bool(options, key, default=false) = haskey(options, key) ?
    lowercase(options[key]) in ("1", "true", "yes", "on") : default

function required_int(options, key)
    haskey(options, key) || throw(ArgumentError("Missing required option --$key=..."))
    return parse(Int, options[key])
end

function usage()
    println("""
Usage:
  julia --project=. scripts/fast_ed.jl plan --config=PROFILE.toml
  julia --project=. scripts/fast_ed.jl prepare --config=PROFILE.toml [--force]
  julia --project=. scripts/fast_ed.jl solve --config=PROFILE.toml --mu-index=N --sector-index=N [--force]
  julia --project=. scripts/fast_ed.jl collect --config=PROFILE.toml [--mu-index=N] [--allow-incomplete]
  julia --project=. scripts/fast_ed.jl compare --config=PROFILE.toml --mu-index=N

Sector index order: 1=(+,+), 2=(+,-), 3=(-,+), 4=(-,-).
""")
end

isempty(ARGS) && (usage(); exit(1))
command = lowercase(ARGS[1])
options = parse_options(ARGS[2:end])
haskey(options, "config") || throw(ArgumentError("Missing required option --config=..."))
config_path = isabspath(options["config"]) ? options["config"] :
    normpath(joinpath(PROJECT_ROOT, options["config"]))
spec = load_spec(config_path)

if command == "plan"
    println(plan(spec))
elseif command == "prepare"
    manifest = prepare_caches(spec; force=option_bool(options, "force"))
    println("cache ready: $(spec.cache_directory)")
    println("sectors: $(length(manifest["sectors"]))")
elseif command == "solve"
    path = solve_sector(
        spec,
        required_int(options, "mu-index"),
        required_int(options, "sector-index");
        force=option_bool(options, "force"),
    )
    println("sector result: $path")
elseif command == "collect"
    if haskey(options, "mu-index")
        result = collect_mu(spec, parse(Int, options["mu-index"]))
        println("score valid=$(result.score.valid) q=$(result.score.q) reason=$(result.score.reason)")
    else
        result = collect_all(spec; allow_incomplete=option_bool(options, "allow-incomplete"))
        println("collected $(length(result.rows)) / $(length(spec.mus)) mu points")
        isempty(result.failures) || foreach(message -> println(stderr, message), result.failures)
    end
elseif command == "compare"
    result = compare_direct(spec, required_int(options, "mu-index"))
    println(result)
    result.passed || error("cached and direct spectra did not pass comparison tolerances")
else
    usage()
    throw(ArgumentError("Unknown command '$command'"))
end
