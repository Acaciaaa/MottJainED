#!/usr/bin/env julia

import Pkg
const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT; io=devnull)

using TOML
using Dates
using MottJainED
include(joinpath(PROJECT_ROOT, "experimental", "FastED.jl"))
include(joinpath(PROJECT_ROOT, "experimental", "FastEDPipeline.jl"))
using .FastED
using .FastEDPipeline

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

function required(options, key)
    haskey(options, key) || throw(ArgumentError("Missing required option --$key=..."))
    return options[key]
end

project_path(path) = isabspath(path) ? normpath(path) : normpath(joinpath(PROJECT_ROOT, path))
config_path(options) = project_path(required(options, "config"))
print_fields(fields) = println(join(string.(fields), '\t'))

function directory_bytes(path)
    isdir(path) || return 0
    total = 0
    for (directory, _, names) in walkdir(path)
        for name in names
            file = joinpath(directory, name)
            isfile(file) && (total += filesize(file))
        end
    end
    return total
end

function inspect_cache(path)
    spec = FastED.load_spec(path)
    present = isdir(spec.cache_directory)
    complete = false
    if present && isfile(FastED.cache_manifest_path(spec))
        manifest = TOML.parsefile(FastED.cache_manifest_path(spec))
        complete = Bool(get(manifest, "complete", false)) &&
                   String(get(manifest, "cache_id", "")) == spec.cache_id
    end
    print_fields((spec.cache_id, present, complete, directory_bytes(spec.cache_directory),
                  spec.cache_directory))
end

function retire_configured_caches(path)
    config = TOML.parsefile(path)
    pipeline = get(config, "pipeline", nothing)
    pipeline isa AbstractDict || throw(ArgumentError("Profile must define [pipeline]"))
    configured = get(pipeline, "retire_before_start", Any[])
    configured isa AbstractVector || throw(ArgumentError(
        "pipeline.retire_before_start must be an array of retirement profiles",
    ))

    for entry in configured
        retirement_path = project_path(String(entry))
        spec = FastED.load_spec(retirement_path)
        retirement = get(spec.config, "cache_retirement", nothing)
        retirement isa AbstractDict || throw(ArgumentError(
            "Configured cleanup is not an approved cache-retirement profile: $retirement_path",
        ))
        expected = String(retirement["expected_cache_id"])
        expected == spec.cache_id || throw(ArgumentError(
            "Configured cleanup identity does not match its retirement profile: $retirement_path",
        ))

        if !isdir(spec.cache_directory)
            print_fields(("cache_absent", expected, spec.cache_directory))
            continue
        end
        manifest_path = FastED.cache_manifest_path(spec)
        isfile(manifest_path) || throw(ArgumentError(
            "Configured cleanup cache has no manifest: $(spec.cache_directory)",
        ))
        manifest = TOML.parsefile(manifest_path)
        Bool(get(manifest, "complete", false)) || throw(ArgumentError(
            "Configured cleanup cache is incomplete: $(spec.cache_directory)",
        ))
        String(get(manifest, "cache_id", "")) == expected || throw(ArgumentError(
            "Configured cleanup cache manifest identity mismatch: $(spec.cache_directory)",
        ))
        retire_cache(retirement_path, Dict("confirm-cache-id" => expected))
    end
end

function retire_cache(path, options)
    spec = FastED.load_spec(path)
    retirement = get(spec.config, "cache_retirement", nothing)
    retirement isa AbstractDict || throw(ArgumentError(
        "The profile is not an approved cache-retirement profile",
    ))
    expected = required(options, "confirm-cache-id")
    configured = String(retirement["expected_cache_id"])
    expected == configured == spec.cache_id || throw(ArgumentError(
        "The typed cache ID, profile cache ID, and retirement ID must match",
    ))
    protected_profile = joinpath(PROJECT_ROOT, "config", "fast_ed", "n7_retained_k20.toml")
    protected_id = FastED.load_spec(protected_profile).cache_id
    spec.cache_id != protected_id || throw(ArgumentError("The retained-center cache is protected"))
    isdir(spec.cache_directory) || throw(ArgumentError("Cache is already absent: $(spec.cache_directory)"))
    islink(spec.cache_directory) && throw(ArgumentError("Refusing to remove a symlinked cache"))
    manifest = FastED.load_cache_manifest(spec)
    String(manifest["cache_id"]) == expected || throw(ArgumentError("Manifest cache ID mismatch"))
    temporary = String[]
    for (directory, _, names) in walkdir(spec.cache_directory)
        append!(temporary, joinpath.(Ref(directory), filter(name -> occursin(".tmp-", name), names)))
    end
    isempty(temporary) || throw(ArgumentError("Temporary cache files exist; inspect before retirement"))
    bytes = directory_bytes(spec.cache_directory)
    receipt = Dict{String,Any}(
        "cache_id" => spec.cache_id,
        "cache_directory" => spec.cache_directory,
        "removed_bytes" => bytes,
        "verified_local_archive_sha256" => String(retirement["verified_local_archive_sha256"]),
        "hamiltonian" => FastED.coupling_dict(spec.couplings; include_mu=false),
        "status" => "approved",
        "approved_at" => string(now()),
    )
    receipt_path = joinpath(PROJECT_ROOT, "output", "fast_ed", "cache_retirements",
                            "$(spec.cache_id).toml")
    MottJainED.atomic_toml(receipt_path, receipt)
    rm(spec.cache_directory; recursive=true)
    receipt["status"] = "complete"
    receipt["retired_at"] = string(now())
    MottJainED.atomic_toml(receipt_path, receipt)
    print_fields((spec.cache_id, bytes, receipt_path))
end

function usage()
    println("""
Usage:
  julia --project=. scripts/fast_ed_pipeline.jl init --config=PROFILE
  julia --project=. scripts/fast_ed_pipeline.jl claim-launch --config=PROFILE
  julia --project=. scripts/fast_ed_pipeline.jl reset-launch --config=PROFILE
  julia --project=. scripts/fast_ed_pipeline.jl record-submission --config=PROFILE --stage=NAME --solve-job=ID --controller-job=ID [--prepare-job=ID]
  julia --project=. scripts/fast_ed_pipeline.jl advance --config=PROFILE
  julia --project=. scripts/fast_ed_pipeline.jl action --config=PROFILE
  julia --project=. scripts/fast_ed_pipeline.jl resources --config=PROFILE
  julia --project=. scripts/fast_ed_pipeline.jl review --config=PROFILE --reason=TEXT
  julia --project=. scripts/fast_ed_pipeline.jl resume-stalled-fit --config=PROFILE
  julia --project=. scripts/fast_ed_pipeline.jl mark-bundled --config=PROFILE --archive=PATH --sha256=HEX
  julia --project=. scripts/fast_ed_pipeline.jl inspect-cache --config=PROFILE
  julia --project=. scripts/fast_ed_pipeline.jl retire-cache --config=PROFILE --confirm-cache-id=ID
  julia --project=. scripts/fast_ed_pipeline.jl retire-configured-caches --config=PIPELINE_PROFILE

Machine-readable action and resource commands print tab-separated fields.
""")
end

isempty(ARGS) && (usage(); exit(1))
command = lowercase(ARGS[1])
options = parse_options(ARGS[2:end])
path = config_path(options)

if command == "init"
    state = FastEDPipeline.initialize(path)
    println(FastEDPipeline.state_path(FastEDPipeline.pipeline_options(FastED.load_spec(path))))
elseif command == "claim-launch"
    FastEDPipeline.claim_launch(path)
    print_fields(FastEDPipeline.action_fields(path))
elseif command == "reset-launch"
    FastEDPipeline.reset_launch(path)
    print_fields(FastEDPipeline.action_fields(path))
elseif command == "record-submission"
    FastEDPipeline.record_submission(
        path; stage=required(options, "stage"), solve_job=required(options, "solve-job"),
        controller_job=required(options, "controller-job"),
        prepare_job=get(options, "prepare-job", ""),
    )
    print_fields(FastEDPipeline.action_fields(path))
elseif command == "advance"
    try
        FastEDPipeline.advance(path)
    catch err
        FastEDPipeline.mark_review(path, "controller_analysis_failed")
        rethrow(err)
    end
    print_fields(FastEDPipeline.action_fields(path))
elseif command == "action"
    print_fields(FastEDPipeline.action_fields(path))
elseif command == "resources"
    print_fields(FastEDPipeline.resource_fields(path))
elseif command == "review"
    FastEDPipeline.mark_review(path, required(options, "reason"))
    print_fields(FastEDPipeline.action_fields(path))
elseif command == "resume-stalled-fit"
    FastEDPipeline.resume_stalled_fit(path)
    print_fields(FastEDPipeline.action_fields(path))
elseif command == "mark-bundled"
    archive = project_path(required(options, "archive"))
    FastEDPipeline.mark_bundled(path, archive, required(options, "sha256"))
    print_fields(FastEDPipeline.action_fields(path))
elseif command == "inspect-cache"
    inspect_cache(path)
elseif command == "retire-cache"
    retire_cache(path, options)
elseif command == "retire-configured-caches"
    retire_configured_caches(path)
else
    usage()
    throw(ArgumentError("Unknown command '$command'"))
end
