#!/usr/bin/env julia

import Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT; io=devnull)

using CSV
using DataFrames
using TOML

function parse_options(args)
    options = Dict{String,String}()
    for argument in args
        startswith(argument, "--") || throw(ArgumentError(
            "Unknown positional argument '$argument'; use --key=value options",
        ))
        parts = split(argument[3:end], "="; limit=2)
        length(parts) == 2 || throw(ArgumentError("Option '$argument' needs a value"))
        options[parts[1]] = parts[2]
    end
    return options
end

function project_path(path)
    return isabspath(path) ? normpath(path) : normpath(joinpath(PROJECT_ROOT, path))
end

options = parse_options(ARGS)
haskey(options, "config") || throw(ArgumentError("Missing --config=PROFILE.toml"))
haskey(options, "source") || throw(ArgumentError("Missing --source=RESULTS.csv"))
config = TOML.parsefile(project_path(options["config"]))
source = project_path(options["source"])
isfile(source) || throw(ArgumentError("FSS results not found: $source"))

expected_sizes = Int.(config["model"]["nm_values"])
fss = config["fss"]
expected_values = Float64.(fss["scan_values"])
parameter = String(fss["scan_parameter"])
data = CSV.read(source, DataFrame)

required_columns = (
    "job_id", "status", "nm1", "scan_parameter", "scan_value", "muc", "q",
    "factor", "delta_s", "delta_o", "score_valid", "completed", "at_boundary",
    "wide_brent_ran", "wide_disagreement",
)
for name in required_columns
    name in names(data) || error("Missing required FSS column '$name' in $source")
end

# A restarted job can append a newer row with the same identity. Audit only the last row.
latest = combine(groupby(data, :job_id), group -> last(group))
issues = String[]
selected_rows = DataFrame[]
for nm1 in expected_sizes, value in expected_values
    matches = filter(row ->
        Int(row.nm1) == nm1 &&
        String(row.scan_parameter) == parameter &&
        Float64(row.scan_value) == value,
        latest,
    )
    if nrow(matches) != 1
        push!(issues, "N=$nm1 $parameter=$value has $(nrow(matches)) latest rows")
        continue
    end
    row = matches[1, :]
    push!(selected_rows, matches)
    String(row.status) == "ok" || push!(issues, "N=$nm1 $parameter=$value status=$(row.status)")
    Bool(row.score_valid) || push!(issues, "N=$nm1 $parameter=$value has an invalid score")
    Bool(row.completed) || push!(issues, "N=$nm1 $parameter=$value search did not converge")
    Bool(row.at_boundary) && push!(issues, "N=$nm1 $parameter=$value best μ is on a boundary")
    all(isfinite, Float64[row.muc, row.q, row.factor, row.delta_s, row.delta_o]) ||
        push!(issues, "N=$nm1 $parameter=$value has a non-finite result")
    if nm1 == maximum(expected_sizes)
        Bool(row.wide_brent_ran) || push!(issues, "N=$nm1 $parameter=$value wide audit did not run")
        Bool(row.wide_disagreement) && push!(
            issues, "N=$nm1 $parameter=$value local/wide μ searches disagree",
        )
    end
end

isempty(selected_rows) || println(vcat(selected_rows...)[
    :, [:nm1, :scan_parameter, :scan_value, :muc, :q, :factor, :delta_s, :delta_o,
        :best_source, :wide_disagreement],
])
isempty(issues) || error("FSS audit failed:\n" * join(issues, "\n"))
println("FSS audit passed: $source")
