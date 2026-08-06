# 本文件集中处理输出目录、CSV 断点续算、可复现元数据和稳定任务编号。
# 物理工作流不应在各处重复手写 CSV/JLD2 读写逻辑。

# 项目根目录（src 的上一级），不依赖用户从哪个当前目录启动 Julia。
const PACKAGE_ROOT = normpath(joinpath(@__DIR__, ".."))

"""建立输出目录并返回绝对路径。"""
function ensure_output(path::AbstractString)
    mkpath(path)
    return abspath(path)
end

"""先写临时文件再原子替换目标 CSV，降低中断时留下半个文件的风险。"""
function atomic_csv(path::AbstractString, table)
    mkpath(dirname(path))
    temporary = path * ".tmp-$(getpid())"
    CSV.write(temporary, table)
    mv(temporary, path; force=true)
    return path
end

"""向结果 CSV 追加一行；文件首次出现时自动写表头。"""
function append_csv(path::AbstractString, row::NamedTuple)
    mkpath(dirname(path))
    CSV.write(path, DataFrame([row]); append=isfile(path), writeheader=!isfile(path))
    return path
end

"""读取结果中 `status=ok` 的 job_id，用于参数扫描断点续跑。"""
function completed_job_ids(path::AbstractString)
    isfile(path) || return Set{String}()
    data = CSV.read(path, DataFrame)
    "job_id" in names(data) || return Set{String}()
    if "status" in names(data)
        return Set(String(data.job_id[i]) for i in 1:nrow(data) if String(data.status[i]) == "ok")
    end
    return Set(String.(data.job_id))
end

"""由任务的物理参数生成可复现的 16 位编号，与运行先后顺序无关。"""
function stable_id(parts...)
    normalized = join((repr(part) for part in parts), "|")
    return bytes2hex(sha1(normalized))[1:16]
end

"""把 `Couplings` 展开成可直接写入 DataFrame 的一行字段。"""
function coupling_namedtuple(c::Couplings)
    return (
        Uf=c.Uf, Uf0=c.Uf0, U0=c.U0, Vf=c.Vf, Vf0=c.Vf0,
        V0=c.V0, t=c.t, mu_initial=c.mu,
    )
end

"""读取一个目录当前的 Git commit；不是 Git repo 时返回 `unknown`。"""
function git_revision(path::AbstractString)
    try
        return readchomp(pipeline(`git -C $path rev-parse --short HEAD`; stderr=devnull))
    catch
        return "unknown"
    end
end

"""
记录命令、配置路径、Julia/线程版本及两个项目的 Git revision。

每个输出目录写 `run_metadata.toml`，用于以后回答“这批数据是怎么跑出来的”。
"""
function write_run_metadata(output::AbstractString; command::AbstractString, config_path=nothing)
    metadata = Dict{String,Any}(
        "command" => String(command),
        "created_at" => string(now()),
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "package_version" => "0.1.0",
        "package_git_revision" => git_revision(PACKAGE_ROOT),
        "fuzzified_path" => normpath(joinpath(PACKAGE_ROOT, "..", "FuzzifiED.jl")),
        "fuzzified_git_revision" => git_revision(normpath(joinpath(PACKAGE_ROOT, "..", "FuzzifiED.jl"))),
    )
    config_path === nothing || (metadata["config_path"] = abspath(config_path))
    path = joinpath(output, "run_metadata.toml")
    open(path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    return path
end

"""同一 job_id 被重跑多次时只保留最后一行，便于后续拟合。"""
function latest_rows(data::DataFrame; key::Symbol=:job_id)
    key_string = String(key)
    key_string in names(data) || return data
    indices = Dict{String,Int}()
    for i in 1:nrow(data)
        indices[String(data[i, key])] = i
    end
    return data[sort!(collect(values(indices))), :]
end

"""把异常压成单行文字，保证能安全写入 CSV。"""
sanitize_error(err) = replace(sprint(showerror, err), '\n' => ' ', '\r' => ' ')
