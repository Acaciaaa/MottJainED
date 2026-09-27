# Conformal v3：服务器更新与预检

只更新 `MottJainED-so3lver` 的 `conformal-algebra-projected` 分支。
**不要取消 `so3-n6-proj-par`；本次不修改它的脚本、配置或结果目录。**
已经完成的旧 conformal job 也不需要再取消。

## 1. 后台拉取，断网自动重试

```bash
cd /public/home/ruiqixu/MottJainED/MottJainED-so3lver
mkdir -p slurm-logs
nohup bash -c '
until git -c http.version=HTTP/1.1 fetch origin \
  refs/heads/conformal-algebra-projected:refs/remotes/origin/conformal-algebra-projected
do
  date
  echo "fetch failed; retrying in 30 seconds"
  sleep 30
done
echo "FETCH_COMPLETE"
' > slurm-logs/fetch-conformal-v3.log 2>&1 &
echo $! > slurm-logs/fetch-conformal-v3.pid
```

用下面的命令查看；**看到 `FETCH_COMPLETE` 后**再进行下一步。

```bash
tail -n 20 slurm-logs/fetch-conformal-v3.log
```

## 2. 切换并快进更新

单独的 bash 子进程会在检查失败时停止，不改变交互式 shell 设置。
有未提交的 tracked 改动时停下检查，不要 reset 或强制覆盖。

```bash
cd /public/home/ruiqixu/MottJainED/MottJainED-so3lver
bash <<'BASH'
set -euo pipefail
test -z "$(git status --porcelain --untracked-files=no)"
git cat-file -e origin/conformal-algebra-projected:experimental/SO3ConformalFit.jl
if git show-ref --verify --quiet refs/heads/conformal-algebra-projected; then
  git switch conformal-algebra-projected
else
  git switch --track -c conformal-algebra-projected origin/conformal-algebra-projected
fi
git merge --ff-only origin/conformal-algebra-projected
test "$(git rev-parse HEAD)" = "$(git rev-parse origin/conformal-algebra-projected)"
git status --short --branch
git log -2 --oneline
BASH
```

核对这里的最新提交与本次回复给出的提交一致。`?? slurm-logs/` 不影响更新。

## 3. 先提交一个 CPU 的固定点预检

```bash
cd /public/home/ruiqixu/MottJainED/MottJainED-so3lver
mkdir -p slurm-logs
sbatch --parsable slurm/so3lver_n6_conformal_pilot.sbatch \
  | tee slurm-logs/latest-so3-n6-conf-check.jobid
```

它只算 anchor、旧 corrected 搜索的最佳点，以及后者的冷启动复算。
**不进行参数搜索，不会自动启动六 CPU 的大任务。** 不设置 `--mem`/`--time`。
预检用于看新内层拟合是否收敛、复算是否一致、每点耗时和内存。
`pilot_passed=true` 不是找到了 CFT，也不是旧最佳点一定通过新的态连续性检查。

完成后从下面的目录下载两个文件（将 `JOBID` 替换为本次 job id）：

```text
/public/home/ruiqixu/MottJainED/MottJainED-so3lver/output/so3lver/conformal_optimization/
  n6_full_algebra_pilot_job-JOBID.tar.gz
  n6_full_algebra_pilot_job-JOBID.tar.gz.sha256
```

放到本地 `parameter_search` 即可。普通失败也会打包，请一并下载日志归档。

## 4. 预检分析后再启动完整搜索

```bash
sbatch --parsable slurm/so3lver_n6_conformal_optimize.sbatch \
  | tee slurm-logs/latest-so3-n6-conf-v3.jobid
```

仍是 **6 个独立 worker 同时搜索，每个 1 CPU**，不是一个进程串行跑 6 个起点。
新目录为 `n6_multistart_03_full_algebra`，不混用旧版本 checkpoint。
完成后下载同一父目录下的：

```text
n6_multistart_03_full_algebra_job-JOBID.tar.gz
n6_multistart_03_full_algebra_job-JOBID.tar.gz.sha256
```

若作业被强杀而没有归档，先保留整个新结果目录与 Slurm 日志，不要删除 checkpoint。
