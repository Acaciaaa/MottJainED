# N=7 有界自动分段流程

`Uf0=2.00, Vf0=0.55, V0=0.34` pilot和`Vf0=0.45`生产点均已完成并通过本地raw复核；
当前生产点是`Uf0=1.834, Vf0=0.65, V0=0.34`。流程复用已经验证过的
“共享矩阵 cache + 独立 `(mu,sector)` array”算法，只自动执行每批结果之后的小文件审计、
下一批局部网格选择和 Slurm 提交。它不会调用旧 `FastMuSearch`、不会做全区间优化，
也不会在等待时占着 CPU。

## 决策顺序

1. 用每个Hamiltonian自己的 N5/N6 阻尼延拓 `mu6 + 0.5*(mu6-mu5)` 作中心，实际计算中心左右
   `[-0.01,-0.005,0,0.005,0.01]` 五个 scout 点。若N5/N6已经确认另有竞争谷，可在首次
   scout中加入最多两个明确的`guard_mus`，这些点与主网格一起实际求解并参与全局比较。
2. 若 scout 最低点在边界，只向下降方向增加至多两个 `0.005` 间隔的点；最低点被包住后，
   用相邻三点的 `q^2` 二次拟合选择网格中心。拟合值只选网格，不会写成实测 `muc`。
3. 以 `0.00025` 间隔实际计算七个 refine 点，中心吸附到 `0.000125` 网格。
4. 只有实测最低点左右都有距离不超过 `0.000251` 的有效点、两侧 q 均回升、二次拟合
   对 q 的改进不超过 `1e-5`，且 factor/五条 raw gap 局部连续时才接受。否则只补一至
   两个点。
5. 最多 16 个不同 mu、最多三轮自适应提交，且 mu 必须留在 `(0.12,0.18)`。
   达到边界/预算、缺谱、身份冲突或连续性失败时立即进入 `review`，不再提交 ED。

接受的是实际算过的最低点。最终审计明确保留“有限网格结果，不证明全局最低点或热力学
临界点”的限定。

## 每轮自动核对

控制任务重新读取原始 CSV/TOML，并检查四个 sector 的 cache/settings/job identity、
Hamiltonian、cold-start k=20 solver、线程数、sector 维数和完整 rank。评分固定为
`ds_s,j,curlj,dj_rank1,t_rank1` 五项 q；S 使用 `(0,0)` raw rank 2，curlJ 使用
`(2,3)` raw rank 3，其余约定保持不变。它还检查基态 singlet、量子数取整误差、最低
对称副本劈裂、factor 和五条 gap 的局部变化。无效分数不会被替换成“大 q”。

这里J和curlJ属于同一个`(L2,C2)=(2,3)`量子数通道：raw rank 1、2是J这一条物理能级在
两个离散对称sector中的等能副本；raw rank 3、4是下一条不同物理能级curlJ的两个副本。
因此curlJ写raw rank 3，去除副本后就是该量子数通道的第二条物理能级。

## 并发和计费

- prepare：一个 8 CPU/8 thread 任务，完成后释放；实测峰值约 23 GiB。
- solve：每个 array element 只算一个 `(mu,sector)`，8 CPU/8 thread；最多四个同时运行，
  每个完成后独立释放。实测峰值约 15.3 GiB。
- controller：1 CPU，仅读取小文件、写下一轮 profile、提交作业并退出。
- 等待中的 Slurm dependency 没有获得节点，不计为驻留的 32 CPU。

当前 8 CPU 已给约 62 GiB 内存，未发现需要扩大申请的证据。若将来真实 OOM，可以在启动
时设置 `PREPARE_CPUS_OVERRIDE` 或 `SOLVE_CPUS_OVERRIDE` 只增加内存额度；计算线程仍按
profile 固定为 8。`MAX_CONCURRENT_OVERRIDE` 可以降低同一时刻的独立任务数。

## 结果和停止点

通过验收后，controller 自动写最终 summary、relations、固定 raw-rank 选态、最佳完整谱、
源文件 SHA-256、含可信 N5/N6 与新 N7 行的 `fss_n567.csv` 和 audit，并生成：

```text
output/fast_ed/archives/<pipeline_name>_final.tar.gz
output/fast_ed/archives/<pipeline_name>_final.tar.gz.sha256
```

归档不再收入打包前的live状态文件。它收入`final/pipeline_state.toml`完成态快照；服务器上的
live状态只在归档SHA验证后才原子更新为`complete`并记录实际SHA。这样归档内部不会再出现
`action=bundle, complete=false`的误导状态，同时也不会在tar失败前把服务器live状态提前标成
完成。正常完成只需下载上面两个文件；日志和sacct仅在失败或资源诊断时需要。

## 服务器入口

Vf0=0.45最终包已在本地逐文件核验，服务器cache ID为`342e15745d81a91a`。Vf0=0.65配置
列出已经核验过的旧cache retirement profiles。用户入口先用纯shell检查队列并只提交一个
1 CPU bootstrap；bootstrap在计算节点完成Julia预编译、旧cache的ID/完整manifest核验与释放，
然后启动prepare、scout和后续controller。登录节点不再启动Julia，也不需要手工执行cache检查。
retained中心cache仍由底层命令独立保护。

该点N5/N6为`mu=0.13608278512369/0.13807344943979`，阻尼中心`0.13906878159784`，
五点主scout范围为`0.12906878159784–0.14906878159784`。N6另有一个已确认但q较高的局部谷
`mu=0.12771774983325`，所以首次额外实算这一个guard点，共6个mu、24个sector task：

```bash
bash scripts/submit_fast_ed_pipeline.sh config/fast_ed/n7_vf0_065_auto.toml
```

状态文件位于
`output/fast_ed/pipelines/n7_vf0_065_auto/pipeline_state.toml`。若状态为 `review`，说明流程
已停止且 cache 保留；不要直接重启或删除文件，应根据 `review_reason` 做一次针对性处理。
