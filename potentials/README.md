# 势函数文件溯源

本目录只保存新 P1–P4 管线使用的版本化副本。生成脚本在 `generated/run_manifest.json` 中记录每个文件的 SHA-256。

- `CH.rebo` 与 `CC.KC`：从本项目既有、已验证的 `论文材料/实验1/case_A_v4000/` 复制；原始文件不作修改。
- `CH.airebo`：来自本机 LAMMPS 安装自带的标准势函数目录，仅用于 P4 的 `airebo 3.0 0 0`（REBO-only）敏感性分析。
- `switched_cc_lj_*.table`：由 `scripts/generate_campaign.py` 按各项参数生成；不是外部势函数文件。

任何替换势函数都必须另存于本目录、在配置中明确命名，并重新生成 manifest；不得静默覆盖文件。
