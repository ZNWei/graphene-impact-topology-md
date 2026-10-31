# P1–P4 数据字典

## `track.dat`

| 列 | 单位 | 定义 |
|---|---:|---|
| `step`, `time_ps` | step, ps | LAMMPS 时间步和模拟时间 |
| `stage` | — | 0 为连续 NVE 记录；冲击/恢复边界由 input metadata 的 `impact_steps` 重建 |
| `ke_all_eV`, `pe_all_eV`, `etotal_eV` | eV | 全系统能量及闭合漂移 |
| `ke_target_eV`, `ke_projectile_eV` | eV | 靶材与弹丸总动能 |
| `ke_projectile_cm_eV`, `ke_projectile_nontrans_eV` | eV | 弹丸质心平动和非平动动能；后者不是完整 internal energy，因为不含势能 |
| `pe_target_assigned_eV`, `pe_projectile_assigned_eV` | eV | 由 `pe/atom` 求和的每原子势能归属；跨界面相互作用按 LAMMPS 原子能量分配，不能与精确 pair partition 混同 |
| `e_reactive_eV`, `e_zbl_eV`, `e_kc_eV`, `e_lj_tail_eV` | eV | 对应 pair style 的全系统势能项 |
| `target_zmin_A`, `target_zmax_A`, `projectile_zmin_A`, `projectile_zmax_A` | Å | 终态全体包络，用于完全分离判据 |
| `projectile_zcm_A`, `projectile_vz_mps` | Å, m/s | 弹丸质心位置和轴向速度 |

## `trajectory.lammpstrj`

`type=1,2,3` 是原始 layer 标签，`type=4` 是弹丸。分析器不得通过当前 z 坐标重标 layer。逐原子 `c_pe_atom` 用于势能归属；它不直接给出 REBO bond order。

## `run_results_all.json`

- `outcome` 仅取 `perforated`、`rebounded`、`unresolved`，由全弹丸/全靶包络、相互作用 cutoff 和速度平台联合决定。
- `projectile_*_loss` 是弹丸动能损失，不是靶材吸收能。
- `structure.projectile_target_*` 与 `structure.layer_layer_*` 是**距离阈值短程接触**和持续时间，不是化学键。
- `structure.target_damage_edge_excluded_net` 是仅靶材同层配位、扣除初始边缘基线后的代理量。
- `local_interlayer_gap_time_series` 是撞击中心局部相邻层的近邻投影间距分布，非整层 COM 距离。
