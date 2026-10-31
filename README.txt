SUPPLEMENTARY MATERIAL
Interaction-topology-dependent energy redirection and projectile-origin
mass-transmission thresholds in hypervelocity impact of multilayer graphene

Package version: 1.0
Prepared: 2026-10-30, Asia/Shanghai
Archive layout: in/, potentials/, structures/, analysis/, results/, README.txt
All portable paths are relative to this directory after extraction.

1. DATA COVERAGE AND LIMITS

This package contains the inputs and 53 archived equilibrated configurations
for the 230 formal runs directly underlying the paper's main figures:
  P1 mechanism:       20 runs (four cells, n=5, 12 km/s)
  P4 sensitivity:    60 run entries (20 settings, n=3)
  P2-F parent:       80 runs
  P2-F round 1:      40 runs
  P2-F round 2:      30 runs

Terminal configurations, energy tracks, initial/final scalar records and
run logs are supplied for 200 runs (P1/P4 80; P2-F parent/round1 120).
For the 20 P1 trajectories, first and last dump frames were extracted
verbatim into first_last_frames.lammpstrj. They support the endpoint
undercoordination metric. Full trajectory time series are not bundled.

The original simulation outputs for all 30 P2-F round-2 runs were NOT found
in the supplied local archive. Their archived processed observations,
inputs, metadata, and ten equilibrated configurations ARE supplied.
results/missing_raw_outputs.csv identifies every affected run.
Thus the pooled threshold fit can be reproduced from all 150 processed
observations, but only 120 terminal fractions can currently be independently
re-counted from supplied original final.data files. This package does not
constitute a complete raw-trajectory archive for all 150 P2-F runs.

The ten round-2 base-preparation input files were also absent in the source
archive. Their code was reconstructed from the archived generate_campaign.py
function and exact run metadata, clearly labelled in archived_inputs and
structures/structure_manifest.csv. The equilibrated data files themselves
were copied unchanged from the archive. No simulation result was fabricated
or replaced by a new simulation.

The preserved source/p1_p4_results_all.json additionally contains 32 historical
screen entries. They do not enter the 20-run main factorial analysis, the
60-entry P4 analysis or the 150-run P2-F fit. Screen inputs/raw trajectories
are outside this compact main-figure package; these observations are not
silently pooled with the formal datasets.

2. SOFTWARE AND FORCE FIELDS

Archived run logs identify LAMMPS 29 Aug 2024 - Update 2, units metal,
atom_style full, boundary p p f. Target types are 1/2/3; projectile-origin
atoms are type 4. All supplied initial structures contain velocities.

Required LAMMPS styles: rebo, airebo (P4 variant), zbl,
kolmogorov/crespi/z, hybrid/overlay and table.
potentials/ contains CH.rebo, CH.airebo, CC.KC and all five switched LJ tables.
Original potential-file headers and the archived potential README are kept.
The switched LJ state is off below 1.90 A, on through 1.90-2.50 A,
and off by 10.20 A, as recorded in the archived specifications.

Numerical analysis scripts require Python >=3.10 and only its standard
library. Executed packaging verification used Python 3.12.14.
Optional figure recreation requires Matplotlib:
  python -m pip install -r analysis/requirements.txt
The original full-trajectory analysis scripts under analysis/original/
additionally use NumPy and SciPy. They retain their original working-tree
assumptions; use the portable top-level analysis scripts for this package.

3. RUN ORDER AND SEEDS

The archived protocol has no independent energy-minimization step.
in.relax is a dispatcher for geometry/velocity preparation, included by
in.equil; in.equil runs the archived NVT equilibration and zeros edge
velocities; in.impact reads an equilibrated structure, launches the projectile,
then performs NVE approach/impact and recovery. The archived baseline uses
300 K, 0.2 ps NVT equilibration, timestep 0.0001 ps and 5.0 ps recovery.
P4 overrides are recorded per run in in/run_specs.json.

Portable input files preserve numerical commands and alter only file paths.
The stage split occurs before fix eqt; in.equil automatically includes
in.relax so geometry and equilibration share a LAMMPS invocation.
Original combined inputs are retained under in/archived_inputs/.

For a single run using a supplied equilibrated structure, from this directory:
  python analysis/run_simulation.py --run-id mechanism_R0_L0_v12000_s51001_9e364d
This prints the command. Add --execute to actually run LAMMPS.
Use --lmp '"/path/to/lmp"' if lmp is not on PATH.

To regenerate a base and run an impact from that newly prepared base:
  python analysis/run_simulation.py --run-id mechanism_R0_L0_v12000_s51001_9e364d --stage equil --execute
  python analysis/run_simulation.py --run-id mechanism_R0_L0_v12000_s51001_9e364d --stage impact --new-base --execute

All new simulation outputs go to results/rerun/; archived observations remain
unchanged. Re-running round 2 creates new results and cannot recover the
missing historical outputs. MPI decomposition, software versions and platform
may affect trajectory-level numerical reproducibility.

Random target-velocity seeds:
  P1 mechanism: 51001-51005, matched across the four cells
  P4 sensitivity: 61001-61003
  P2-F parent: 71001-71010
  P2-F round 1: 91001-91010
  P2-F round 2: 101001-101010
The original base generator uses projectile seed = target seed + 100003.
Exact offsets, geometry, seed, timestep, force-field state and input mapping
are in in/run_specs.json and results/run_manifest.csv. There are 53 unique
base configurations across campaigns, shared where the protocol matches
initial states across R/L cells or velocities.

4. REPRODUCING THE PUBLISHED NUMBERS

  python analysis/verify_package.py
  python analysis/energy_loss.py
  python analysis/vt50_fit.py

New numerical output is written to results/recomputed/ and can be compared
with the supplied publication CSV files directly under results/.

energy_loss.py re-reads the original initial scalar and last track record for
all 80 P1/P4 runs and recomputes loss = initial projectile COM KE - final COM
KE. For the 20 P1 runs it also re-counts the endpoint undercoordination from
the two retained dump frames: same-layer neighbors within 1.85 A, excluding
edge atom IDs defined in the initial frame, final count minus baseline count.
The portable spatial-hash implementation uses the original non-periodic
Euclidean neighbor rule. Mean/error bars use sample SD (denominator n-1).

vt50_fit.py re-counts 120 available final configurations, pools these with
the 30 archived round-2 processed records, groups events at each speed, and
uses the unmodified logistic_fit function extracted from the original
round-2 analysis. The ordinary grouped-binomial model is:
  logit[p(v)] = alpha + beta * (v_mps - 21500)/1000
  V_T,50_mps = 21500 - 1000 * alpha/beta
Group proportions use Wilson 95% intervals. Threshold intervals are nominal
Wald 95% intervals from the fitted covariance matrix and delta method.
The point and interval calculations are conditional on the sequentially
chosen velocities, endpoint and model; no pairwise significance/equivalence
test is introduced.

Primary endpoint: fraction of original projectile atoms (type 4; n=471)
with terminal z <= -500 A after the preregistered approach plus 5.0 ps NVE
recovery. An event requires fraction >=0.50. The baseline target has 4116
atoms. This is an origin-resolved mass-transmission threshold under
fragmentation. It is not conventional ballistic V50, an intact-projectile
perforation criterion, or evidence for chemical bonds.
Projectile COM KE loss is not identical to target-absorbed energy.

5. FIGURE AND TABLE MAP

Fig.1 / Table 3 energy and structure columns:
  results/Fig1_individual_runs.csv          20 individual trajectories
  results/Fig1_summary.csv                   four means and sample SDs
  results/P1_factorial_contrasts.csv         descriptive energy contrasts
  analysis/fig1.py

Fig.2 / Table 2 / Table 3 threshold column:
  results/Fig2_pooled_150_runs.csv            150 event/fraction records
  results/Fig2_grouped_events.csv             15 groups and Wilson intervals
  results/Fig2_VT50_fit.csv                    four fits and Wald intervals
  results/Fig2_logistic_curves.csv             plotted fitted curves
  results/threshold_descriptive_comparisons.csv
  analysis/fig2.py

Fig.3 sensitivity:
  results/Fig3_individual_runs.csv             60 run entries
  results/Fig3_sensitivity_summary.csv         20 settings, mean and SD
  analysis/fig3.py
The seven axes are sheet size, edge width, timestep, recovery duration,
LJ epsilon, LJ sigma and reactive model. The 60 A reference uses three
sensitivity baseline trajectories (mean 4197.2036 eV), a different sample
from the five-trajectory P1 R1L1 mean (4188.6743 eV).

Other tables and checks:
  results/Table1_design.csv                   factor design
  results/Table3_combined.csv                  cross-endpoint summary
  results/P3_energy_drift.csv                  20 signed and absolute drifts
  results/P3_summary.csv                       maximum absolute drift
  results/P1_P4_raw_checks.csv                 scalar consistency checks
  results/DATA_DICTIONARY.txt                  fields, units and definitions

Optional plotting:
  python analysis/fig1.py
  python analysis/fig2.py
  python analysis/fig3.py
Recreated PNG/PDF files go to results/recreated_figures/. Original figure
PNG/PDFs, identical to the three manuscript embedded images, are kept in
results/figures_original/. Plot functions retain the original plotting logic
and have portable paths. Font availability can change figure appearance.

6. SOURCE INTEGRITY AND VERIFICATION SCOPE

results/source/ keeps the original unmodified JSON results, locked parent
snapshots, round-2 gate, input validation and figure-source manifests.
The historical round-2 gate says pass for threshold estimability; it is not
a new assertion that the currently missing original raw outputs are supplied.
results/source_file_manifest.csv records original relative archive paths,
source SHA-256, package SHA-256 and the preparation method. CHECKSUMS.sha256
covers every delivered file except itself. Source copies are not additional
independent simulation observations.

Numerical scripts were executed: all 80 P1/P4 energy results, all 20 P1
structural endpoints, 120 P2-F terminal fractions, 15 pooled event groups,
four threshold point estimates and their intervals match the archived
results. LAMMPS zero-step compatibility checks were executed with LAMMPS
4 Feb 2025 for separated REBO+LJ, shared REBO+LJ and the AIREBO model variant.
The unchanged switched LJ table produced a warning at 2 of 20000 table
entries about force/energy consistency; no potential file was altered.
These checks do not constitute new full-duration MD simulations or universal
force-field validation. Plot scripts were syntax-checked; Matplotlib rendering
was not executed on the packaging host. The exact original figures are supplied.
See results/package_validation.json for the machine-readable verification scope.
