# Changes to the host model and coupling since `tend-nudging2`

Branch `fix-coupling-residual`, five commits on top of `tend-nudging2` (a760842).
Written for Kairui, 2026-09-14. Everything stays on this branch; nothing is merged to
`main` (`origin/main` is at 2025-09-06 and predates most of the host-model work).

Repo: <https://github.com/Qiyu-Song/SAM-MMF-Walker>

---

## READ THIS FIRST — one change breaks existing `prm` files

[`a85416b`](https://github.com/Qiyu-Song/SAM-MMF-Walker/commit/a85416b7075bed16d5cb29a9827e0a2c1e37119f)
**removed** two namelist variables and **added** one:

| removed from `KUANG_PARAMS` | added |
|---|---|
| `inverse_prefilter_k1_fraction` | `suppress_k_start` |
| `inverse_prefilter_k2_fraction` | |

Fortran's namelist read **fails on an unrecognised name**. If you build this branch and
run an old `prm`, SAM stops during `setparm` with exit status 19 and no useful message.
The reverse also fails: a new `prm` against an old executable. This cost us two runs, once
in each direction.

**How to convert.** `suppress_k_start` is the lowest wavenumber the coupling filter
suppresses, in host-column index space (1 … nsx/2):

    suppress_k_start = nint( inverse_prefilter_k1_fraction * nsx/2 )

`-1` means "resolve to nsx/2 in `setparm`", i.e. suppress only the Nyquist wavenumber.
At nsx = 160: `k1 = 0.95` → `76`; `k1 = 1.00` → `80`, or equivalently `-1`.
`setparm` validates the range and aborts with a clear message if it is outside 1 … nsx/2.

`a85416b` also moved `dx_hm` out of the namelist into `domain.f90` as `dx_hm_km`, next to
`nx_gl` and `nsubdomains_x`, so the three numbers that have to be consistent now sit
together. Remove `dx_hm` from your `prm` and set `dx_hm_km` in `domain.f90` instead.

**Worth adopting**: before submitting, check every key in a `prm` against the executable:

```bash
strings <exe> | tr 'A-Z' 'a-z' > /tmp/exe.txt
for k in $(sed -nE 's/^[[:space:]]*([A-Za-z_][A-Za-z0-9_]*)[[:space:]]*=.*/\1/p' prm.X \
           | tr 'A-Z' 'a-z' | sort -u); do
  grep -q "$k" /tmp/exe.txt || echo "UNKNOWN KEY: $k"
done
```

---

## `domain.f90` now carries `dx_hm_km` — set it, and check it every time

[`a85416b`](https://github.com/Qiyu-Song/SAM-MMF-Walker/commit/a85416b7075bed16d5cb29a9827e0a2c1e37119f)
moved the host grid spacing out of the namelist into `domain.f90`, so **three numbers
that must be consistent with each other now live together**:

```fortran
integer, parameter :: nx_gl          = 2560   ! CRM points in x, whole domain
integer, parameter :: nsubdomains_x  = 160    ! = nsx, the number of host columns
real,    parameter :: dx_hm_km       = 64.    ! host grid spacing [km]   <-- NEW
```

with, at `dx = 4 km`,

    subdomain width = (nx_gl / nsubdomains_x) * dx
    host domain     = nsubdomains_x * dx_hm_km

`setparm` prints all of these at startup. The configurations we have used:

| case | `nx_gl` | `nsubdomains_x` | `dx_hm_km` | subdomain | host domain |
|---|---|---|---|---|---|
| Walker | 2560 | 160 | 64. | 64 km | 10240 km |
| shear, dx_hm = 128 | 1024 | 32 | 128. | 128 km | 4096 km |
| shear, dx_hm = 64 | 2048 | 64 | 64. | 128 km | 4096 km |
| shear, dx_hm = 32 | 4096 | 128 | 32. | 128 km | 4096 km |
| shear, L_sub = 256 | 2048 | 32 | 128. | 256 km | 4096 km |

**Why this needs watching.** It replaced `dx_hm = dx * nx / 4.0` in `setparm`, which
silently tied the host spacing to a quarter of the subdomain width. Stating it explicitly
is the point of the change — but it also means **`dx_hm_km` will not follow when you edit
`nx_gl` or `nsubdomains_x`.** Change one, check all three, and read the geometry block
`setparm` prints.

**And the merge hazard.** `domain.f90` is per-experiment configuration, not code, so git
line-merges it *without raising a conflict* and produces a hybrid of two configurations
that compiles and runs. Merging `tend-nudging2` into this branch gives
`nx_gl = 4096` and `nsubdomains_x = 128` from one side with `dx_hm_km = 64.` from the
other — an 8192 km host domain that is nobody's intended setup, reported with no warning
beyond the startup print.

**This branch keeps the Walker configuration** (2560 / 160 / 64.). That discards nothing:
the only `domain.f90` change on `tend-nudging2` is `nx_gl` 2560 -> 4096 and
`nsubdomains_x` 160 -> 128, i.e. the dx_hm = 32 km shear setup, which you set per
experiment anyway. **Set the three numbers deliberately after any pull or merge; never
accept git's version of this file.**

## Merging into `tend-nudging2`

Already done from our side: `15b4990` is merged into `fix-coupling-residual`, with
`setparm.f90` resolved to `dx_hm = dx_hm_km * 1000.`, the `public ::` list in
`module_hostmodel.f90` keeping every name from both sides, and `domain.f90` kept as
above. So you can take this branch without untangling anything — just set `domain.f90`
for whatever you are running.


---

## The five commits

Four of the five are drop-in: at their default settings the model reproduces
`tend-nudging2` exactly. Only `a85416b` requires action, and only `fef8700` changes any
number (by 0.2%, in shear runs only).

### 1. [`79450ac`](https://github.com/Qiyu-Song/SAM-MMF-Walker/commit/79450ac90f52a3b0632778626efa15c9482e8252) — remove the coupling residual orphaned in subdomain U

The U round trip (host faces → subdomain centres → back) has spectral gain
`|cos(pi*m/nsx)|`, which is **exactly zero at m = nsx/2**. So part of the CRM-side U
adjustment cannot be handed to the host at all and stays stranded in the subdomain — this
is what produced the persistent stripes in the subdomain-mean U.

The fix forms, at every coupling step,

    u_resid = u_sub - return( send( u_sub ) )

using the **same** send/return operators the coupling already uses, and subtracts it from
the CRM-mean U. Its spectral weight is exactly `1 - h^2` where `h` is the coupling's own
transfer: identically zero wherever the host can respond, one where it is blind. No tuned
parameter and no new timescale.

The reason this is better than "remove the 2*dx mode" is that the blind set is not only
Nyquist. With the inherited shoulder at `k1 = 0.95` the partly-invisible band is
k = 77–80, and the weight follows the shoulder automatically
(`1 - h^2` = 0.27 / 0.75 / 0.98 / 1.00 at k = 77–80).

| namelist | default | effect at default |
|---|---|---|
| `do_remove_coupling_residual` | `.false.` | identical to `tend-nudging2` |
| `do_remove_nyquist_u` | `.false.` | identical |

**Recommended**: `do_remove_coupling_residual = .true.` with `suppress_k_start = -1`
(Nyquist only). The projection makes the wider shoulder unnecessary, and the shoulder
costs real signal across its width.

### 2. [`a85416b`](https://github.com/Qiyu-Song/SAM-MMF-Walker/commit/a85416b7075bed16d5cb29a9827e0a2c1e37119f) — coupling filter as a wavenumber; `dx_hm` to `domain.f90`

See the section above. `setparm` now also prints the grid geometry at startup — `dx_hm`,
`nsx`, host domain width, subdomain width, CRM extent and the host Nyquist wavelength —
which makes a mis-specified configuration visible in the log immediately.

### 3. [`fef8700`](https://github.com/Qiyu-Song/SAM-MMF-Walker/commit/fef87005df5c2f857cda71ecec39272cbf70b034) — `tau_damp_mean` as a namelist parameter, in seconds

`damping_hm` applied a weak drag relaxing the domain-mean wind toward **zero** with a
hard-coded 20-day timescale. In a run with `apply_hm_u_external_nudging = .true.` that
fights the nudging, which is pulling the same quantity toward the prescribed shear.

Now the drag relaxes toward `u_damp_ref`, which is the external profile when the nudging
is on and zero otherwise, so the two terms pull the same way.

| namelist | default | effect at default |
|---|---|---|
| `tau_damp_mean` | `1728000.` (20 days, **seconds**) | same timescale as the old hard-coded value |
| `do_damp_hm_mean` | `.true.` | drag on, as before |

**This is the only commit that changes a number.** Walker runs
(`apply_hm_u_external_nudging = .false.`) are bit-identical. Shear runs differ by the
equilibrium offset the old code produced,
`1/(1 + tauls_large_scale/tau_damp_mean)` = **0.2%** at the default settings — i.e. the
bug it fixes was real but immaterial. Do not expect your shear results to move.

### 4. [`b7d829e`](https://github.com/Qiyu-Song/SAM-MMF-Walker/commit/b7d829e84b33cc6ff229c8deeddd70bcaf1d5d72) — `do_fix_u_halo`

On coupling steps the subdomain U is modified in place by `modify_U_for_subdomain` /
`remove_nyquist_U_for_subdomain` / `remove_residual_U_for_subdomain`, but the MPI halo
is not refreshed, so `advect_mom()` that step sees a stale halo column.

| namelist | default | effect at default |
|---|---|---|
| `do_fix_u_halo` | `.false.` | bit-identical to `tend-nudging2` |

**Tested and it does not matter statistically.** One coupling injects max |dU| = 5.5e-3 m/s;
a within-subdomain composite of `|on - off|` is flat across all 16 cells at every output
time for both USFC and U200, i.e. no seam signature. It changes a trajectory, not the
statistics. Left off by default so old runs stay reproducible.

**Trap if you turn it on**: the refresh must be `periodic(1)`, not `boundaries(1)`.
`main.f90` sets `dompi = .true.` for `dompimmf` before the coupling block, so
`boundaries()` takes the MPI path and deadlocks. See the comment on `do_fix_u_halo` in
`vars.f90`.

### 5. [`b37cf5d`](https://github.com/Qiyu-Song/SAM-MMF-Walker/commit/b37cf5de393f4a2fdd3b087d73a173398e331dcc) — `hm_smoother`: grad^4 and Smagorinsky

The host's horizontal smoother was a grid-**index** Laplacian with no `dx^2`:

    dudt += diffuse_intensity * (u(i+1) - 2u(i) + u(i-1)) / dt_hm_subcycle

so the implied physical diffusivity is `nu = diffuse_intensity * dx_hm^2 / dt_hm_subcycle`
and scales as **dx_hm^2** — 9.1e5 / 2.3e5 / 5.7e4 m2/s at dx_hm = 128 / 64 / 32 km, a 15x
confound across a resolution sweep.

`hm_smoother` selects the operator. The legacy loops are unchanged and became option 1, so
`hm_smoother = 1` with your existing `diffuse_intensity` reproduces previous runs exactly.

| `hm_smoother` | operator | coefficient | implied viscosity |
|---|---|---|---|
| 0 | none | — | — |
| 1 | grad^2 (legacy, default) | `diffuse_intensity` | `I*dx_hm^2/dt_hm_subcycle` |
| 2 | grad^4 | `hyper_intensity` | `H*dx_hm^4/dt_hm_subcycle` |
| 3 | Smagorinsky | `smag_cs` | `(Cs*dx_hm)^2*|du/dx|` |

Other knobs: `smag_max_diff_vel` (default `10.`, WRF's diffusive-velocity cap
`nu <= 10 m/s * dx_hm`) and `smag_nu_max_frac` (default `0.1`, stability clamp
`nu <= frac*dx_hm^2/dt_hm_subcycle`; AB3 needs it below 0.136).

Damping rate of a wave with `k*dx = theta`:
grad^2 gives `(4I/dt_sub) sin^2(theta/2)`, grad^4 gives `(16H/dt_sub) sin^4(theta/2)`.
Both are **exactly dx_hm-independent at 2*dx_hm**, which is the one property of the legacy
operator worth keeping. `hyper_intensity = diffuse_intensity/4` matches the legacy 2*dx
damping exactly (1.25 h at I = 5e-3) while being far weaker at every resolved scale.

`setparm` prints the implied viscosity and the e-folding times at 2*dx and 1000 km at
startup, and warns if a coefficient exceeds the AB3 real-axis limit (0.545).

---

## Recommended settings, and the evidence

```
&KUANG_PARAMS
  suppress_k_start            = -1        ! Nyquist only
  do_remove_coupling_residual = .true.
  hm_smoother                 = 2         ! grad^4
  hyper_intensity             = 1.25e-3   ! at dt_hm_subcycle = 90 s
  diffuse_intensity           = 0.
```

Judged against a pure-SAM benchmark run in the same domain with the same forcing
(`puresam60`, Walker, dx_hm = 64 km, day 30–60), as rms amplitude of U by spectral band
with the ratio to pure SAM in brackets:

| 1003 hPa | 488–148 km | 146–130 km | 2*dx = 128 km |
|---|---|---|---|
| pure SAM | 0.968 | 0.227 | 0.046 |
| no smoother | 0.988 (1.02x) | 0.387 (1.70x) | 0.444 (**9.6x**) |
| grad^2 I = 5e-3 | 0.598 (0.62x) | 0.174 (0.76x) | 0.061 (1.3x) |
| **grad^4 H = 1.25e-3** | **0.740 (0.76x)** | **0.183 (0.80x)** | **0.059 (1.3x)** |

| 418 hPa | 488–148 km | 146–130 km | 2*dx |
|---|---|---|---|
| pure SAM | 0.597 | 0.113 | 0.022 |
| no smoother | 0.858 (1.44x) | 0.346 (3.07x) | 0.076 (3.4x) |
| grad^2 I = 5e-3 | 0.503 (0.84x) | 0.173 (1.53x) | 0.033 (1.5x) |
| **grad^4 H = 1.25e-3** | **0.626 (1.05x)** | 0.186 (1.65x) | 0.035 (1.6x) |

grad^4 reaches the **same 2*dx result** as the legacy grad^2 while leaving the resolved
mesoscale much closer to truth. It also takes the 2*dx mode's persistence from **32.7 h**
(no smoother) to **2.1 h**, which is what the stationary stripes in the lowest-level U
were.

**Caveats, stated honestly.** At 1003 hPa the un-smoothed run is already within 2% of
pure SAM above ~150 km, so the surface case for any smoother rests on the last two or
three wavenumbers; it is aloft that the smoother clearly earns its place. Every smoother
we tried, including this one, leaves the quiescent flanks too quiet (0.65x pure SAM).
Neither grad^2 nor grad^4 fixes the 146–130 km band aloft.

This matches standard practice. Skamarock (2004, MWR 132, 3019) tunes model dissipation
exactly this way — compare the KE spectrum to a reference and judge the tail — and reports
the same grad^2-vs-grad^4 result: matched at 2*dx, the second-order filter "removes a
significant amount of energy at the lower wavenumbers". He also argues that a decaying
tail is the **correct** target, not a failure, because 2*dx–4*dx modes in finite-difference
models have zero or negative group velocity. WRF's effective resolution is ~7*dx.

---

## Things that are not code changes but will cost you time

- **Never run two SAM jobs concurrently from the same model directory.** `<case>/prm` and
  the top-level `CaseName` are shared; the second job's startup overwrites them and both
  processes read the same namelist. This silently happened twice here (jobs 44482857 and
  45705667 ran under the wrong subcase). Only the *startup* window needs protecting —
  `CaseName` is read once in `task_init.f90:54` and `prm` only in the `setparm()` calls at
  `main.f90:45` and `main.f90:80`, all before the time loop.
- **SLURM exit code 9 with a complete log is an MPI-finalize teardown artefact**, not a
  failure. Check for "Finished with SAM, exiting..." and the timing table. Exit 19 is a
  namelist error; exit 24 is a truncated or missing restart file.
- **A job can hang while `squeue` still says RUNNING.** If rank 0 dies (e.g. a bad restart
  file) the allocation is held until the walltime. Watch whether the output file is still
  growing, not whether the job is listed.
- **`dosavemultirestart = .true.` writes one indexed restart set per day.** To branch from
  day N, symlink a new caseid onto index N — see how `shrL64d30` is built in `RESTART/`.
  Pointing `caseid_restart` at the run's own name looks for an *un-indexed* file, creates
  it empty, and dies with `forrtl: severe (24)`.
- **Adding a `use` statement to an existing file silently breaks the parallel build.**
  The Makefile gets its build order from `include Depends`, which is generated under the
  rule `Depends: Srcfiles Filepath` — so it is regenerated only when the *list* of source
  files changes, never when a file's *contents* change. Your `15b4990` added
  `use module_hostmodel` to `forcing.f90` and `setdata.f90`, both existing files, so
  `Depends` stayed stale and `make -j8` compiled `forcing.o` before `module_hostmodel.o`,
  against the previous build's `.mod`:

      forcing.f90(7): error #6580: Name in only-list does not exist or is not accessible.
                                   [SET_UG0_FROM_EXTERNAL_PROFILE]

  **Fix: delete `$SAM_OBJ/Depends` (or the whole `OBJ`) and rebuild.** After a clean
  rebuild the edges are right:

      forcing.o : forcing.f90 simple_ocean.o module_hostmodel.o vars.o params.o microphysics.o
      setdata.o : setdata.f90 vars.o simple_ocean.o module_hostmodel.o params.o microphysics.o sgs.o

  This one failed loudly because the symbol did not exist. The same stale `Depends` can
  fail **silently**: if a module's *data* changes — a default in `vars.f90`, a variable's
  kind — and a dependent file is not recompiled because the edge is missing, you link two
  inconsistent copies of the same module, and it compiles, links and runs. Delete
  `Depends` whenever you add or change a `use`.
- **`nstop` must be at least `nstat`.** `printout.f90:44` has
  `if(nstop-nstep.lt.nstat) call task_abort()`. Setting a short `nstop` for a quick test
  without lowering `nstat` aborts during startup — and the job then **hangs holding the
  whole allocation** until the walltime, while `squeue` still reports it RUNNING.
- **`hm_only = .true.` could not do a wind-shear experiment — your `15b4990` fixes this.**
  On this branch as of `b37cf5d`, the host develops the shear (the external-profile
  nudging is applied inside the subcycle loop, after the `hm_only` branch) but `main.f90`
  takes the `nudging()` path instead of `nudging_hm()`, which relaxes the CRM mean wind
  toward `ug0` from `snd` — whose u column is zero — giving a sheared host over unsheared
  CRMs. Your `set_ug0_from_external_profile`, called every step at the end of `forcing`,
  makes `ug0` *be* the shear profile, so the `nudging()` path now targets the shear too.
  After the merge this trap is gone.

## What we deliberately did not change

- **How the shear profile reaches the model — you have already solved this better than we
  proposed.** On `fix-coupling-residual` the profile is read from `large_u_profile_filename`
  and applied only as a host-side nudging at `tauls_large_scale`; `snd`'s u column is zero,
  so the CRMs never see it except through the host increment. We considered moving the
  profile into `snd` and rejected it, because it would invalidate every existing shear
  `prm` and we measured that it does **not** change convective organization (`shrI_L128`,
  which has the full profile in `snd`, is indistinguishable from `shr_L128`:
  neighbour-column precip correlation 0.896 vs 0.903, phase speed +6.8 vs +7.1 m/s).

  Your `15b4990` takes the better third option: keep the text file as the single source of
  truth, and overwrite `ug0` from it every step in `forcing` plus the CRM initial wind in
  `setdata`. That avoids hand-editing `snd` per experiment, removes the spin-up, and fixes
  the `hm_only` case — none of which the `snd` route would have done as cleanly. Adopt
  yours; we have nothing to add here.
- **`diffuse_TQ` remains commented out** at the `host_model_evolve` call site.
- **No CFL check in the host.** `dx_hm = 32 km` with `hm_subcycle = 5` blew up at day 29 of
  a 60-day Walker run; `hm_subcycle = 10` fixed it. The failure is silent until it happens.
- **`nudging_hm` applies `ug0_hm`, not `ug0`.** In coupled MMF mode `donudging_uv` gates
  the block but what flows through is the host increment, not the sounding. A comment in
  one of our `prm` files claims `snd` is "the only path by which shear reaches the CRMs" —
  that is wrong for coupled runs.
