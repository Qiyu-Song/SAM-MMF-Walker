# Changes to the host model and coupling since `tend-nudging2`

Branch `fix-coupling-residual`, five commits on top of `tend-nudging2` (a760842).
Written for Kairui, 2026-09-14. Everything stays on this branch; nothing is merged to
`main` (`origin/main` is at 2025-09-06 and predates most of the host-model work).

Repo: <https://github.com/Qiyu-Song/SAM-MMF-Walker>

---

## READ THIS FIRST — one change breaks existing `prm` files

[`5b718df`](https://github.com/Qiyu-Song/SAM-MMF-Walker/commit/5b718dfafad3bac517e89aeb98e87964b36a182f)
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

`5b718df` also moved `dx_hm` out of the namelist into `domain.f90` as `dx_hm_km`, next to
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

[`5b718df`](https://github.com/Qiyu-Song/SAM-MMF-Walker/commit/5b718dfafad3bac517e89aeb98e87964b36a182f)
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
`tend-nudging2` exactly. Only `5b718df` requires action, and only `7ad87f6` changes any
number (by 0.2%, in shear runs only).

### 1. [`a6c055f`](https://github.com/Qiyu-Song/SAM-MMF-Walker/commit/a6c055f13e54fff0d473d127c42d8203052977e3) — remove the coupling residual orphaned in subdomain U

This is the fix for the **persistent stripes in the subdomain-mean U**. The mechanism is
worth spelling out, because the fix only makes sense once you see where the stripes come
from, and because the same structure decides which field is affected.

#### The geometry

With `subdomain_center_at_hm_u_center = .true.` (the setting we use), the subdomain
centres coincide with the host **cell centres**, where T and Q live. Host `u` lives on
the **faces** of the C grid. So U — and only U — has to be interpolated in both
directions, while T, Q and the condensates are collocated and pass through untouched:

    t_hm_map = t_hm_map_save + t0_in - t_sub_map_save      ! direct increment, no operator

If you flip that switch to `.false.`, U becomes the collocated field and T/Q become the
interpolated ones, so everything below would apply to T/Q instead.

#### The forward operators are 2-point averages, and they are lossy

    center2face_U:  u_face(i)   = 0.5*( u_center(i) + u_center(i-1) )
    face2center_U:  u_center(i) = 0.5*( u_face(i)   + u_face(i+1)   )

Both have spectral transfer `|cos(theta/2)|` with `theta = pi*m/(nsx/2)`. That is
**exactly zero at theta = pi**, i.e. at `m = nsx/2` — the wave that alternates from one
subdomain to the next. Averaging annihilates it, in both directions.

#### Why an inverse, and why the inverse alone cannot exist

Because plain averaging loses information, the round trip host -> CRM -> host would
smear the field twice. The inherited scheme instead **deconvolves**: given the average
and one endpoint, solve for the other. Both inverses are the same marching recursion:

    face2center_U_inverse:  u_center(i+1) = 2*u_face_filtered(i)   - u_center(i)
    center2face_U_inverse:  u_face(i+1)   = 2*u_center_filtered(i) - u_face(i)

That recursion has a homogeneous solution `x(i) = C*(-1)^i` — the alternating mode —
which it **cannot determine**, because it is precisely the null space of the forward
average. So each inverse projects the Nyquist component out **twice**:

* **before** the recursion, from the *input*: an input carrying that component is
  inconsistent, since no forward average could have produced it;
* **after** the recursion, from the *output*: the marching starts from an arbitrary
  value (the code uses 0) and therefore injects an arbitrary amount of it.

You can see both projections in `face2center_U_inverse` and `center2face_U_inverse` as
the `nyq` / `cnyq` / `fnyq` alternating-sign sums. **This happens in both directions.**

#### The prefilter handles the shoulder, not the null

Just below Nyquist, `1/cos(theta/2)` is finite but large, so the raw inverse amplifies.
`damp_for_target_inverse_prefilter` multiplies by `h_target*cos(theta/2)` *before* the
recursion divides by `cos(theta/2)`:

    prefilter = h_target * cos(0.5*theta)

so the **net transfer of each filtered operator is exactly `h_target`** — a raised
cosine that is 1 below `suppress_k_start`, tapers through a shoulder, and is 0 at
Nyquist. The original setting was `k1 = 0.95*Nyquist`. Both
`face2center_U_inverse_filtered` and `center2face_U_inverse_filtered` call this same
prefilter, so **the filter is applied in both directions too.**

#### Where the stripes come from

Each coupling step the code isolates the part of the change that the CRM generated,
as opposed to the part the host had just handed it:

    call face2center_U_inverse_filtered(u_hm_updated_map_save - u_hm_map_save, tmp1)
    u_adj_cs = u0_in - u_sub_map_save - tmp1        ! CRM-generated, in subdomain space
    call center2face_U_inverse_filtered(u_adj_cs, tmp2)     ! what the host is handed

Sending it to the host costs one factor of `h`; mapping it back costs another. So the
round trip passes `h^2`, and **`1 - h^2` of what the CRM produced never reaches the host
at all**. At `m = nsx/2`, `h = 0`, so `1 - h^2 = 1` — that component is *entirely*
untransferable. Through the shoulder it is partly untransferable.

The host is the only sink for it, so it has none. Every coupling step the CRMs generate
more, it stays in the subdomain-mean U, and it accumulates into a standing pattern whose
structure is exactly the modes where `1 - h^2` is large — alternating from one subdomain
to the next. **That is the striping.**

The earlier mitigation was `diffuse_intensity_subdomain_large_scale`, extra diffusion on
the subdomain-mean fields. It did suppress the stripes, but by damping the symptom, at a
cost to the resolved scales — which is why switching it off is what exposes the free band
cleanly.

#### The fix: measure the orphan and remove it

Rather than damping it, compute exactly what the host could not accept, using the **same
send/return operators the coupling itself uses**, and subtract it from the CRM's U:

    if (do_resid) then
      call face2center_U_inverse_filtered(tmp2, tmp3)   ! what the host effectively accepted
      u_resid_map = u_adj_cs - tmp3                      ! the orphan
    end if

`remove_residual_U_for_subdomain` then subtracts `u_resid` from each subdomain's `u`.
Its spectral weight is exactly `1 - h^2`: identically zero wherever the host can respond,
one where the host is blind. There is **no coefficient and no timescale** — it is a
projection, not a damping, and it removes only what provably cannot be communicated.

Because the weight is built from the coupling's own transfer, it follows the shoulder
automatically. We verified this across three settings without retuning anything: with
`k1 = 0.85` the subdomain-mean U spectrum lifts off the resolved-scale background from
k ~ 69, with `k1 = 0.95` from k ~ 77, and the projection flattens both back down
(`1 - h^2` = 0.27 / 0.75 / 0.98 / 1.00 at k = 77-80 for `k1 = 0.95`, and a pure step at
Nyquist for `k1 = 1.00`).

| namelist | default | effect at default |
|---|---|---|
| `do_remove_coupling_residual` | `.false.` | identical to `tend-nudging2` |
| `do_remove_nyquist_u` | `.false.` | identical — the earlier, narrower version that removes only Nyquist |

**Recommended**: `do_remove_coupling_residual = .true.` with `suppress_k_start = -1`
(Nyquist only). Once the projection is doing the work, the wide shoulder is unnecessary,
and the shoulder costs real signal across its width.

#### One thing this does *not* address

This removes what the CRM produced and the host cannot accept. It says nothing about
what the **host itself** generates at 2*dx — a separate problem, handled by the smoother
in section 5, and diagnosed in `diag/SMOOTHER_TEST_PLAN.md`.

### 2. [`5b718df`](https://github.com/Qiyu-Song/SAM-MMF-Walker/commit/5b718dfafad3bac517e89aeb98e87964b36a182f) — coupling filter as a wavenumber; `dx_hm` to `domain.f90`

See the section above. `setparm` now also prints the grid geometry at startup — `dx_hm`,
`nsx`, host domain width, subdomain width, CRM extent and the host Nyquist wavelength —
which makes a mis-specified configuration visible in the log immediately.

### 3. [`7ad87f6`](https://github.com/Qiyu-Song/SAM-MMF-Walker/commit/7ad87f699191e29ea8e2c6bd277ff05c3586a293) — `tau_damp_mean` as a namelist parameter, in seconds

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

### 4. [`b2d4232`](https://github.com/Qiyu-Song/SAM-MMF-Walker/commit/b2d4232b0f0c9f7998ea1e7fb2b8fb809878747c) — `do_fix_u_halo`

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

### 5. [`77d245d`](https://github.com/Qiyu-Song/SAM-MMF-Walker/commit/77d245da9ed948722184497fa1fbcc9fbffd12b2) — `hm_smoother`: grad^4 and Smagorinsky

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
&PARAMETERS
  dolargescale                = .false.

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
- **`OBJ/Filepath` is written once and never refreshed.** `Build.csh:93` is
  `if ( !(-e Filepath) ) then ... cat >! Filepath`, and the paths it writes are
  absolute. So if you put this code in a *new* `SRC` directory but reuse an existing
  `OBJ/`, the build quietly keeps compiling the **old** tree — no warning, and the exe
  looks fine. `Srcfiles` and `Depends` derive from `Filepath` (`Makefile:92-96`), so
  everything downstream is stale too. Relevant only when you relocate the source; a
  fresh `OBJ/` avoids it.
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
- **How big a difference counts as real.** Two 1-day Walker runs from the same day-30
  restart, differing only at the model top, decorrelate to rms 0.20 K at the surface
  with a mean of 0.005 K (ratio 0.03) — zero-mean, grown from zero. Anything at or
  below that in a 1-day comparison is chaos, not signal.
- **`nstop` must be at least `nstat`.** `printout.f90:44` has
  `if(nstop-nstep.lt.nstat) call task_abort()`. Setting a short `nstop` for a quick test
  without lowering `nstat` aborts during startup — and the job then **hangs holding the
  whole allocation** until the walltime, while `squeue` still reports it RUNNING.
- **`hm_only = .true.` could not do a wind-shear experiment — your `15b4990` fixes this.**
  On this branch as of `77d245d`, the host develops the shear (the external-profile
  nudging is applied inside the subcycle loop, after the `hm_only` branch) but `main.f90`
  takes the `nudging()` path instead of `nudging_hm()`, which relaxes the CRM mean wind
  toward `ug0` from `snd` — whose u column is zero — giving a sheared host over unsheared
  CRMs. Your `set_ug0_from_external_profile`, called every step at the end of `forcing`,
  makes `ug0` *be* the shear profile, so the `nudging()` path now targets the shear too.
  After the merge this trap is gone.

- **`dolargescale = .false.` is the standard from now on.** Historically our Walker
  prms set `.true.` and the shear prms `.false.`; cgmacdonald's pure-SAM runs are
  `.false.` throughout. For our configuration the switch does nothing except at the
  model top: the `lsf` tendencies are zero, and the `ug0` overwrite at
  `forcing.f90:200` is a no-op because the `snd` u column is zero (and after
  `15b4990`, `set_ug0_from_external_profile` at `forcing.f90:344` wins anyway). What
  is left is `upperbound.f90:17` — whether the top one or two levels are relaxed to
  `tg0`/`qg0` on 1 h or left to extrapolate. **Do not flip it mid-run**: an A/B from a
  day-30 restart (`mg_dls` vs `mg_wk`) put the top level 15 K off within 3 h and it was
  still drifting at 0.12 K/h after 24 h, nowhere near settled. Change it at
  initialisation. Walker runs made before this are `.true.` and are not directly
  comparable; that is accepted.
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

  Two small things we noticed while merging, neither of which bites today:

  1. **`set_initial_U_from_external_profile` is called before the `ug` subtraction, not
     after.** `setdata.f90:243` makes the call, but `u0(k) = u0(k) - ug` is at
     `setdata.f90:261`, 18 lines *below* it — contrary to the comment you put at
     `setdata.f90:242` ("上面的 u0 = u0 - ug 已经做完") and the one at
     `module_hostmodel.f90:124`. Because that loop subtracts `ug` from `u0` and `ug0` but
     not from `u_domain_avg`, a non-zero `ug` would leave the three inconsistent:
     `u0` and `ug0` at `profile - ug`, `u_domain_avg` at `profile`. Dormant right now —
     `ug` defaults to 0 (`params.f90:43`) and no `prm` here sets it — but a translating
     frame is a natural thing to want in a shear run, which is exactly this feature's
     use case. The fix is to move the call below the loop that ends at
     `setdata.f90:268`, still before `u(i,j,k) = u0(k)` at `:278`.

  2. **`add_initial_bubble` on a restart taken before host-model init.** The host-side
     call sits inside `if (.not. wsub_inited)` in `host_model_init`, and `wsub_inited`
     is restored from the restart file (`restart.f90:400`), so a normal restart correctly
     does not re-add the bubble. The one gap: restarting from a checkpoint written
     *before* the host model initialised would run `add_initial_bubble_to_hm` against a
     `t0` that never went through `set_initial_bubble_in_crm` (that only runs in
     `setdata`, which restarts skip), so the `- dt_local(k)` would subtract something
     that was never added. Only reachable with `hm_spinup_step > 0`; it is 0 in
     `vars.f90:264` and in every `prm` here.
- **`diffuse_TQ` is gone.** It was a grad^2 smoother for the host `t_hm_map` /
  `q_hm_map`; the subroutine existed but its only call site was commented out, so it
  had never run. Removed rather than left looking like an available option. Three
  reasons:

  1. Not needed. The 2dx problem the `hm_smoother` work addresses is a property of the
     **U staggering**. With `subdomain_center_at_hm_u_center = .true.` (the default in
     `vars.f90:354`, and set explicitly in all 61 `prm` files here) the CRM subdomain
     centres are collocated with the host T/Q cell centres, so T and Q are coupled by a
     plain increment, `t_hm_map = t_hm_map_save + t0_in - t_sub_map_save`, with no
     averaging and no inverse recursion. Only U goes through
     `center2face_U_inverse_filtered` / `face2center_U_inverse_filtered`, whose forward
     transfer `|cos(theta/2)|` vanishes at 2dx and whose inverse carries the `(-1)^i`
     null space. Measured 2dx spike index (m=80 amplitude over the median of m=20..60)
     in `wk_h4a` day 40-60: `T_Out` 0.67-0.74, `Q_Out` 0.46-0.70, against pure SAM
     coarse-grained to the same 160 columns at 0.33-0.43 — and those are lower bounds,
     since a 16-point block average attenuates the coarse-grid Nyquist by
     sinc(pi/2) = 0.637, which puts pure SAM at roughly 0.52-0.68. No clear excess.
     **If anyone ever sets `subdomain_center_at_hm_u_center = .false.`**, T/Q take the
     `else` branch and inherit exactly the same machinery as U, and this question
     reopens.
  2. It was broken as written. `t_map` was `intent(inout)` and updated in place, so the
     loop read the already-updated `i-1` value — not a symmetric Laplacian, and
     direction-dependent. It also modified the field instead of a tendency and did not
     divide by `dt_hm_subcycle`, so its `diffuse_intensity` meant something different
     from the same parameter in `diffuse_u_lap`, and it looped to `nzm` where
     `diffuse_u_lap` stops at `nzm-2`.
  3. It was never wired into the `hm_smoother` dispatcher, and it read
     `diffuse_intensity`, which the recommended configuration (`hm_smoother = 2`) sets
     to 0. Uncommenting it would have been a no-op anyway.
- **The host now has a CFL guard — `hm_cfl_max`, default 0.7.** `kurant_hm` computes
  `sqrt(cflh^2 + cflz^2)` in the convention of SAM's own `kurant.f90`, but with
  `dt_hm_subcycle` and `dx_hm`, once per subcycle. Above the limit it prints the split,
  the level where the vertical term peaks, and the `hm_subcycle` that would be needed,
  and the run aborts. Set `hm_cfl_max <= 0.` to keep the diagnostic and never abort.

  The host has no adaptive subcycling — `hm_subcycle` is fixed, unlike SAM's `ncycle`,
  which `kurant.f90` raises on the fly — so before this the failure was silent until the
  fields blew up. Two things worth knowing:

  - **The constraint is vertical, not horizontal.** Earlier notes here blamed `dx_hm`,
    but `|u|` is tiny against `dx_hm ~ 1e4 m`. Measured over `dxhm32`'s whole life: the
    horizontal term never exceeded 0.033, even in the record where it died. `hm_subcycle`
    is the lever because it shortens `dt_hm_subcycle`, which is what the vertical term
    scales with.
  - **It detects, it does not predict.** `dxhm32` sat in 0.08-0.40 for 28 days, then went
    0.255 (day 28.5) -> 0.723 (day 29.0) and the run ended at 29.15. A healthy
    `dx_hm = 64 km` run (`wk_h4a`, day 30-60) peaks at 0.352. So 0.7 is well clear of the
    healthy range and fires when the run is genuinely failing — but it will not warn you
    a week ahead. The sub-limit warning line (printed on each new record above half the
    limit) is there to show the trend.

  Implementation note if you touch it: `kurant_hm` runs inside `host_model_evolve`, which
  `hm_coupling.f90:176` calls on **masterproc only**. It therefore must not call
  `task_abort` itself — `task_abort` -> `task_stop` -> `MPI_FINALIZE` is collective, so
  aborting from rank 0 alone leaves the others waiting and holds the allocation until the
  walltime. It sets `hm_cfl_abort`; `hm_couple_step` broadcasts that with
  `task_bcast_integer` and every rank aborts together.
- **`nudging_hm` applies `ug0_hm`, not `ug0`.** In coupled MMF mode `donudging_uv` gates
  the block but what flows through is the host increment, not the sounding. A comment in
  one of our `prm` files claims `snd` is "the only path by which shear reaches the CRMs" —
  that is wrong for coupled runs.
