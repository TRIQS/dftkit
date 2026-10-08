# Charge self-consistent DFT+DMFT for SrVO3 with Wien2k

A charge self-consistent DFT+DMFT calculation using Wien2k and ModEST.

Here, we project onto the V-t2g shell in a wide energy window (-10 eV to +10 eV) around the Fermi energy.

## What is here

| file | role |
|---|---|
| `SrVO3.struct` | the Wien2k structure file |
| `SrVO3.indmftpr` | the dmftproj projector definition |
| `wien2k_modest_csc.py` | the CSC DFT+DMFT script |

## Prerequisites

* Wien2k, with `WIENROOT` set and `$WIENROOT/x` on disk
* `dmftproj` on `$PATH` (built and installed by this package)
* `triqs_modest`, and an impurity solver — the script uses `triqs_cthyb`

## Running it

Wien2k takes the case name from the directory name, so the directory must be
called `SrVO3` to match `seedname` in the script.

```bash
mkdir SrVO3 && cd SrVO3
cp /path/to/wien2k_csc_svo/SrVO3.struct .

# generate the Wien2k inputs (in0, in1, in2, inm, inc, klist, the starting
# density, ...). Accept the defaults unless you know you want otherwise.
init_lapw -b -numk 1000 -rkmax 7.0

# the dmftproj input: either copy the one provided
cp /path/to/wien2k_csc_svo/SrVO3.indmftpr .
# ...or generate your own interactively
init_dmftpr

cp /path/to/wien2k_csc_svo/wien2k_modest_csc.py .
mpirun -n 16 python wien2k_modest_csc.py
```

The driver converges the DFT SCF cycle itself on the first call, so there is no
separate `run_lapw` step. An already-converged `case.scf` is reused; pass
`force_scf=True` to `one_body_elements_from_dft()` to redo it.

The run is checkpointed every iteration into a `svo_csc_beta<beta>_U<U>_J<J>.ckpt`
directory, so it can be killed and restarted: the DMFT state is restored from the
checkpoint and the Wien2k state from the files on disk. Restarting is just
re-running the same command.

To check the DFT part alone first, without any DMFT:

```python
from triqs_dftkit.wien2k import Driver
Driver("SrVO3").run_dft_only()          # or run_dft_only(n_iter=1) for one cycle
```

## The projection window

Line 15 of `SrVO3.indmftpr` is the correlated energy window relative to E_F, in
Rydberg.

Every k-point must have bands inside this window: `lapw2 -qdmft` indexes its
density matrices per k-point but only sets their dimension for k-points that
dmftproj marked as included, so an excluded k-point would be read against
uninitialised memory. The driver checks `case.oubwin` and refuses up front
rather than letting that happen, so if it reports excluded k-points, widen the
window.

## Notes

`lapw2 -qdmft` runs serially regardless of how the rest of the calculation is
parallelised. The impurity solver still uses all available ranks.

The entire Wien2k chain runs on the master rank alone. Typically, the job size
is determined by the impurity solver, as we expect the DFT part to take a small
wall time.

The driver takes a `verbosity` argument: `1` (the default) prints one line per SCF
cycle, `2` adds one line per Wien2k program launched, and `0` leaves only warnings.
Warnings ignore the setting and go to stderr, so quietening the progress output
cannot hide a problem.

Only non-magnetic calculations (`SP=0`, `SO=0`) are supported at present.
