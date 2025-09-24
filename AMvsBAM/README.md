# AM vs BAM (no-HL baseline)

This repo quantifies and explains the separation between AM-BSH (abiotic-only suitability) and BAM-BSH (abiotic + prey requirement) before any habitat loss.

- `src/` provides generation of metawebs, climatic niches, BAM logic, metrics, and Makie plotting helpers.
- `scripts/run_baseline.jl` runs two 2-D sweeps:
  1) (connectance C, alignment)
  2) (diet redundancy R95, niche breadth σ)

Outputs: CSV summaries in `data/` and PNG figures in `data/figs/`.

Run:
```bash
julia --project=. scripts/run_baseline.jl
```