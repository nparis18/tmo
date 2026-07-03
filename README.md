# tmo — Thresholding Multiple Outcomes in Stata

Stata implementation of **Thresholding Multiple Outcomes (TMO)**, the procedure for
inference with spatial dependence proposed by DellaVigna, Imbens, Kim and Ritzwoller
(2025, [arXiv:2504.13295](https://arxiv.org/abs/2504.13295)), together with the
companion Stata Journal article *"A Practitioner's Guide to Inference with Spatial
Dependence"* (DellaVigna, Imbens, Kim, Ritzwoller, Clarke and Parı́s).

## Repository structure

```
├── src/                 Stata package (what gets installed)
│   ├── tmo.ado          Current tmo command
│   ├── tmo_new.ado      Development version (tmo_NP): plots, SCPC/Conley combos
│   ├── tmo.sthlp        Help file
│   ├── tmo.pkg          Package file
│   └── stata.toc        Table of contents for net install
├── example/             Example datasets (county_differences.dta, county_panel.dta)
├── scpc_tmo/            Müller–Watson SCPC package (dependency, vendored copy)
├── SJ/                  Stata Journal paper
│   ├── paper/           spatial.tex + style files, compiled spatial.pdf
│   │   ├── examples/    sjlog output (*.tex.log) input by spatial.tex
│   │   └── figures/     All figures included by spatial.tex
│   ├── code/            codeSj.do (example logs + tmo figures), maps.do (map figures)
│   ├── data/maps/       Census cartographic boundary files for maps.do
│   └── temp/            Scratch space (test.do, logs) — not part of the build
└── supporting_material/ Reference papers and replication material
```

## Reproducing the paper

1. Run `SJ/code/codeSj.do` — regenerates every example log in `SJ/paper/examples/`
   and the figures `_hist.png` / `_qt.pdf` in `SJ/paper/figures/`.
2. Run `SJ/code/maps.do` — regenerates the method-comparison maps
   (`state_clustering.pdf`, `Conley.pdf`, `SCPC.pdf`, `tmo.pdf`) in `SJ/paper/figures/`.
3. Compile `SJ/paper/spatial.tex` (pdflatex + bibtex).

All figures used by the paper live in a single folder, `SJ/paper/figures/`; do-files
write there directly.

## Requirements

Stata 16+. User-written commands used by the examples: `reghdfe`, `ivreg2`,
`ivreghdfe`, `geoplot`/`geoframe` (maps), and `scpc` (included in `scpc_tmo/`).

## Related repositories

- Original TMO code: <https://github.com/wjnkim/tmo>
- SCPC (Müller–Watson 2022): <https://www.princeton.edu/~umueller/>
