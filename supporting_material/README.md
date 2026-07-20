# supporting_material

The replication packages and article supplements used by this project live
once, outside the repo, in `~/Dropbox/literature/` (one folder per paper), so
other projects can reuse them without duplicating data. The entries in this
folder are local symlinks into that location and are not tracked by git.

The small final datasets that SJ/code/codeSj.do actually loads are committed
as real files under SJ/data/replication/, so the paper replicates from a
clone of this repo alone:
  - dataset_wide_1.dta        (Bernini et al. 2023 package)
  - Republican_vote_data.dta  (Bazzi et al. 2023; rebuilt from the package's
    raw data -- see README_generated.txt alongside it)

Packages currently in ~/Dropbox/literature:
  acemoglu et al (2019)  bazzi et al (2023)  bernini et al (2023)
  muller watson (2022)   muller watson (2024)
