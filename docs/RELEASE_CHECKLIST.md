# Release checklist (v0.2.0 and later)

A stable, citable release is required by the reviewers (version number, Git tag, archived snapshot
with a DOI, DOI cited in the manuscript). Work through the list in order.

## 1. Before tagging

- [ ] All analyses used in the manuscript were produced with the code at the commit that will be tagged
      (`docs/revision/RERUN_GUIDE.md`); the numbers in the manuscript and the response letter were filled in
      from the re-run output (no `[TODO: ...]` placeholder left).
- [ ] `pytest` passes (`pip install -e ".[dev]" && pytest`).
- [ ] `scripts/` and `src/gloom/pipeline/` are identical (the packaged copy used by the CLI):
      `diff -rq scripts src/gloom/pipeline` should only list the helper tools that exist in `scripts/` only
      (`fetch_gdc_tcga_luad.py`, `benchmark_runtime.py`).
- [ ] Version is `0.2.0` in `pyproject.toml`, `src/gloom/__init__.py` and `CITATION.cff`
      (`conda.recipe/meta.yaml` still points at the v0.1.0 archive and its sha256; update it only after the tag exists).
- [ ] `CHANGELOG.md` has the release entry with today's date; remove "Unreleased" items that shipped.
- [ ] README: installation instructions verified on a clean machine (`git lfs install`, `git lfs pull`,
      `mamba env create -f environment.yml`, `pip install -e .`, `gloom --version`).
- [ ] No claim of Bioconda availability anywhere (README, manuscript) unless the package is really published.
- [ ] `git lfs ls-files` lists the data files and they are pushed (`git lfs push --all origin`).

## 2. Tag the release

```bash
git checkout main && git pull
git tag -a v0.2.0 -m "GLOOM 0.2.0: revision release"
git push origin v0.2.0
```

Create a GitHub Release from the tag (Releases -> Draft a new release -> choose `v0.2.0`), with the
CHANGELOG entry as the description.

## 3. Archive on Zenodo (GitHub-Zenodo integration)

1. Log in to <https://zenodo.org> with the GitHub account that owns (or can administer) the repository.
2. Zenodo -> *Account* -> *GitHub*: switch the toggle **on** for the `omicscodeathon/gloom` repository
   (do this **before** publishing the GitHub Release; an organization repository needs an administrator to enable it).
3. Publish the GitHub Release created above. Zenodo archives it automatically and mints a DOI
   (a *version DOI* for this release and a *concept DOI* that always points to the latest version).
4. Open the Zenodo record, check title / authors / license (MIT) / description; edit the metadata if needed
   (new version of the record is **not** required for metadata fixes).
5. Copy the DOI.

## 4. After archiving

- [ ] Put the DOI in `CITATION.cff` (replace `10.5281/zenodo.XXXXXXX`) and in the README (badge / citation line); commit.
      (Tag again only if the archived code itself changed; metadata-only commits do not need a new release.)
- [ ] Cite the DOI in the manuscript (Data/Code availability section and reference list) and in the response to reviewers.
- [ ] Optional: add the Zenodo DOI badge to the README:
      `[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)`.
- [ ] If a Bioconda recipe is submitted later, update `conda.recipe/meta.yaml` (source URL + sha256 of the tagged
      release archive) and only then restore a Bioconda installation line in the README.
