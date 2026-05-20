# How to contribute to this repository

This document explains how Advik, Arockiaraj sir, and any future
collaborator can work on the project alongside the maintainer
(`Singati2` on GitHub).

## Option A — push access (recommended for active collaborators)

1. Send your GitHub handle to the maintainer.
2. The maintainer adds you with `gh api -X PUT
   repos/Singati2/topological-index-orthogonality/collaborators/<handle>
   -f permission=push`. You'll receive an email invitation; accept it.
3. Clone:
   ```
   git clone https://github.com/Singati2/topological-index-orthogonality.git
   cd topological-index-orthogonality
   ```
4. Create a feature branch for your work:
   ```
   git checkout -b paper1-section3-bounds
   ```
5. Make changes (edit files, add new scripts, etc.). Use commits in
   small logical units; meaningful commit messages help reviewers.
   ```
   git add docs/paper1_wuzi_manuscript_skeleton.md
   git commit -m "Section 3: add Theorem 3.1 with proof"
   ```
6. Push the branch and open a pull request:
   ```
   git push -u origin paper1-section3-bounds
   gh pr create --base main --fill
   ```
7. Maintainer reviews and merges, or requests changes.

## Option B — pull requests from a fork (lighter-weight)

If you don't want push access:

1. Fork the repository on GitHub (use the **Fork** button).
2. Clone *your* fork, make changes on a branch, push.
3. Open a pull request from your fork to
   `Singati2:topological-index-orthogonality:main`.

## Option C — share content directly (no Git)

If Git is friction, you can send LaTeX / Markdown / images directly to
the maintainer. The maintainer will commit them with appropriate
attribution in the commit message.

## Repository conventions

- **Branch names:** `paper1-section<N>-<topic>` for Paper 1 (Wuzi)
  sections; `paper2-<topic>` for Paper 2 (orthogonality screening)
  work; `fix-<thing>` for bug fixes.
- **No raw data in git:** datasets are auto-downloaded by `src/load_data.py`.
  Don't commit `data/`.
- **Reproducibility:** every new script should run end-to-end from a
  fresh checkout with only `pip install -r requirements.txt`. If your
  script depends on a script that runs first, document the order in
  `README.md`.
- **Citations:** placeholder citations go to `docs/literature_notes.md`
  with `[REFERENCE NEEDED]` markers. Filling them requires reading the
  actual paper, not guessing.
- **Section labels in manuscript skeleton:** sections marked
  `[NEEDS DERIVATION BY ...]` are explicit gaps; if you fill one, remove
  the marker and add the content.

## Maintainer

Ganesh Shiwakoti (`mpcrlab@gmail.com`).
