# SiCell.jl repo updates — what to do with each file

Read this whole thing before touching git. Three of the items below are files
you can drop straight in; two are things only you can do (I don't have push
access to your repo, and one requires a judgment call I shouldn't make for
you).

---

## 1. Drop-in files (safe to copy over as-is)

**`CITATION.cff`** → repo root, overwrite the existing one.
Complete rewrite: proper title/version/date, your ORCID
(0000-0002-8602-8092), MIT license, a `preferred-citation` block for the
software paper, and a `references` entry for the TUF preprint (DOI already
filled in). One thing left for you: once the software paper is on bioRxiv,
add a `doi:` line under `preferred-citation` (marked with a `# TODO` comment
in the file so you don't lose track of it).

**`test/runtests.jl`** → `test/`, overwrite the existing one.
Your original file is untouched except for one insertion: a new `@testset
"Trajectory Uncertainty (TES/TDS/LPS)"` block, slotted in as item 6 (which
was oddly missing between your 5 and 7 — this fills that gap). It builds a
small synthetic object, runs the full pipeline through `trajectory_uncertainty!`,
and checks: the three output columns exist, TES and TDS are bounded in
[0,1] as documented, and the pseudotime-out-of-range guard actually throws.
Verified byte-for-byte CRLF-consistent with your original file, so `git diff`
will show a clean, minimal, single-block insertion — not a wall of red from
line-ending noise.

**`benchmarks/`** → new folder at repo root.
Contains `synthetic_scaling_benchmark.jl` (the final, working version) and
both result CSVs. This is what the paper's Data and Code Availability
section now points to — it needs to actually exist in the repo before you
submit anywhere.

## 2. Requires your judgment call — the TDS question

I have **not** touched `trajectory_uncertainty!` in `src/trajectory_inference.jl`.
The shipped code computes a forward-only TDS; the TUF preprint you already
submitted describes an omnidirectional one. The software paper now explains
this as a deliberate refinement — but that's true only if it actually *was*
deliberate. Two honest paths, and only you know which one is right:

- **If forward-only was an intentional improvement after the TUF preprint
  was written:** leave the code as-is. The paper's current wording already
  covers this correctly.
- **If it wasn't intentional — e.g., it drifted during later development
  without you deciding TUF's own math should change:** the code should
  match the published, validated formula (omnidirectional — average over
  *all* neighbors, not just forward-pseudotime ones), and the paper's
  explanatory sentence should be removed since it would no longer apply.

I didn't guess here because it's a scientific claim about what TUF actually
measures, not a bug fix — that call is yours. Tell me which way and I'll
write the actual code change.

## 3. Requires GitHub directly — I can't do this for you

**Delete the malformed `"0.1.1"` tag** (the one with literal quote characters
in its name, confirmed via `git ls-remote --tags` — sitting alongside the
correct `v0.1.1`). Easiest: GitHub web UI → repo → Releases/Tags → find the
entry literally showing `"0.1.1"` with quotes → delete it. Leave `v0.1.1`
alone. Terminal alternative (PowerShell or Git Bash — the single quotes
matter, they keep the double quotes literal):
```
git push origin --delete 'refs/tags/"0.1.1"'
```

## 4. Your local Project.toml — read this before editing

I do **not** have your actual current local `Project.toml` — only the
`[compat]` block you pasted in chat. I'm not handing you a full-file
replacement, because I can't see your local `[deps]` UUIDs for
`CommunityDetection`, `igraph_jll`, `ProfileView`, or `BenchmarkTools`, and
guessing at those risks silently breaking something I can't see.

`reference/Project.toml.github-reference` in this package is your **current
clean GitHub version**, included only so you can diff against it — do not
copy it over your local file, it's missing whatever you've been building
locally.

What to actually do, in your real local `Project.toml`:

1. Before deleting anything, run `grep -r "CommunityDetection" src/` in your
   local checkout. I confirmed it's unused in the GitHub version, but I've
   never seen your local `src/` changes — if you've started wiring it into
   something new, don't delete it blindly.
2. **If unused:** delete the `CommunityDetection = "aaaa..."` line from
   `[deps]` and the `CommunityDetection = "0.2.0"` line from `[compat]`.
   Two lines, nothing else.
3. **If it is used:** instead loosen `Clustering = "0.15"` to
   `Clustering = "0.14, 0.15"` in `[compat]` — but first check nothing you
   wrote depends on something added in Clustering 0.15 specifically.
4. Either way, sanity-check the result before pushing:
   ```julia
   using Pkg; Pkg.activate(temp=true); Pkg.develop(path="."); Pkg.resolve()
   ```
   If that resolves cleanly, you're done. This is the same check worth
   running before any future push, not just this once.

## 5. Suggested commit

Once 1–4 are done:
```
git add CITATION.cff test/runtests.jl benchmarks/ Project.toml
git commit -m "Add CITATION.cff, TUF test coverage, synthetic benchmark suite; fix Project.toml compat conflict"
git push
```
Then delete the malformed tag (step 3) separately — it's not a file change,
so it won't be part of that commit.
