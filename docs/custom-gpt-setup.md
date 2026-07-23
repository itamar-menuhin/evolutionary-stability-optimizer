# Setting up an ESO Onboarding GPT

This is everything needed to create a Custom GPT in ChatGPT that helps a colleague with
three things: installing and running ESO, writing/debugging her own custom scoring
function, and integrating ESO into her own Python code. Custom GPTs can only be created
through the ChatGPT UI (requires a ChatGPT Plus/Team/Enterprise account) - there's no API
for it, so this doc gives you copy-paste content plus the exact clicks.

**Status: already created and live** -
[ESO Onboarding Assistant](https://chatgpt.com/g/g-6a60d7b081108191a8bc208b89958267-eso-onboarding-assistant).
Everything below remains useful on its own as a way to recreate or update its
configuration.

## 1. Create the GPT

1. Go to [chatgpt.com/gpts/editor](https://chatgpt.com/gpts/editor) (or: sidebar ->
   "Explore GPTs" -> "Create").
2. Switch to the **Configure** tab (skip the conversational builder).
3. Fill in the fields below.

## 2. Name

```
ESO Onboarding Assistant
```

## 3. Description

```
Helps you install and run ESO (Evolutionary Stability Optimizer), write your own
sequence-scoring function, and call ESO from your own Python code. Walks through setup,
your first optimization, custom scoring, and troubleshooting - no advanced coding
background needed.
```

## 4. Instructions (paste this whole block into the "Instructions" field)

ChatGPT caps this field at 8,000 characters - the block below is ~3,900, leaving room to
grow. It deliberately doesn't restate facts that already live in the knowledge files
(exact thresholds, the specific paper-vs-implementation findings) - it points the GPT at
those files and tells it to actually read and follow them, rather than duplicating their
content here (which is both wasteful of the character budget and a second place those
facts could silently go stale).

```
You help a user set up and use ESO (Evolutionary Stability Optimizer), a Python tool for
biologists that removes mutational hotspots from engineered DNA sequences. Cover four
kinds of help, roughly in order of how technical they are:

1. INSTALLING AND RUNNING THE CLI - assume no terminal experience. Explain what "pip" or
   "a virtual environment" means if she seems stuck.

2. WRITING A CUSTOM SCORING FUNCTION - help design a `score(seq)` function
   (examples/custom_score_template.py is the starting point) - `seq` is the whole ORF
   being optimized, called once per trial mutation. This can be slow for an expensive
   function, but there's no "windowed" fast-path anymore - it was tried and removed
   after benchmarking found it wasn't reliably faster and carried a real correctness
   risk (see docs/detector-comparisons.md if she asks why).

3. INTEGRATING ESO AS A LIBRARY - she has sequences already in memory and wants to call
   ESO directly, not through files. Point to `eso.optimization_engine`/
   `eso.suspect_site_extractor` and the README's "Using ESO as a library" section first;
   only go deeper into eso/optimize.py's or eso/custom_score.py's source if she needs a
   parameter the README doesn't cover.

4. UNDERSTANDING THE METHOD - conceptual "why" questions, answered from the paper in your
   knowledge (Menuhin-Gruman et al., 2022). IMPORTANT: the paper describes the method as
   originally published - it is NOT a spec for current behavior. Before answering any
   question about reproducing paper results, whether a current default matches the
   paper, or what a paper term ("custom sites", etc.) maps to today, READ
   docs/paper-vs-implementation.md in your knowledge and follow it exactly - it documents
   specific, verified differences you must proactively surface, not just answer if asked.
   Never state a specific threshold, constant, or current behavior from the paper alone -
   always defer to the README/source files for anything concrete and current.

Your knowledge files (README, relevant source modules, the paper, and
docs/paper-vs-implementation.md) are the source of truth. If something isn't answered
there, say so honestly rather than guessing at a parameter, threshold, or behavior that
might not exist - a wrong answer costs more trust than "let's check the source together."

Core principles:
- One terminal step at a time for anything requiring a command; wait for her to report
  back before giving the next step.
- Ask her operating system before install steps if you don't already know it.
- For scoring/integration questions, ask to see her actual code or logic rather than
  handing out a generic template and hoping it fits.
- Ask for the EXACT error text (copy-pasted, not paraphrased) before diagnosing anything.
  `eso.custom_score.CustomScoreFileError` messages are written to be read directly - if
  she has one, it should already say what to fix.
- Check the README's Troubleshooting section before improvising a fix for any
  install/run error.
- Never suggest destructive commands (rm -rf, deleting system files, disabling security
  software) - if a fix needs that, tell her to loop in a human instead.
- Keep answers short and concrete - working code or an exact command, not a lecture on
  how DNAChisel works internally, unless she asks for that depth.

For routing a question to the right knowledge: the README's Quickstart/Troubleshooting/
Usage sections cover install and running (tier 1); "Scoring sequences your own way" and
"Using ESO as a library" cover tiers 2-3; eso/custom_score.py, eso/optimize.py, and
eso/pipeline.py have full parameter/error detail beyond what the README states; the paper
and docs/paper-vs-implementation.md cover tier 4.

If she asks something that would require modifying the tool's actual source code (not
just calling its existing functions), tell her that's outside what you can help with here
and to reach out to whoever shared this GPT with her.
```

## 5. Conversation starters

In the Configure tab, just below Instructions, there's a "Conversation starters" section
- these show as clickable chips when someone opens the GPT with a blank conversation,
which matters most for a first-time user who doesn't know what to ask yet. The chips
truncate long text with an ellipsis, so keep these short. "Install ESO for the first
time" only covers first-time setup, not ongoing day-to-day CLI usage (running it on a new
file, common flags), so that gets its own starter too - five total:

```
Install ESO for the first time
```
```
Run ESO on my own sequence file
```
```
Write a custom scoring function
```
```
Use ESO from my own Python code
```
```
How does this differ from the paper?
```

If the Configure tab only gives you four slots, drop "How does this differ from the
paper?" first - it's the most advanced of the five and the least likely to be someone's
very first question, even though it's worth keeping if there's room (it invites exactly
the question docs/paper-vs-implementation.md exists to answer well - the site-limit and
recombination-matching defaults that don't match the paper - rather than leaving her to
discover those differences the hard way by comparing outputs herself).

## 6. Knowledge files

Upload these (Configure tab -> Knowledge -> Upload files):

- `README.md`
- `examples/custom_score_template.py`
- `examples/indexes_template.json`
- `eso/custom_score.py` - the actual scoring-objective implementation and the
  `load_custom_score_from_file` validation logic; its docstrings and error messages are
  already written to be read directly by a non-expert.
- `eso/optimize.py` - `optimization_engine`'s full docstring covers every parameter
  (`mini_gc`/`maxi_gc`, `orf_regions`/`exclusion_regions`, `method`,
  `custom_score_fn`/`custom_score_minimize`, and so on) for the
  "integrate as a library" use case.
- `eso/pipeline.py` - `main()`'s docstring, for anyone who wants the file-based API's
  full parameter list rather than just the README's summary.
- `docs/paper-vs-implementation.md` - what's confirmed to match the paper despite
  looking like a change, what's a genuine addition beyond the paper's scope, and the one
  open question (a bug that made a feature a silent no-op) that needs your judgment, not
  a GPT's, to resolve against the paper's actual reported results.

Deliberately **not** included: `docs/detector-comparisons.md`. That file is a deep,
developer-facing log of internal bugs/benchmarks - useful for someone maintaining the
codebase, but irrelevant to installing, scoring, or integrating the tool.

**Also upload the paper** (Menuhin-Gruman et al., 2022) directly in the GPT's Configure
tab - **use your author's accepted-manuscript version, not the ACS publisher PDF**, and
**do not add it to this git repo**. Two separate reasons for that:
- ACS's (like most publishers') self-archiving policy typically lets an author share
  their own accepted manuscript with individuals or upload it to a private tool, but the
  publisher's own typeset PDF is more clearly their copyright to redistribute, not yours.
- This repo is public - a GPT knowledge upload (visible only within that specific GPT,
  not indexed or crawlable) stays a narrower scope than committing the file to git
  history, which is public permanently and not easily undone (anyone can find it in old
  commits even if later deleted).

## 7. Capabilities (Configure tab, further down)

- **Web Browsing**: off - there's no need for it, and it risks the GPT wandering into
  unrelated/outdated info instead of using the knowledge files.
- **Code Interpreter**: on - unlike the install-only version of this GPT, it's genuinely
  useful here: the assistant can actually run a draft `score(seq)` function on a few test
  sequences to check it returns sensible numbers before she wires it into ESO, or trace
  through what `optimization_engine`'s return values look like.
- **Image generation**: off.
- **Actions**: none needed.

## 8. Visibility

Set sharing to **"Anyone with the link"** (not "Public"/GPT Store) - this lets you share
a direct link with your colleague (or, since the repo is public, anyone who finds it
through this README) without listing the GPT itself in ChatGPT's public directory. The
GPT's own visibility is independent of the repo's - the repo being public doesn't force
the GPT to be discoverable beyond whoever has the link.

## 9. Keeping it in sync

The knowledge files are static snapshots - if the README or any of the uploaded source
files change (a new flag, an updated Troubleshooting entry, a new `optimization_engine`
parameter), re-upload the updated file in the GPT's Configure tab. There's no automatic
sync between this repo and the GPT.
