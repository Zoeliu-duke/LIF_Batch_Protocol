# Batch image analysis protocol — Leica `.lif` confocal immunofluorescence

**Pipeline files:** `process_lif.py` (Python, does everything) · `fiji_macro_template.ijm`
(template the script turns into a ready-to-run Fiji macro)
**Works for:** any Leica `.lif` file with any number of channels, any antibodies and colours,
tile scans and/or single fields, any objectives.
**Worked example used throughout:** `ZL01_20260925.lif` — mouse lymph node, CD8a titration,
frozen sections; C1 DAPI · C2 CD8a-AF594 · C3 CD3-FITC · C4 CD4-AF647; 10x tile scans + 63x fields.

The pipeline turns one `.lif` file into figures and analysis-ready TIFFs with **identical
processing and identical display settings for every sample of the same magnification**, so
samples can be compared channel by channel.

---

## Contents

1. [What the pipeline does](#1-what-the-pipeline-does)
2. [Before you start — acquisition checklist](#2-before-you-start--acquisition-checklist)
3. [One-time setup](#3-one-time-setup)
4. [Protocol A — a new file with the same antibody panel](#4-protocol-a--a-new-file-with-the-same-antibody-panel)
5. [Protocol B — a different antibody panel / channels / colours](#5-protocol-b--a-different-antibody-panel--channels--colours)
6. [Processing several `.lif` files (a batch)](#6-processing-several-lif-files-a-batch)
   - [Running again: overwrite or keep](#running-again-overwrite-or-keep-the-previous-results)
   - [Merged z-stack instead of one plane](#merged-z-stack-projection-instead-of-one-plane)
7. [Quality control checklist](#7-quality-control-checklist)
8. [Fiji: the generated macro and doing it by hand](#8-fiji-the-generated-macro-and-doing-it-by-hand)
9. [Output files reference](#9-output-files-reference)
10. [Settings reference](#10-settings-reference)
11. [Notes on the ZL01 example](#11-notes-on-the-zl01-example)
12. [Fix the image — what to change](#12-fix-the-image--what-to-change-when-the-figure-does-not-look-right)
13. [Troubleshooting](#13-troubleshooting)
14. [Methods text (template)](#14-methods-text-template)
15. [Glossary](#15-glossary)

---

## 1. What the pipeline does

For every image in the `.lif` file, in this order:

| Step | What happens | Why |
|---|---|---|
| 1. Sort | Images are grouped by objective (read from the file, e.g. `10x`, `63x`). Tile scans use Leica's own `…_Merging` image by default (`STITCH` method `"leica"`). Samples matching `EXCLUDE`, or with a different number of channels than `CHANNELS`, are skipped. | Each group gets its own display settings. |
| 2. Stitch (tile scans only) | Default: Leica's merged image. With method `"registered"`: raw tiles are stitched from scratch (start at 10 % overlap, measure the real offset of every neighbouring tile pair, place the tiles, blend seams linearly). | Check `stitching_review/` (one sheet per scan) and pick the better option — for all scans (`STITCH` method) or scan by scan (`STITCH_OVERRIDE`). ZL01 (10x LN): stitching here was better (Leica's merge misaligned 5 of 7 scans). FE71 (63x brain, sparse seams): Leica's merge was better. Registered = Fiji *Grid/Collection stitching* with "Compute overlap". |
| 3. Z | Default: keep the **sharpest single plane** (focus score = edge strength of the `FOCUS_CHANNELS`). Alternatives: average or max projection. | Thick sections with widely spaced planes make projections muddy. |
| 4. Denoise | Each channel: 3×3 median filter, then Gaussian blur σ = 0.8 px. | Removes photon "grain" without blurring cell outlines. |
| 5. Display ranges | Per channel, one min/max pair is calculated from the group's **STANDARD** sample (plus background rules) and applied unchanged to every sample in the group. | Samples are directly comparable; background haze is shown as black. |
| 6. Save | Panel figures, merges, TIFFs (raw + denoised + full z-stack), comparison sheets, tables, a **ready-to-run Fiji macro**, and check figures. | See [section 9](#9-output-files-reference). |

Nothing is auto-contrasted per image. Pixel values in the `_raw` TIFFs are never changed.

---

## 2. Before you start — acquisition checklist

Processing cannot fix acquisition differences. For samples to be comparable:

- [ ] **Same microscope settings for every sample you want to compare:** laser power, detector
      gain, zoom, pinhole, scan speed, bit depth, averaging/accumulation, z-step.
- [ ] **Same channel order** in every sample of the file.
- [ ] **Consistent, informative sample names** (e.g. `LN02_CD8_1-50_10X`). Controls share one
      word (e.g. `noPrimary`, `FITCsec`).
- [ ] **Decide the STANDARD sample** for each objective — a well-stained, representative sample.

General recommendations (learned from ZL01):

- **12-bit rather than 8-bit** — more room for dim signals, fewer saturated pixels.
- **Dim channels:** raise that laser and add **line accumulation ×2–4**; check saturation with
  the QLUT/glow lookup table. **Line average ×2** reduces grain in all channels.
- **Thinner sections / one well-focused plane** instead of a few widely spaced planes.
- **Flat-field reference:** image a uniform fluorescent slide once per session per channel with
  the same settings (removes the faint tile grid in stitched scans; not yet built into the script).
- **Clean controls:** secondary-only or isotype controls in an animal without fluorescent
  reporters in that channel.

---

## 3. One-time setup

### 3.1 On this Mac (already done)

The Python environment with all image libraries is at `~/.venvs/imaging`. Check it by opening
**Terminal** and pasting:

```bash
~/.venvs/imaging/bin/python -c "import readlif, skimage, tifffile; print('ready')"
```

It should print `ready`.

> **Important:** the normal Python/IDLE on the Mac does **not** have these libraries. Always use
> `~/.venvs/imaging/bin/python`, or open IDLE with it (Protocol A, step 3).

### 3.2 On a different computer

1. Install Python 3.12 from python.org (Windows: tick "Add python.exe to PATH").
2. In Terminal (Mac) or Command Prompt (Windows), create the environment:

```bash
python3 -m venv ~/.venvs/imaging
```

3. Install the libraries:

```bash
~/.venvs/imaging/bin/pip install readlif numpy scipy scikit-image tifffile matplotlib pandas pillow PyWavelets
```

Versions used for ZL01: Python 3.12.8, readlif 0.6.6, numpy 2.5.3, scipy 1.18.1, scikit-image
0.26.0, tifffile 2026.9.20, matplotlib 3.11.2, pandas 3.0.6, pillow 12.3.0, PyWavelets 1.10.0.
On Windows use `%USERPROFILE%\.venvs\imaging\Scripts\python.exe` wherever this README says
`~/.venvs/imaging/bin/python`.

### 3.3 Fiji (optional)

Standard Fiji; no extra plugins. Only needed for [section 8](#8-fiji-the-generated-macro-and-doing-it-by-hand).

---

## 4. Protocol A — a new file with the same antibody panel

Time: ~10 min hands-on + ~2 min computing per ~1 GB `.lif`.

### Step 1 — Create a project folder

1. In Finder, create a folder, e.g. `Documents/Claude code/ZL02`.
2. Copy the new `.lif` into it. **One `.lif` per folder** (outputs would otherwise overwrite each other).
3. Copy into the same folder: `process_lif.py`, `fiji_macro_template.ijm`, and this `README.md`.

```
ZL02/
├── ZL02_2026xxxx.lif
├── process_lif.py
├── fiji_macro_template.ijm
└── README.md
```

### Step 2 — Write down your samples

In your notebook: the **STANDARD** sample for each objective, the word that marks **controls**
(or "none"), and samples to **exclude** (test scans, failed sections).

### Step 3 — Open the script in IDLE (with the imaging Python)

**Do not use the IDLE app from Applications / the Dock** — that is a plain Python without the
imaging packages (`No module named 'numpy'`). Easiest: double-click **`Open IDLE for imaging.command`**
in the `image` folder (the first time, macOS may ask: right-click → Open → Open). Or:

1. Open **Terminal** (Cmd + Space → `Terminal` → Enter).
2. Paste and press Enter:

```bash
~/.venvs/imaging/bin/python -m idlelib
```

3. An **IDLE Shell** window opens (leave Terminal open in the background).
4. **File › Open…** → your project folder → `process_lif.py`. A code window opens; edit there.

> **Paste into a new, empty window (File › New File), not the "IDLE Shell" window** that opens
> first (the one showing `Python 3.12…` and `>>>`). A file saved from the Shell contains that text and
> `...` prompts, and running it gives **invalid syntax** on line 1.

### Step 4 — Edit the SETTINGS block

Everything you normally change is between `SETTINGS` and `END OF SETTINGS` at the top.
For the same panel, change only:

| Setting (section) | What to type | Example |
|---|---|---|
| `SRC` (1) | Exact `.lif` file name, in quotes. | `SRC = "ZL02_20261015.lif"` |
| `STANDARD` (5) | Leave it — the first run prints the exact names. | — |
| `CONTROL_TAG` (5) | Word shared by control names, or `""`. | `CONTROL_TAG = "FITCsec"` |
| `EXCLUDE` (5) | Words identifying samples to skip, or `[]`. | `EXCLUDE = ["test"]` |

**Editing rules:** keep quotes `" "`, commas and brackets exactly as they are; only change text
inside quotes or numbers. Save with **Cmd + S**.

### Step 5 — First run (gets the sample names)

1. Click in the code window → **Run › Run Module** (or **fn+F5**; on a Mac keyboard F5 alone is the
   Dictation key — if a "Do you want to enable Dictation?" box appears, click **Not Now**).
2. The IDLE Shell shows `loaded 10x …` lines (tile scans are slow: they are being stitched).
3. The run stops with:

```
SETTINGS PROBLEM: STANDARD for group '10x' is 'ZL01_LN02_CD8_1-50_10X', but the samples in this group are:
  ZL02_LN01_CD8_1-100_10X
  ZL02_LN01_CD8_1-50_10X
Copy one of these names into STANDARD, e.g.
STANDARD = {'10x': 'ZL02_LN01_CD8_1-100_10X', '63x': 'ZL02_LN01__CD8_1-100_63X'}
```

4. Copy the last line into the script, then change each name to the sample you chose as
   standard (from the list). Copy names exactly — tile-scan names have single underscores,
   single-field names keep the original ones. Save (**Cmd + S**).

### Step 6 — Full run

Run it again (**Run › Run Module** or **fn+F5**). Expected messages:

```
loaded 10x ...                  (one per sample)
saved 10x ...                   (one per sample)
wrote processed/fiji_macro.ijm
extra: extra_coloc_check
extra: extra_stitching_review
extra: extra_display_style
done
```

`done` = finished. `extra skipped: …` = one check figure could not be made (main outputs are
complete). A `SETTINGS PROBLEM` or red error = stopped; see [section 13](#13-troubleshooting).

### Step 7 — Quality control

Go through the [QC checklist](#7-quality-control-checklist) **before using any figure**.

### Step 8 — Record (the web page can write the lab note)

In your notebook: `.lif` name, date, the whole SETTINGS block (copy-paste it), and the display
ranges (`processed/display_ranges.json`, also printed at the bottom of each `_panels.png`). Keep
the exact `process_lif.py` you ran next to the outputs. The LIF Batch Protocol web page has a
**Lab note** section (bottom of the left column): fill in experiment name, date, purpose, animals
(genotype, sex, age, treatment), tissue, and **panels** — each staining condition with its own
antibodies (clone, fluorophore, dilution), linked samples, and optionally its own animals and
organ ("Same for all panels" switches) — plus results and summary, choose **First analysis** or **Re-analysis** (adds previous analysis, reason, what
changed), attach `processed/run_settings.json` as the **Run record**, and press **Generate lab
note** — it copies a plain-text note with the exact settings of that run (readable, plus a
machine-readable block at the end). Keep it in your notebook or as a .txt file: the page's **Load
settings** (or **Load a previous lab note**) can read it back, and so can Claude.

**Not finished in one sitting?** Press **Save session** (top right of the web page, or Cmd+S).
It saves one small `.json` file (settings, the script including hand edits, what was read from the
`.lif`, the half-written lab note, and the run history) — online, confirm the download; from the
double-clicked page it goes to Downloads. Keep it in the project folder. Next time, open the page and
press **Open session** to continue exactly where you stopped. The note next to the buttons says
whether there are changes you have not saved to a file yet.

---

## 5. Protocol B — a different antibody panel / channels / colours

Do Protocol A, but in step 4 also edit sections 2–4 and 8 of the SETTINGS block. Work through
this in order.

### B1 — Find out what is in the file

Channel order and number are fixed at the microscope. Check them in LAS X (or open the `.lif`
in Fiji with Bio-Formats and look at the channel slider). Write down C1, C2, C3… with antibody,
fluorophore and colour you want.

### B2 — `CHANNELS` (section 2): one line per channel, in order

```python
CHANNELS = [
    dict(name="DAPI",        color="blue",    brightness=1.8, background=("pct", 1)),
    dict(name="F4/80-AF488", color="green",   brightness=1.0, background=("pct", 1)),
    dict(name="CD11c-AF647", color="magenta", brightness=1.0, background=("pct", 1)),
]
```

- **Number of lines = number of channels in the file.** Samples with a different channel count
  are skipped with a message.
- `name`: any label (it appears in figures and file columns). Every other setting refers to
  channels by this exact name.
- `color`: `"blue"`, `"red"`, `"green"`, `"magenta"`, `"cyan"`, `"yellow"` or `"gray"`.
- `brightness`: start with `1.0` for everything; adjust after looking at the first results
  (higher number = dimmer).
- `background`: start with `("pct", 1)` for everything (almost no background removal); refine
  later (B6).

### B3 — `FOCUS_CHANNELS` and `LOW_REGION` (section 3)

- `FOCUS_CHANNELS`: 1–3 **bright, structured** channels (nuclei, membranes). They are used to
  pick the sharpest plane and to align tiles. Avoid dim or diffuse channels.
  Example: `FOCUS_CHANNELS = ["DAPI", "F4/80-AF488"]`.
- `LOW_REGION`: only needed if you use `("region", …)` backgrounds (B6). List the channels that
  mark the cells of interest; the "marker-poor" tissue is where these are weakest. Example for a
  T-cell panel: CD3 + CD4; for a macrophage panel: F4/80 + CD11c.

### B4 — `PAIRS` (section 4): the overlay panels

Each line = one panel in the top row of every figure (after the full merge). Use 0–4 pairs.

```python
PAIRS = [
    dict(channels=["F4/80-AF488", "CD11c-AF647"], title="F4/80 / CD11c"),
]
```

Optional `colors={"CD11c-AF647": "cyan"}` recolours a channel in that panel only (useful when two
channels have similar colours, e.g. red + magenta). `PAIRS = []` = merge only.

On the web page the title fills itself in from the channel names in section 2 (e.g.
`CD3-FITC / CD8a-AF594`, with the colours added when a panel is recoloured) and follows along if you
rename a channel. Type over it to use your own title; **Use channel names** switches back.

### B5 — `STANDARD`, `CONTROL_TAG`, `EXCLUDE` (section 5) — as in Protocol A

### B6 — Background rules (the black point), after the first full run

**Recipe for a new panel** (web page: **Neutral first run** / **Suggested starting values** buttons
above the channel rows; each needs a second click and keeps names, colours and controls):

1. **Neutral first run:** brightness `1.0` and `background=("pct", 1)` for every channel. Run once
   and look at `comparison_<group>.png` — the data as recorded.
2. **Pick one rule per channel:** nuclear stain → `("pct", 1)`; channel with a *clean* negative
   control → `("control", 95)` (about control mean + 2 SD; 99 = stricter); no control, marker
   restricted to areas → `("region", 50)`, then adjust 20–75; no control, marker everywhere →
   `("pct", 25)`.
3. **Check:** the negative control / known-negative tissue looks black, and the dimmest cells you
   believe are truly positive are still visible in the STANDARD. Raise the number if hazy, lower it
   if staining disappears; brightness up if washed out, down if faint. One channel at a time.
4. **Lock it** for the whole study and state it in the methods. These are display choices, not
   positivity thresholds — measure on the unfiltered `_raw` TIFFs.

(ZL01's values — CD3 20, CD4/CD8a 75 — were tuned by eye this way; its FITCsec control is not clean
because of eGFP, so CD3 uses the marker-poor rule, not the control rule.)

#### The rules in detail

Look at `comparison_<group>.png` and `display_style_comparison_<group>.png`, then choose per
channel:

| Rule | Effect | Use when |
|---|---|---|
| `("pct", 1)` | 1st percentile of the STANDARD's tissue — almost no removal. | Nuclear stains; bright, clean channels; first run. |
| `("pct", 75)` | 75th percentile of the STANDARD's tissue — strong removal. | Very sparse, dim markers on a hazy background. |
| `("region", 75)` | 75th percentile of this channel in marker-poor tissue (`LOW_REGION`) of the non-control samples. Lower number = milder (ZL01 CD3 used 20). | Markers restricted to certain areas (e.g. T-cell zone). |
| `("control", 95)` | 95th percentile of this channel in **its own** negative-control samples (`control="…"` on that channel). | You included a clean control for this antibody (no primary, secondary only, isotype, leave-one-out). |
| `12.5` | Fixed value. | Matching a previously published display. |
| a list, e.g. `[("pct", 1), ("region", 20)]` | The highest of the rules is used. | Combine a floor with a background measurement. |

Then set `brightness` per channel so the STANDARD looks right (not saturated, structures
visible). Re-run (Run › Run Module, fn+F5) after each change; it takes ~2 minutes.

### B6b — Negative controls (optional, per channel)

If you included a negative control for an antibody — no primary, secondary only, isotype, or a
section stained with everything except that antibody — add a word that appears only in those
sample names to that channel:

```python
dict(name="CD8a-AF594", color="red", brightness=1.5, background=("control", 95),
     control="noCD8"),
```

That control is used **only for that channel** (a CD3 control can never set CD8a's background).
Channels without `control=` simply cannot use the `("control", ..)` rule; for antibodies you trust
without a control, use `("pct", ..)` / `("region", ..)`. All control samples are left out of the
`region` background measurement and listed last in the comparison sheets. `CONTROL_TAG` remains for
controls that are not tied to one channel (optional). On the web page: the **Negative control** box
in each channel row.

### B7 — Extras (section 8)

- `CHECK_CHANNEL`: a bright channel with structure everywhere, shown in the stitching review sheets.
- `COLOC_CHECK`: an optional crop testing whether a `marker` sits on `positive` cells and avoids
  `negative` cells (ZL01: does CD8a sit on CD3+ and avoid CD4+?). Set `COLOC_CHECK = None` if not
  meaningful for your panel.

### B8 — Keep the settings for the whole study

Once a panel's settings are final, **use the same SETTINGS block (except `SRC`, `STANDARD`,
`EXCLUDE`) for every file of that study**, and save a copy of it with your notebook.

---

## 6. Processing several `.lif` files (a batch)

1. **One folder per `.lif`**, each with its own copy of `process_lif.py` + `fiji_macro_template.ijm`.
2. Use the **same SETTINGS** in every copy (only `SRC`, `STANDARD`, `EXCLUDE` differ).
3. Run each (IDLE: File › Open… the next script, Run › Run Module), or from Terminal, per folder:

```bash
cd "/Users/zoeliu/Documents/Claude code/ZL02"
```

```bash
~/.venvs/imaging/bin/python process_lif.py
```

4. **Comparing across files:** display ranges are recalculated from each file's own STANDARD —
   same method, different numbers — so brightness is comparable **within** a file, not
   **between** files. To compare across files: identical microscope settings **and** fixed
   display numbers (use fixed-number backgrounds, e.g. `background=12.5`, and the same
   brightness; or ask for a "reuse ranges from another file" option to be added).

### Running again: overwrite or keep the previous results

`IF_OUTPUT_EXISTS` (web page: Settings › 1) decides what happens when the output folder already
holds results of a previous run:

| Value | What happens |
|---|---|
| `"ask"` (default) | The IDLE shell asks: type **o** = overwrite, **n** = new folder, **c** = cancel, then Enter. |
| `"overwrite"` | The previous results are **deleted** first, then this run is written — no old files mixed in. Only a folder that looks like this script's output is ever deleted; otherwise the run stops with a message. |
| `"new"` | Old results stay untouched; this run writes to `processed_2` (then `_3`, …). The shell prints the name. |

Without this, files with the same name would be replaced but files whose names changed (another
`Z_MODE`, an excluded sample, another STANDARD) would stay behind and mix with the new run.
Overwriting also deletes Fiji results in `processed/fiji/`. Record which settings made which folder.

### Merged z-stack (projection) instead of one plane

By default each image is reduced to its **sharpest single plane**. To combine all planes
instead, set `Z_MODE = "average"` or `Z_MODE = "max"` (web page: Settings › 6) and run as usual.

| | Average | Max |
|---|---|---|
| Each pixel becomes | the mean of its planes | the brightest of its planes |
| Grain | less | more (noise peaks win) |
| Cells from other depths | faint haze | all shown; crowded areas look busier |
| Good for | cleaner figures, comparisons | where positive cells are across the section |
| Files | `*_zavg_raw.tif`, `*_zavg_denoised.tif` | `*_zmax_raw.tif`, `*_zmax_denoised.tif` |

- **Display ranges are recalculated** from the projected STANDARD with the same rules, so samples
  stay comparable within that mode. Do not compare a projection with a single plane.
- **Figure titles** say "average of 3 planes" / "max of …"; `image_summary.csv` › `z_mode` records it.
- **Keep both results:** every run replaces the output folder, so give the projection run its own,
  e.g. `OUT = "processed_average"`.
- **Clear `Z_OVERRIDE`:** a sample forced to one plane is not projected and becomes the odd one
  out; the script prints a `NOTE` at the start if this happens.
- **Sum projection is not offered:** samples with more planes would look brighter regardless of staining.
- **Thick sections:** if samples differ in the number or spacing of planes (ZL01: 3.4–17 µm apart,
  covering 7–22 µm), a projection mixes a different depth of tissue per sample; single planes
  compare more fairly.
- **Fiji:** the generated macro always includes ranges for single plane, Average and Max — choose
  the same mode in its dialog.

---

## 7. Quality control checklist

**Stitching (tile scans)** — `processed/stitching_review/` has one sheet per tile scan: Leica's
merge, the fixed overlap and the registered stitch (plus your own numbers, if set) side by side —
the whole scan and zooms on two seam junctions, all with the same display. Each column shows its
seam score and the word to use; the one in your figures is green (USED).
- [ ] Look through every sheet. Scans that should use a different option go into
      `STITCH_OVERRIDE` (web page: Settings › 7 › *Stitching for individual tile scans*), e.g.
      `STITCH_OVERRIDE = {"scan A": "fixed", "scan B": dict(overlap_LR=13.1, overlap_TB=13.8, drift_LR=9, drift_TB=-10)}`.
      For your own numbers, start from the *Registered* numbers printed on that scan's sheet
      (overlap in %, drift in pixels) and nudge them by 1–2 % or 5–10 px. Doubled cells side by
      side → change `overlap_LR` / `drift_TB`; doubled one above the other → `overlap_TB` / `drift_LR`.
      Re-run and check again. Choosing per scan only moves tiles; pixel values are unchanged, and
      the run record lists every choice.
- [ ] `stitching_seam_scores.csv`: seam scores of every option per scan (higher = better-matching
      seams). `STITCH` method `"best"` picks the highest one per scan automatically.
      (ZL01 10x: registered 0.64–0.86 > Leica. FE71 63x brain: Leica 0.66–0.87 > registered.)
- [ ] Only with method `"registered"` — `stitching.csv`: overlaps ~8–20 %, not at the edge of
      `overlap_search`; `pair_ncc_min` above ~0.5 (ZL01: 0.91–0.93). Low values = seams with
      little structure, where the measured overlap can be wrong.
- [ ] No doubled structures or steps at tile borders in `comparison_<group>.png`.

**Z** — `processed/image_summary.csv`:
- [ ] `z_used` (0 = first plane) looks sensible; blurry images or black holes → force another
      plane with `Z_OVERRIDE`.

**Display** — `comparison_<group>.png`:
- [ ] The STANDARD (top row) looks good: structures visible, background dark, not saturated.
- [ ] Controls look appropriately empty in the channel they control for.
- [ ] `image_summary.csv` → `*_pct_saturated` low (ideally < 1 %).

**Figures** — a few `<group>/*_panels.png`:
- [ ] Scale bars present; title shows z-plane/projection and denoising; bottom line shows ranges.

**Marker check (if used)** — `coloc_check_<group>_<marker>.png`: marker on positive cells, not on
negative cells. (Auto-contrasted — inspection only, not for comparing samples.)

---

## 8. Fiji: the generated macro and doing it by hand

### 8.1 The generated macro (recommended)

Every run writes **`processed/fiji_macro.ijm`** with this file's channel names, colours, display
ranges (single plane, average and max projection) and overlay panels already filled in —
nothing to copy by hand.

1. Fiji: **Plugins › Macros › Run…** → `processed/fiji_macro.ijm`.
2. Dialog:
   - **Mode:** `Batch: processed folder` (all samples) or `Active image` (the open image).
   - **Z handling:** `Sharpest plane (auto)` (same as Python), `Single plane (I choose)`
     (**Fiji counts slices from 1**), `Average projection`, or `Max projection`.
   - **Group:** only for Active mode.
3. Batch mode asks for the `processed` folder; output goes to `processed/fiji/`
   (`*_fiji_processed.tif`, `*_fiji_merge.png`, `*_fiji_panels.png`; names include
   `_zAuto`, `_z2`, `_zAVG`, `_zMAX`). Focus scores appear in Fiji's **Log** window.
4. If you change `DENOISE` in the script, also change `medianRadius` / `gaussSigma` near the top
   of `fiji_macro_template.ijm` (Fiji radius 1 = 3×3).

> The macro was written without access to Fiji and has not yet been test-run. If it stops,
> note the error message and the line number Fiji shows.

### 8.2 By hand in Fiji (one image)

| Step | Fiji menu | Setting |
|---|---|---|
| Open | File › Import › Bio-Formats | Hyperstack, Composite, **Autoscale off** |
| Tile scan stitching | Plugins › Stitching › Grid/Collection stitching | Tile overlap **10**, **Compute overlap** ticked, **Linear Blending** |
| One z-plane | Image › Duplicate | Duplicate hyperstack, one slice |
| — or projection | Image › Stacks › Z Project | Average or Max Intensity (**not** Sum if samples have different plane counts) |
| 32-bit | Image › Type › 32-bit | — |
| Denoise 1 | Process › Filters › Median | Radius **1** |
| Denoise 2 | Process › Filters › Gaussian Blur | Sigma **0.8**, "Scaled units" **unticked** |
| Colours | Image › Lookup Tables | As in `CHANNELS` |
| Display | Image › Adjust › Brightness/Contrast › **Set** | Min/max from `processed/display_ranges.json`. **Never** Auto or Apply. |
| Scale bar | Analyze › Tools › Scale Bar | As in the Python figures |
| Export | Image › Color › Stack to RGB → File › Save As › PNG | — |

---

## 9. Output files reference

All outputs go to `processed/` (see [Running again](#running-again-overwrite-or-keep-the-previous-results)
for what happens to a previous run's results). `<group>` = objective, e.g. `10x`; `<name>` = sample name.

### Per sample — `processed/<group>/`

| File | Contents | Use for |
|---|---|---|
| `<name>_panels.png` | Top row: Merge + `PAIRS` overlays. Bottom row: each channel. Scale bars, z info in title, display ranges at bottom. | Viewing, presentations |
| `<name>_merge.png` | Full-resolution merge (no scale bar). | Figures (add a scale bar) |
| `<name>_zsingle_raw.tif` | Chosen plane, all channels, **unfiltered**, calibrated (`_zavg_` / `_zmax_` in projection mode). | **Intensity measurements** |
| `<name>_zsingle_denoised.tif` | Same after median + Gaussian (`_zavg_` / `_zmax_` in projection mode). | Viewing, segmentation |
| `<name>_zstack_raw.tif` | **All** z-planes (stitched for tile scans), original bit depth, calibrated. | Fiji macro, re-analysis |

### Per file — `processed/`

| File | Contents |
|---|---|
| `comparison_<group>.png` | All samples side by side; STANDARD first, controls last. |
| `display_ranges.json` | Display min/max per channel and group. |
| `display_ranges_projections.json` | Same for single plane / average / max (used by the Fiji macro). |
| `fiji_macro.ijm` | Ready-to-run Fiji macro with this file's settings. |
| `run_settings.json` | **Run record**: the exact settings this run used, run date/time, samples (z-plane used, control), display ranges, stitching, software versions. Attach it to the web page's lab note; the page's **Load settings** can restore it. |
| `image_summary.csv` | Per sample: group, z used / sharpest / count, pixel size, size, control flag, tissue %, per-channel median / 99th percentile / saturation. |
| `stitching.csv` | Tile scans: measured overlaps, drift, worst tile-pair match. |
| `stitching_seam_scores.csv` | Seam scores per tile scan: Leica merge, fixed overlap, registered; which was used; the registered numbers. |
| `stitching_review/<scan>.png` | One sheet per tile scan: every stitching option side by side (whole scan + seam zooms). |
| `coloc_check_<group>_<marker>.png` | Marker check crop (if `COLOC_CHECK` is set). |
| `display_style_comparison_<group>.png` | STANDARD (+ one control): background rules off vs current. |

---

## 10. Settings reference

All in the SETTINGS block of `process_lif.py`. **Anything other than `SRC`, `STANDARD`,
`EXCLUDE`, `Z_OVERRIDE` changes the processing** — record it and keep it constant within a study.

| Section | Setting | ZL01 value | Meaning |
|---|---|---|---|
| 1 | `SRC` | `"ZL01_20260925.lif"` | Input file (same folder as the script). |
| 2 | `CHANNELS` | 4 channels (see file) | Per channel: `name`, `color`, `brightness` (display max = 99.8th pct of STANDARD × brightness; higher = dimmer), `background` (black-point rule(s), [B6](#b6--background-rules-the-black-point-after-the-first-full-run)). |
| 3 | `FOCUS_CHANNELS` | DAPI, CD3, CD4 | Channels for focus scoring and tile alignment. |
| 3 | `LOW_REGION` | CD3 + CD4, 0.15 | Marker-poor tissue = lowest 15 % of these (smoothed) within tissue. |
| 4 | `PAIRS` | 3 overlays | Overlay panels; optional per-panel `colors`. |
| 5 | `STANDARD` | 1:50 at 10x and 63x | Reference sample per group. |
| 2 | `control=` (in `CHANNELS`) | CD3-FITC: `"FITCsec"` | Optional per-channel negative-control word; required for that channel's `("control", ..)` rule. |
| 5 | `CONTROL_TAG` | `""` | Optional other controls not tied to one channel (left out of `region` backgrounds, listed last). |
| 5 | `EXCLUDE` | `["LN01", "Series001"]` | Skipped samples. |
| 6 | `Z_MODE` | `"sharpest"` | `"sharpest"` (files `*_zsingle_*`), `"average"` (`*_zavg_*`) or `"max"` (`*_zmax_*`); see [Merged z-stack](#merged-z-stack-projection-instead-of-one-plane). |
| 6 | `Z_OVERRIDE` | `{}` | Force a single plane per sample, **0 = first plane** (Fiji slice 1). Also applies in projection mode. |
| 7 | `DENOISE` | median 3, gaussian 0.8 | Filters (median 1 / gaussian 0 = off). |
| 7 | `HI_PCT` | 99.8 | Percentile for the display max. |
| 7 | `STITCH` | method `"leica"`; start 0.10, search 0.05–0.22, drift 30 px | `"leica"` = Leica's merged image; `"fixed"` = raw tiles at the start overlap; `"registered"` = overlap measured here with the search settings; `"best"` = highest seam score per scan. |
| 7 | `STITCH_OVERRIDE` | `{}` | Per tile scan: one of the words above, or `dict(overlap_LR=…, overlap_TB=…, drift_LR=…, drift_TB=…)` (% and pixels; `overlap=` sets both). |
| 7 | `SCALEBAR_UM` | `{}` (auto: 200 µm at 10x, 20 µm at 63x) | Fixed scale bars per group, e.g. `{"10x": 200}`. |
| 7 | `SCALEBAR_POS` | `"lower right"` | Scale bar position in all figures (and the Fiji macro): `"lower right"`, `"lower left"`, `"upper right"`, `"upper left"` or `"center"`. |
| 8 | `MAKE_EXTRAS`, `CHECK_CHANNEL`, `COLOC_CHECK` | on, CD3, CD8a vs CD3/CD4 at 63x | Check figures. |
| 8 | `OUT` | `"processed"` | Output folder name; use a different one (e.g. `"processed_average"`) to keep results of another `Z_MODE`. |
| 8 | `IF_OUTPUT_EXISTS` | `"ask"` | `"ask"`, `"overwrite"` (delete previous results first) or `"new"` (write to `OUT_2`, …); see [Running again](#running-again-overwrite-or-keep-the-previous-results). |

Fixed choices in the code: focus score = variance of the Laplacian of the Gaussian-blurred
(σ 1.5 px) sum of focus channels; tissue mask = smoothed sum of channels above 0.35 × Otsu;
marker-poor smoothing σ 20 px for coarse pixels (> 0.5 µm) or 25 px for fine pixels; seam
blending linear; groups from the objective magnification in the file (or pixel size if missing).

Previous (ZL01-only) versions of the script and macro are kept in `archive/`.

---

## 11. Notes on the ZL01 example

- **Foxp3-IRES-eGFP mice:** eGFP is detected in the FITC channel, so CD3-FITC always contains
  some Treg eGFP. The `FITCsec` samples (CD3 secondary-only) contain eGFP and bright non-T cells
  and are **not** a clean FITC background — they are assigned as CD3-FITC's negative control
  (`control="FITCsec"`, which also keeps them out of the `region` background) but no `control` rule
  uses them (it would give a CD3 black point of ~61 instead of ~10). There are no negative
  controls for CD8a-AF594 or CD4-AF647 in ZL01.
- **CD8a-AF594 (53-6.7, frozen):** specific but dim (on CD3+CD4− membranes, not on CD4+ cells;
  0–3 photons/pixel at ≈ 3 % 561 nm, no averaging). 1:50 ≈ 1:100. Fix at acquisition.
- **Tile-grid pattern** in dim channels at 10x = vignetting, not stitching; needs flat-field correction.
- **Thick sections:** planes 3.4–17 µm apart covering 7–22 µm, differing per sample.
- **DAPI differs between slides** at identical settings (staining/section, not processing).
- **Excluded:** LN01 (different laser/gain) and Series001 (unnamed 63x).

---

## 12. Fix the image — what to change when the figure does not look right

Find what you see, change that one setting, re-run, compare. **Change one thing at a time**, keep
every change identical for all samples, and write down what you changed. "Settings › n" = section n
of the web page; the script names are in brackets.

### Brightness and background — Settings › 2 (`CHANNELS: brightness, background`)

| You see | Change | Notes |
|---|---|---|
| One channel too bright: positive cells blown out, white blobs | **brightness ↑** for that channel, e.g. 1.5 → 1.8 | Display max = 99.8th percentile of the STANDARD × brightness, so a **higher number = dimmer**. Change in steps of 0.2–0.3. |
| One channel too dim: positive cells hard to see | **brightness ↓**, e.g. 1.6 → 1.3 | Below ~1 the brightest cells saturate (flat white). If it is still faint at ~1, the signal is weak at acquisition (see the last table). |
| Haze or glow over the whole tissue; background looks grey or grainy | **Raise the black point**: a higher percentile (`("pct", 1)` → `("pct", 25)`, or `("region", 20)` → `("region", 50)`), or add a `("region", …)` rule | Higher percentile = more background set to black. With two rules the **higher** value wins. See **display_style_comparison_…png** for the effect. |
| Dim real cells vanish; image looks too clean, with black holes in tissue | **Lower the black point**, e.g. `("pct", 75)` → `("pct", 50)`, or remove the rule that gives the high value | Check the dimmest real positives in the STANDARD: they must stay visible. |
| Negative-control sample not black | Channel with a control: `("control", 95)` → `("control", 99)`. No control rule yet: add `control="word"` and `("control", 95)` | Only if the control is a clean negative (e.g. not an eGFP mouse in the FITC channel). |
| DAPI too strong, nuclei are flat blobs | DAPI **brightness ↑**, e.g. 1.8 → 2.2 | Keep DAPI's black point low (`("pct", 1)`). DAPI only needs to show where the tissue is. |
| Two channels look alike in the merge (red + magenta, blue + cyan) | Change a channel's **colour** (Settings › 2), or recolour it only in an overlay panel (Settings › 4) | Colour does not change any numbers, only the display. |

### All channels at once — Settings › 7 (`HI_PCT`)

| You see | Change | Notes |
|---|---|---|
| Every channel slightly too bright, lots of saturated pixels | **display max percentile ↑**, 99.8 → 99.9 | Prefer per-channel brightness. This changes all channels together. |
| Every channel slightly too dim | **display max percentile ↓**, 99.8 → 99.5 |  |

### One sample looks different — Settings › 5 (`STANDARD, EXCLUDE`)

| You see | Change | Notes |
|---|---|---|
| One sample is much brighter or dimmer than the others | Usually **nothing**: identical display settings show real differences | Never adjust one sample separately, or comparisons are lost. |
| Every sample of a group looks dim (or all washed out) | Choose a more typical **STANDARD** (Settings › 5) | Every range comes from the STANDARD. A STANDARD with a bright artifact (dust, fold, edge) makes everything dim. An unusually faint one makes everything too bright. |
| A damaged or failed sample distorts the ranges or panels | Add a word from its name to **EXCLUDE** (Settings › 5) | Excluded samples are not processed at all. |

### Noise and grain — Settings › 7 (`DENOISE, Z_MODE`)

| You see | Change | Notes |
|---|---|---|
| Salt-and-pepper speckles (single bright pixels) | **median size ↑**, 3 → 5 | Odd numbers only; 1 = off. Larger values erase thin processes and small dots. |
| Fine grain over everything | **Gaussian σ ↑**, 0.8 → 1.2 | Softens the image. Above ~1.5 it looks blurry. |
| Image too soft; thin processes or small puncta lost | **Gaussian σ ↓** (0.8 → 0.5) and/or **median size** 3 → 1 | Accept more grain in exchange for detail. |
| One dim channel stays grainy whatever you do | **Z handling → Average projection** (Settings › 6) | Averaging planes reduces grain, but adds some haze. The real fix is at the microscope: line averaging or accumulation, or more laser or gain. |
| Grainy specks in the background | **Raise that channel's black point** (Settings › 2) | Background noise below the black point becomes black. |

### Focus and z-planes — Settings › 6 (`Z_MODE, Z_OVERRIDE, FOCUS_CHANNELS`)

| You see | Change | Notes |
|---|---|---|
| One sample is out of focus: the wrong plane was picked | **Z_OVERRIDE**: add the sample and the plane you want (Settings › 6) | Planes count from **0** here. Fiji's slice n = plane n − 1. The focus scores are printed while the script runs. Look through **<name&gt;_zstack_raw.tif** in Fiji to choose a plane. |
| The wrong plane is picked in many samples | Change **FOCUS_CHANNELS** (Settings › 3): use channels with sharp, small structures everywhere (DAPI, membrane markers) | Drop diffuse or sparse channels; they judge focus poorly. |
| Thick section: different cells in focus in different planes | **Z handling → Average** (smooth) or **Max** (brighter, grainier) | Display ranges are recalculated for projections. Change the output folder (Settings › 1) to keep both versions. |
| Max projection too grainy or too bright | Use **Average** instead, or raise the black points | Max picks the brightest noise too. |

### Stitching (tile scans) — Settings › 7 (`STITCH`)

| You see | Change | Notes |
|---|---|---|
| Doubled cells or steps at tile borders | Switch **Tile scans**: Leica's merged image ↔ stitch the raw tiles here | Compare them in **stitching_review/<scan>.png**. The one marked USED (green) is in your figures. Higher seam score = better seams; method `"best"` picks it per scan. |
| Some scans look best with one option, others with another | `STITCH_OVERRIDE`: add the scan and its option (Settings › 7 › Stitching for individual tile scans) | Review every sheet first, then set only the scans that differ from the default. |
| A scan is almost right with every option | Your own numbers for that scan, starting from the Registered numbers on its sheet | Doubled side by side → `overlap_LR` / `drift_TB`; doubled one above the other → `overlap_TB` / `drift_LR`. Steps of 1–2 % or 5–10 px. |
| Raw-tile stitching: the overlap in stitching.csv sits at the edge of the search (e.g. exactly 5 % or 22 %) | **Widen** search from/to (e.g. 0.03–0.30) | The true overlap was outside the range that was searched. |
| Raw-tile stitching: overlaps differ a lot between scans, or pair_ncc_min < 0.5 | **Narrow** search from/to around the overlap set in LAS X (10 % → 0.07–0.15), or use Leica's merge | Seams with little structure (sparse tissue, 63x) fool the matching. A narrow range limits the damage. |
| Raw-tile stitching: tiles shifted sideways | **sideways (px) ↑**, 30 → 50 | Allows a larger sideways offset between neighbouring tiles. |
| Regular grid of brighter and darker squares, but no doubling | Not a stitching problem: uneven illumination (vignetting) | Needs flat-field correction (not built in). A slightly higher black point hides it in dim channels. |

### Marker-poor tissue rule — Settings › 3 (`LOW_REGION (only for ("region", …) black points)`)

| You see | Change | Notes |
|---|---|---|
| The region rule removes too much: real positives get dimmer | **Lower that region percentile**, or **fraction ↓** (0.15 → 0.10) | The "marker-poor" area is not truly negative. Also check that the LOW_REGION channels really are absent from part of the tissue. |
| The region rule barely changes anything | **Region percentile ↑**, or **fraction ↑** (0.15 → 0.25) |  |

### Figure layout — Settings › 7 (`SCALEBAR_UM, SCALEBAR_POS, PAIRS`)

| You see | Change | Notes |
|---|---|---|
| Scale bar covers something | **Scale bar position** (Settings › 7): another corner or the centre |  |
| Scale bar length awkward | Set a length per group (Settings › 7), e.g. 63x → 20 µm | Empty = chosen automatically. |
| Want other channel combinations on top | **Overlay panels** (Settings › 4): add, remove or recolour | Titles follow the channel names unless you type your own. |

### Cannot be fixed by processing (fix at the microscope next time)

| You see | Change | Notes |
|---|---|---|
| Flat white areas in the raw data (pixels at 255) | Lower the gain or laser power | Saturated pixels hold no information. Processing cannot bring detail back. |
| Photon-starved channel: grainy even with good black points | Line averaging or accumulation, more laser, slower scan | Denoising only hides grain. |
| One channel shows another channel's pattern (bleed-through) | Sequential scanning; single-stain controls |  |
| Signal that does not belong to the antibody (autofluorescence, a reporter like eGFP) | Pick another fluorophore or channel, or include a matching control | E.g. Foxp3-eGFP sits in the FITC channel. |

---

## 13. Troubleshooting

| Message / problem | Cause | Fix |
|---|---|---|
| `No module named readlif` (numpy, skimage…) / `WRONG PYTHON` | IDLE was opened from Applications/Dock (plain Python). | Close it; double-click `Open IDLE for imaging.command` in the `image` folder, or `~/.venvs/imaging/bin/python -m idlelib` in Terminal. |
| `SETTINGS PROBLEM: cannot find XXX.lif` | `SRC` misspelled / file elsewhere. | The message lists the `.lif` files it sees. |
| `SETTINGS PROBLEM: STANDARD for group …` | Name doesn't match. | Copy the suggested line, pick your standard from the list. |
| `SETTINGS PROBLEM: 'X' in PAIRS is not a channel name` | A name in `FOCUS_CHANNELS`/`LOW_REGION`/`PAIRS`/extras differs from `CHANNELS`. | Use exactly the `name` from `CHANNELS` (the message lists them). |
| `SETTINGS PROBLEM: colour …` | Unsupported colour. | Use blue, red, green, magenta, cyan, yellow, gray. |
| `skipped …: N channels, but CHANNELS lists M` | That image has a different channel count. | Fine if intended; otherwise fix `CHANNELS`. |
| `SETTINGS PROBLEM: no images left` | Everything excluded or channel count wrong. | Check `EXCLUDE` and `CHANNELS`. |
| `… uses ('control', ..) but has no negative control` | `control` rule on a channel without `control=`. | Add `control="word"` to that channel, or use another rule. |
| `…: no processed sample … contains its control word` | The channel's control word matches no sample. | Check the spelling against the sample names. |
| `extra skipped: …` | A check figure can't be made (e.g. no tile scans). | Main outputs are fine. |
| `Fiji macro NOT written: …` | `fiji_macro_template.ijm` missing from the folder. | Copy it next to the script. |
| `SyntaxError` after editing | Quote, comma or bracket deleted. | IDLE highlights the line; compare with the examples here. |
| Seams / doubled structures | Stitching failed. | Compare the options in `stitching_review/<scan>.png`; change `STITCH` method or set single scans in `STITCH_OVERRIDE`; with registered, check `stitching.csv` / widen the search; or stitch in Fiji (8.2). |
| Figures too dark / bright / hazy | Display settings. | Adjust `brightness` / `background` (B6); check the STANDARD is well stained. |
| Very slow / freezes | Large file, all in memory. | Close other programs; one file at a time; `MAKE_EXTRAS = False`. |
| Stop a run | — | IDLE: Shell › Interrupt Execution (Ctrl + C); Terminal: Ctrl + C. |

---

## 14. Methods text (template)

Adapt bracketed parts.

> Confocal images were acquired on a Leica [SP-series] microscope with HyD detectors ([8]-bit,
> sequential acquisition) using [10×/0.30 dry] (tile scans) and [63×/1.40 oil] objectives. Images
> were processed with a custom Python pipeline (readlif, scikit-image, tifffile). Tile scans were
> merged in LAS X [or, with method "registered": stitched from the raw tiles by normalized
> cross-correlation of tile overlaps (initial overlap 10%) with linear blending]. For each image, the sharpest optical section was selected (maximal
> variance of the Laplacian of the combined [DAPI, CD3 and CD4] signal). Each channel was denoised
> with a 3 × 3 median filter followed by a Gaussian filter (σ = 0.8 pixels). Display ranges were set
> identically for all samples of a given magnification, derived from a reference sample; black
> points were set at [the 75th/20th/75th percentile of CD8a/CD3/CD4 signal in T-cell-poor regions].
> Intensity measurements were performed on unfiltered images.

---

## 15. Glossary

| Term | Meaning |
|---|---|
| **`.lif`** | Leica's file format; one file holds many images ("series"). |
| **Group** | All images taken with the same objective (e.g. `10x`); each group has its own display settings. |
| **Tile scan / stitching** | Large area imaged as overlapping tiles, then joined; "compute overlap" = measuring the real overlap. |
| **z-plane / z-stack / projection** | One optical section / sections at several depths / combining them (average or brightest per pixel). |
| **Median filter** | Each pixel → middle value of its 3×3 neighbourhood; removes isolated noisy pixels. |
| **Gaussian blur (σ)** | Weighted average of neighbours; σ = reach in pixels. |
| **Display range (min/max)** | Values shown as black (min, "black point") and full colour (max). Changes the look, not the data. |
| **STANDARD** | Reference sample whose pixel values set the display ranges for its group. |
| **Percentile** | Value below which that % of pixels fall (median = 50th). |
| **NCC** | Normalized cross-correlation, a similarity score up to 1 (1 = identical), used for tile alignment. |
| **Flat-field correction** | Dividing by the microscope's illumination pattern to remove uneven brightness (vignetting). |
| **Saturation** | Pixels at the maximum value (255 in 8-bit, 4095 in 12-bit); true brightness lost. |
