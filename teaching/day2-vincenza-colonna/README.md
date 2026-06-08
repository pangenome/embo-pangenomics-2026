# Pangenomics Sequencing Primer

Source material for a 60-minute levelling lecture on genome sequencing technologies and data types for the EMBO Practical Course on Pangenomics.

## Contents

- `handout.qmd`: handout entry point, metadata, and ordered section includes.
- `handout/`: editable long-form handout sections.
- `slides.qmd`: concise Reveal.js presentation with speaker notes.
- `styles.css`: handout-specific formatting.
- `figures/`: images shared by the handout and slides.
- `_quarto.yml`: Quarto project configuration.

The `.qmd` files are plain-text Markdown with Quarto metadata and extensions, so they remain readable and reviewable in Git.

Files under `handout/` follow the document structure: `00` contains the learning objectives, `01` through `11` match the numbered handout sections, and `12` through `13` contain the sources and production notes.

## Rendering

Install [Quarto](https://quarto.org/docs/get-started/) and run:

```bash
quarto render
```

In this workspace the Quarto cache needs to stay inside the repository, so the
local wrapper is the reliable command:

```bash
bash render-local.sh
```

Generated files are written to `_output/`. To render one document:

```bash
quarto render handout.qmd
quarto render slides.qmd
```

The default build creates an HTML handout and Reveal.js slides. In the slide presentation, press `S` to open speaker view.

To generate a PDF handout separately:

```bash
quarto install tinytex
quarto render handout.qmd --to pdf
```

## Editing workflow

1. Edit detailed teaching content in the relevant file under `handout/`; update the include order in `handout.qmd` when adding or moving sections.
2. Adapt the key messages in `slides.qmd`; do not copy full handout paragraphs onto slides.
3. Store image source files and final assets in `figures/` with descriptive names.
4. Run `quarto render` and inspect both outputs before committing.

## Authorship and AI assistance

Vincenza Colonna is the author and subject-matter reviewer. OpenAI Codex assisted with drafting, editing, slide adaptation, and initial Quarto configuration. AI assistance is disclosed in the handout and is not treated as authorship; the named author retains responsibility for checking and using the material.


## 60-minute structure

| Time | Section | Goal |
|---:|---|---|
| 0–5 min | Why sequencing technology matters for pangenomics | Motivate the lecture with examples of variants invisible to a single linear reference. |
| 5–12 min | From molecule to digital read | Define library, read, basecalling, quality score, coverage, and metadata. |
| 12–23 min | Short-read sequencing | Compare Illumina and Ultima short-read sequencing, and explain what short reads are excellent at. |
| 23–37 min | Long-read single-molecule sequencing | Compare PacBio HiFi and Oxford Nanopore sequencing. |
| 37–45 min | Complementary data types | Discuss Hi-C, optical maps, linked reads, RNA, and single-cell data. |
| 45–54 min | Data formats as biological evidence | Walk through the ladder from raw signal to graph. |
| 54–60 min | Choosing data for a pangenomics project | Summarize decision rules and transition to practical sessions. |

---
