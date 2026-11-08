# `figures/` — Thesis Figures, Code & Caption Registers

This directory contains the reproducible figure-generation scripts (`R` and `Python`), numbering registers, and exported figure plates (`PDF`) for **Chapter 3** and **Chapter 4** of the thesis.

```text
figures/
├── README.md
├── chapter3/
│   ├── chapter3_figure_numbering.md      # Sequential register: Fig 3.1–3.15 & Supp Fig 3.1–3.44
│   ├── Chapter3_figure_captions.md       # Authoritative Chapter 3 figure legends
│   ├── code/                             # Python & R figure renderers for Chapter 3
│   └── exports/                          # Exported Chapter 3 PDF figure plates
└── chapter4/
    ├── chapter4_figure_numbering.md      # Sequential register: Fig 4.1–4.8 & Supp Fig 4.1–4.16
    ├── code/                             # Python & R figure renderers for Chapter 4
    └── exports/                          # Exported Chapter 4 PDF figure plates
```

---

## Naming & Numbering Conventions

1. **Sequential First-Mention Thesis Numbering:**
   - **Chapter 3** comprises **15 main figures (`Fig 3.1–3.15`)** and **44 supplementary figures (`Supp Fig 3.1–3.44`)**. See [`chapter3/chapter3_figure_numbering.md`](chapter3/chapter3_figure_numbering.md) for the exact crosswalk between sequential thesis figure numbers, section anchors, descriptive filename stems (`ch03_main_*` / `ch03_supp_*`), and legacy provisional IDs (`CH3-S3.*`).
   - **Chapter 4** comprises **8 main figures (`Fig 4.1–4.8`)** and **16 supplementary figures (`Supp Fig 4.1–4.16`)**. See [`chapter4/chapter4_figure_numbering.md`](chapter4/chapter4_figure_numbering.md) for the full register.
2. **Descriptive Filename Stems:**
   - Main figure exports carry `ch03_main_<stem>` or `ch04_main_<stem>`; supplementary figure exports carry `ch03_supp_<stem>` or `ch04_supp_<stem>` (with methods QC figures prefixed `ch03_methods_*`).
3. **FullUpstream Source Data:**
   - Complete high-resolution `PDF`/`SVG`/`PNG` sets and figure-specific `source_data/` folders reside in `A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/figures/chapter3/` and `chapter4/`.
