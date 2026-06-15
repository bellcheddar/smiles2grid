# ⚗️ smiles_grid

> **Turn a JSON of SMILES into a paginated, property-annotated PDF grid with RDKit.**

![python](https://img.shields.io/badge/python-3.9+-3776AB?logo=python&logoColor=white) ![RDKit](https://img.shields.io/badge/RDKit-depiction-00897B) ![output](https://img.shields.io/badge/output-PDF%20+%20CSV-467FF7) ![input](https://img.shields.io/badge/input-JSON%20SMILES-9b51e0) ![author](https://img.shields.io/badge/author-Marc%20C.%20Deller%2C%20D.Phil.-1C244B)

<table>
<tr>
<td>🌐 <b>Website</b></td><td><a href="https://marcdeller.com" target="_blank" rel="noopener noreferrer">marcdeller.com</a></td>
<td>✉️ <b>Contact</b></td><td><a href="mailto:marc@marcdeller.com">marc@marcdeller.com</a></td>
<td>🐙 <b>GitHub</b></td><td><a href="https://github.com/bellcheddar/smiles2grid" target="_blank" rel="noopener noreferrer">bellcheddar/smiles2grid</a></td>
</tr>
</table>

---

Render SMILES strings from a JSON file into a paginated PDF grid using RDKit.

Why it matters: after a generative or screening run you often have hundreds of SMILES and no fast way to actually look at them together, compare scaffolds, or pull a labelled subset for a slide. smiles_grid renders them into clean 5x5 PDF pages with structure labels tied to the original input order, computes a full property panel (MW, cLogP, HBD/HBA, TPSA, and more) into a companion CSV, and can highlight near-neighbours of a query structure in red. It is useful for medicinal chemists and computational scientists who need a quick, shareable visual of a compound set with the numbers attached: ideal for triage, SAR discussions, and decision-making meetings.

## What it does

- Extracts SMILES strings from a nested JSON file.
- Preserves the original input order for labeling and selection.
- Renders molecules into a boxed 5x5 grid, 25 compounds per page.
- Supports optional similarity search against a query SMILES.
- Supports range-based selection such as `all`, `1-25`, `50-75`, or `1,5,10-20`.
- Writes a CSV summary for every rendered record.
- Highlights closely related structures within a page in red.

## Example input

```bash
python smiles_grid.py scratch_Baricitinib_sim.out.20260530_0219.json
```

The script prompts by default. You can also run it non-interactively with `--no-prompt` and selection arguments.

## Features

### PDF output

- Single PDF output.
- 5x5 bordered panels per page.
- Structure labels based on original JSON order.
- Clean footer text with predicted properties.
- Structure drawing scaled to fit each panel.

### CSV output

A CSV file is written alongside the PDF and includes:

- output index.
- original index.
- source index.
- label.
- SMILES and canonical SMILES.
- formula.
- MW, cLogP, HBD, HBA, TPSA.
- rotatable bonds, rings, heavy atoms, fraction sp3, aromatic rings, hetero atoms.
- input JSON path.
- PDF filename.
- selection mode.

## Usage

### Interactive mode

```bash
python smiles_grid.py input.json
```

At runtime the script prompts for:

1. An optional query SMILES for similarity search.
2. A similarity cutoff percentage.
3. If no query is provided, a range selection such as `all`, `1-25`, or `50-75`.

### Similarity search

```bash
python smiles_grid.py input.json --query-smiles "CCO" --similarity 75
```

### Range selection

```bash
python smiles_grid.py input.json --range 1-25
python smiles_grid.py input.json --range all
python smiles_grid.py input.json --range 1,5,10-20
```

### Non-interactive mode

```bash
python smiles_grid.py input.json --no-prompt --range 1-25
```

## Output files

By default the script writes:

- `input_grid.pdf`
- `input_summary.csv`

You can override the PDF filename with `--output`.

## Requirements

- Python 3.9+
- RDKit
- ReportLab
- Pillow
- cairosvg, if your RDKit build does not include Cairo drawing support

## Installation

### Conda

```bash
conda install -c conda-forge rdkit reportlab pillow cairosvg
```

### Pip-style environments

If you already have RDKit available in your environment:

```bash
pip install reportlab pillow cairosvg
```

## Command-line options

- `--query-smiles`: similarity search query.
- `--similarity`: similarity cutoff as `0-100` or `0-1`.
- `--range`: selection by range or comma-separated ranges.
- `--no-prompt`: disable interactive prompting.
- `--max-pages`: render only the first N pages for debugging.
- `--label-prefix`: change the label prefix, default `J`.
- `--output`: choose the PDF filename.

## Notes

- Numbering is tied to the original JSON order, not to whether a SMILES can be rendered.
- Failed renderables are skipped visually but do not shift labels.
- The script keeps the input order for selection and CSV reporting.
- Similarity highlighting is page-local and uses red bond highlighting.

## Example workflow

1. Place the JSON file in your working directory.
2. Run the script.
3. Enter a query SMILES or press Enter.
4. Enter a similarity cutoff or press Enter.
5. Or press Enter again and provide a range such as `all` or `1-25`.
6. Open the generated PDF and CSV.

## License

Add your preferred license here.

---

## 👤 Author

**Marc C. Deller, D.Phil.**  
Structural biologist & drug discovery scientist  

<table>
<tr>
<td>🌐</td><td><a href="https://marcdeller.com" target="_blank" rel="noopener noreferrer">marcdeller.com</a></td>
<td>✉️</td><td><a href="mailto:marc@marcdeller.com">marc@marcdeller.com</a></td>
<td>🐙</td><td><a href="https://github.com/bellcheddar/smiles2grid" target="_blank" rel="noopener noreferrer">github.com/bellcheddar/smiles2grid</a></td>
</tr>
</table>
