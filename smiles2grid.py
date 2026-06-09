#!/usr/bin/env python3
"""
smiles_grid.py

Render SMILES extracted from a JSON file into a single PDF with boxed 5x5 grids
(25 compounds per page) using RDKit and ReportLab.

Usage modes
-----------
Prompt mode:
    python smiles_grid.py scratch_Baricitinib_sim.out.20260530_0219.json

Range mode:
    python smiles_grid.py scratch_Baricitinib_sim.out.20260530_0219.json --range 1-25
    python smiles_grid.py scratch_Baricitinib_sim.out.20260530_0219.json --range all

Similarity mode:
    python smiles_grid.py scratch_Baricitinib_sim.out.20260530_0219.json \
        --query-smiles CCO --similarity 75

Default mode:
    If no selection arguments are provided, the script prompts the user by
    default (you do not need to pass --prompt).

Prompt behavior
---------------
- Interactive prompts are the default behavior.
- Use --no-prompt to suppress prompts and rely on command-line selection.
- First prompt: optional query SMILES for similarity search.
- Second prompt: similarity cutoff percentage.
- If a query SMILES is supplied, range selection is skipped.
- If the query prompt is left blank, the script asks for a range string.
- Range examples: all, 1-25, 50-75, 1,5,10-20.
- Original input-file order is preserved for selection, numbering, and labels.

Command-line behavior
---------------------
- --query-smiles activates similarity selection.
- --similarity accepts either 0-100 or 0-1.
- --range accepts all / ranges / comma-separated ranges.
- --no-prompt disables prompts and uses command-line selection only.
- --max-pages is a debug option to render only the first N pages.
- --label-prefix changes labels, default J -> J0001, J0002, ...
- --output sets the PDF filename.

Customization notes
-------------------
- Change GRID_ROWS / GRID_COLS to alter page layout.
- Change PAGE_SIZE to reportlab.lib.pagesizes.A4 / letter as desired.
- Change DEFAULT_SIMILARITY_THRESHOLD for the default query cutoff.
- Change PAGE_LOCAL_SIMILARITY_THRESHOLD to adjust red similarity highlighting
  within each displayed 25-compound page.
- Change HIGHLIGHT_COLOR and VARIABLE_BOND_WIDTH_MULTIPLIER to alter highlight style.
- Change descriptor_lines() to alter the properties shown at the bottom of each box.
- Set CLEAN_FRAGMENTS=False if you want to keep salts/counterions instead of using
  the largest fragment for descriptors, similarity, and drawing.
- Atom numbering is disabled.
- If RDKit lacks MolDraw2DCairo, the script falls back to MolDraw2DSVG and will
  use cairosvg if available to convert SVG to PNG for ReportLab.

Recommended install:
    conda install -c conda-forge rdkit reportlab pillow cairosvg
"""

from __future__ import annotations

import argparse
import io
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Any

from PIL import Image
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdFMCS, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas

RDLogger.DisableLog("rdApp.*")

try:
    import cairosvg  # type: ignore
except Exception:
    cairosvg = None

PAGE_SIZE = A4
GRID_ROWS = 5
GRID_COLS = 5
ITEMS_PER_PAGE = GRID_ROWS * GRID_COLS
MARGIN = 24
CELL_GAP = 8
LABEL_FONT = "Helvetica-Bold"
TEXT_FONT = "Helvetica"
LABEL_SIZE = 10
TEXT_SIZE = 5.8
BORDER_WIDTH = 0.8
DEFAULT_LABEL_PREFIX = "J"
DEFAULT_SIMILARITY_THRESHOLD = 0.70
PAGE_LOCAL_SIMILARITY_THRESHOLD = 0.68
VARIABLE_BOND_WIDTH_MULTIPLIER = 40
HIGHLIGHT_COLOR = (1.0, 0.0, 0.0)
DRAW_SIZE = (700, 520)
CLEAN_FRAGMENTS = True
TOP_PAD = 16
INNER_PAD_X = 6
INNER_PAD_Y = 4
TITLE_BLOCK_H = 18
FOOTER_BLOCK_H = 34
PADDING = 8
FOOTER_TEXT_GAP = 6
CSV_SUFFIX = "_summary.csv"


@dataclass
class MolEntry:
    original_index: int
    source_index: int
    smiles: str
    mol: Chem.Mol
    fp: DataStructs.cDataStructs.ExplicitBitVect
    record: Dict[str, Any]


def largest_fragment(mol: Chem.Mol) -> Chem.Mol:
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) <= 1:
        return mol

    def score(m: Chem.Mol) -> Tuple[bool, int, int]:
        heavy = m.GetNumHeavyAtoms()
        carbons = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 6)
        return (carbons > 0, heavy, carbons)

    return sorted(frags, key=score, reverse=True)[0]


def prepare_mol_from_smiles(smiles: str) -> Optional[Chem.Mol]:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        if CLEAN_FRAGMENTS:
            mol = largest_fragment(mol)
        Chem.SanitizeMol(mol)
        mol = Chem.RemoveHs(mol)
        AllChem.Compute2DCoords(mol)
        return mol
    except Exception:
        return None


def iter_smiles_candidates(obj) -> Iterable[str]:
    if isinstance(obj, dict):
        if isinstance(obj.get("candidates"), dict):
            for key in obj["candidates"].keys():
                if isinstance(key, str):
                    yield key

        preferred_keys = (
            "smiles",
            "SMILES",
            "canonical_smiles",
            "canonicalSmiles",
            "candidate_smiles",
            "candidate",
            "structure",
            "molecule",
            "mol",
        )
        for key in preferred_keys:
            value = obj.get(key)
            if isinstance(value, str):
                yield value

        for value in obj.values():
            yield from iter_smiles_candidates(value)
    elif isinstance(obj, list):
        for item in obj:
            yield from iter_smiles_candidates(item)


def extract_smiles_from_json(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)

    seen: set[str] = set()
    smiles_list: List[str] = []

    for candidate in iter_smiles_candidates(data):
        if candidate in seen:
            continue
        seen.add(candidate)
        smiles_list.append(candidate)

    return smiles_list


def mol_properties(mol: Chem.Mol) -> Dict[str, Any]:
    return {
        "mw": round(Descriptors.MolWt(mol), 3),
        "clogp": round(Crippen.MolLogP(mol), 3),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 3),
        "rot_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "rings": int(rdMolDescriptors.CalcNumRings(mol)),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "fsp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        "aromatic_rings": int(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "hetero_atoms": int(rdMolDescriptors.CalcNumHeteroatoms(mol)),
    }


def build_entries(smiles_list: Sequence[str]) -> List[MolEntry]:
    entries: List[MolEntry] = []
    for source_idx, smiles in enumerate(smiles_list, start=1):
        mol = prepare_mol_from_smiles(smiles)
        if mol is None:
            continue
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        record = mol_properties(mol)
        entries.append(MolEntry(original_index=len(entries) + 1, source_index=source_idx, smiles=smiles, mol=mol, fp=fp, record=record))
    return entries


def page_selection_chunks(entries: Sequence[MolEntry], page_size: int = ITEMS_PER_PAGE) -> List[List[MolEntry]]:
    return [list(entries[i:i + page_size]) for i in range(0, len(entries), page_size)]


def descriptor_lines(mol: Chem.Mol) -> List[str]:
    props = mol_properties(mol)
    return [
        f"MW {props['mw']:.1f}  cLogP {props['clogp']:.2f}",
        f"HBD {props['hbd']}  HBA {props['hba']}  TPSA {props['tpsa']:.1f}",
        f"RotB {props['rot_bonds']}  Rings {props['rings']}  Het {props['hetero_atoms']}",
    ]


def format_label(prefix: str, n: int) -> str:
    return f"{prefix}{n:04d}"


def parse_range_selection(text: str, max_index: int) -> List[int]:
    text = text.strip().lower()
    if not text or text == "all":
        return list(range(1, max_index + 1))

    selected = set()
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            start = int(a)
            end = int(b)
            if start > end:
                start, end = end, start
            start = max(1, start)
            end = min(max_index, end)
            for i in range(start, end + 1):
                selected.add(i)
        else:
            idx = int(part)
            if 1 <= idx <= max_index:
                selected.add(idx)
    return sorted(selected)


def prompt_similarity_mode() -> Tuple[Optional[Chem.Mol], Optional[float]]:
    try:
        query_smi = input("Optional query SMILES for similarity search (press Enter to skip): ").strip()
    except EOFError:
        return None, None
    if not query_smi:
        return None, None

    query_mol = prepare_mol_from_smiles(query_smi)
    if query_mol is None:
        print("Invalid query SMILES. Falling back to range selection.", flush=True)
        return None, None

    try:
        raw_threshold = input(
            f"Similarity percentage cutoff (0-100, default {int(DEFAULT_SIMILARITY_THRESHOLD * 100)}): "
        ).strip()
    except EOFError:
        raw_threshold = ""
    if not raw_threshold:
        return query_mol, DEFAULT_SIMILARITY_THRESHOLD

    try:
        pct = float(raw_threshold)
        if pct > 1:
            pct = pct / 100.0
        pct = max(0.0, min(1.0, pct))
    except ValueError:
        print("Invalid percentage. Using default threshold.", flush=True)
        pct = DEFAULT_SIMILARITY_THRESHOLD

    return query_mol, pct


def prompt_range_selection(max_index: int) -> List[int]:
    try:
        raw = input("Range to render (all, 1-25, 50-75, 1,5,10-20; Enter=all): ").strip()
    except EOFError:
        raw = ""
    if not raw:
        return list(range(1, max_index + 1))
    try:
        return parse_range_selection(raw, max_index)
    except Exception:
        print("Invalid range input. Using all molecules.", flush=True)
        return list(range(1, max_index + 1))


def page_highlight_info(entries: List[MolEntry]) -> List[Tuple[List[int], Dict[int, Tuple[float, float, float]]]]:
    n = len(entries)
    flags = [False] * n
    for i in range(n):
        for j in range(i + 1, n):
            sim = DataStructs.TanimotoSimilarity(entries[i].fp, entries[j].fp)
            if sim >= PAGE_LOCAL_SIMILARITY_THRESHOLD:
                flags[i] = True
                flags[j] = True

    info: List[Tuple[List[int], Dict[int, Tuple[float, float, float]]]] = []
    for i in range(n):
        if not flags[i]:
            info.append(([], {}))
            continue

        best_j = None
        best_sim = -1.0
        for j in range(n):
            if i == j:
                continue
            sim = DataStructs.TanimotoSimilarity(entries[i].fp, entries[j].fp)
            if sim > best_sim:
                best_sim = sim
                best_j = j

        bond_ids: List[int] = []
        atom_colors: Dict[int, Tuple[float, float, float]] = {}
        if best_j is not None and best_sim >= PAGE_LOCAL_SIMILARITY_THRESHOLD:
            try:
                res = rdFMCS.FindMCS(
                    [entries[i].mol, entries[best_j].mol],
                    timeout=2,
                    ringMatchesRingOnly=True,
                    completeRingsOnly=True,
                )
                if res and res.smartsString:
                    patt = Chem.MolFromSmarts(res.smartsString)
                    if patt is not None:
                        match = entries[i].mol.GetSubstructMatch(patt)
                        if match:
                            aset = set(match)
                            atom_colors = {a: HIGHLIGHT_COLOR for a in aset}
                            for b in entries[i].mol.GetBonds():
                                if b.GetBeginAtomIdx() in aset and b.GetEndAtomIdx() in aset:
                                    bond_ids.append(b.GetIdx())
            except Exception:
                pass

        info.append((bond_ids, atom_colors))
    return info


def _draw_with_options(
    drawer,
    mol: Chem.Mol,
    bond_ids: List[int],
    atom_colors: Dict[int, Tuple[float, float, float]],
):
    opts = drawer.drawOptions()
    opts.variableBondWidthMultiplier = VARIABLE_BOND_WIDTH_MULTIPLIER
    opts.addAtomIndices = False
    opts.addBondIndices = False
    opts.padding = 0.03
    if hasattr(opts, "minFontSize"):
        opts.minFontSize = 6
    if hasattr(opts, "maxFontSize"):
        opts.maxFontSize = 28
    if hasattr(opts, "bondLineWidth"):
        opts.bondLineWidth = 2.0
    if hasattr(opts, "useBWAtomPalette"):
        opts.useBWAtomPalette()

    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        mol,
        highlightAtoms=list(atom_colors.keys()),
        highlightAtomColors=atom_colors,
        highlightBonds=bond_ids,
        highlightBondColors={b: HIGHLIGHT_COLOR for b in bond_ids} if bond_ids else {},
    )
    drawer.FinishDrawing()


def render_mol_image(
    mol: Chem.Mol,
    bond_ids: List[int],
    atom_colors: Dict[int, Tuple[float, float, float]],
) -> Image.Image:
    if hasattr(rdMolDraw2D, "MolDraw2DCairo"):
        drawer = rdMolDraw2D.MolDraw2DCairo(DRAW_SIZE[0], DRAW_SIZE[1])
        _draw_with_options(drawer, mol, bond_ids, atom_colors)
        return Image.open(io.BytesIO(drawer.GetDrawingText())).convert("RGB")

    if hasattr(rdMolDraw2D, "MolDraw2DSVG"):
        drawer = rdMolDraw2D.MolDraw2DSVG(DRAW_SIZE[0], DRAW_SIZE[1])
        _draw_with_options(drawer, mol, bond_ids, atom_colors)
        svg = drawer.GetDrawingText()
        if cairosvg is None:
            raise RuntimeError(
                "Your RDKit build does not provide MolDraw2DCairo. Install cairosvg "
                "or use an RDKit build with Cairo support."
            )
        png_bytes = cairosvg.svg2png(bytestring=svg.encode("utf-8"))
        return Image.open(io.BytesIO(png_bytes)).convert("RGB")

    raise RuntimeError("This RDKit build provides neither MolDraw2DCairo nor MolDraw2DSVG.")


def draw_cell(
    c: canvas.Canvas,
    entry: MolEntry,
    label_prefix: str,
    x: float,
    y: float,
    w: float,
    h: float,
    bond_ids: List[int],
    atom_colors: Dict[int, Tuple[float, float, float]],
):
    c.setLineWidth(BORDER_WIDTH)
    c.setStrokeColor(colors.black)
    c.rect(x, y, w, h, stroke=1, fill=0)

    label = format_label(label_prefix, entry.original_index)
    label_y = y + h - TOP_PAD
    c.setFont(LABEL_FONT, LABEL_SIZE)
    c.drawCentredString(x + w / 2.0, label_y, label)

    text_lines = descriptor_lines(entry.mol)
    footer_bottom = y + PADDING
    c.setFont(TEXT_FONT, TEXT_SIZE)
    line_gap = 6.1
    first_line_y = footer_bottom + 2
    text_center_x = x + w / 2.0
    for i, line in enumerate(text_lines):
        yy = first_line_y + i * line_gap
        c.drawCentredString(text_center_x, yy, line)

    image_x = x + PADDING
    image_y = first_line_y + len(text_lines) * line_gap + FOOTER_TEXT_GAP
    image_w = w - 2 * PADDING
    image_h = h - TITLE_BLOCK_H - (len(text_lines) * line_gap + FOOTER_TEXT_GAP + PADDING) - 2 * PADDING

    img = render_mol_image(entry.mol, bond_ids, atom_colors)
    iw, ih = img.size
    scale = min(image_w / iw, image_h / ih)
    draw_w = iw * scale
    draw_h = ih * scale
    dx = image_x + (image_w - draw_w) / 2.0
    dy = image_y + (image_h - draw_h) / 2.0
    c.drawImage(ImageReader(img), dx, dy, width=draw_w, height=draw_h, preserveAspectRatio=True, mask="auto")


def output_pdf_name(input_path: str) -> str:
    p = Path(input_path)
    return str(p.with_suffix("")) + "_grid.pdf"


def output_csv_name(input_path: str) -> str:
    p = Path(input_path)
    return str(p.with_suffix("")) + CSV_SUFFIX


def resolve_selection(entries: Sequence[MolEntry], args) -> Tuple[List[MolEntry], str]:
    if args.query_smiles is not None:
        query_mol = prepare_mol_from_smiles(args.query_smiles)
        if query_mol is None:
            raise SystemExit("Query SMILES is invalid")
        threshold = args.similarity if args.similarity is not None else DEFAULT_SIMILARITY_THRESHOLD
        if threshold > 1:
            threshold = threshold / 100.0
        qfp = rdMolDescriptors.GetMorganFingerprintAsBitVect(query_mol, radius=2, nBits=2048)
        chosen = [e for e in entries if DataStructs.TanimotoSimilarity(qfp, e.fp) >= threshold]
        return chosen, f"similarity >= {threshold:.2f}"

    if args.range_spec is not None:
        indices = set(parse_range_selection(args.range_spec, len(entries)))
        chosen = [e for e in entries if e.original_index in indices]
        return chosen, "range/all selection"

    if not args.no_prompt:
        query_mol, threshold = prompt_similarity_mode()
        if query_mol is not None and threshold is not None:
            qfp = rdMolDescriptors.GetMorganFingerprintAsBitVect(query_mol, radius=2, nBits=2048)
            chosen = [e for e in entries if DataStructs.TanimotoSimilarity(qfp, e.fp) >= threshold]
            return chosen, f"similarity >= {threshold:.2f}"
        indices = set(prompt_range_selection(len(entries)))
        chosen = [e for e in entries if e.original_index in indices]
        return chosen, "range/all selection"

    indices = set(range(1, len(entries) + 1))
    chosen = [e for e in entries if e.original_index in indices]
    return chosen, "range/all selection"


def write_csv(entries: Sequence[MolEntry], input_json: str, out_csv: str, pdf_name: str, mode: str) -> None:
    import csv
    fieldnames = [
        "output_index",
        "original_index",
        "source_index",
        "label",
        "smiles",
        "can_smiles",
        "formula",
        "mw",
        "clogp",
        "hbd",
        "hba",
        "tpsa",
        "rot_bonds",
        "rings",
        "heavy_atoms",
        "fsp3",
        "aromatic_rings",
        "hetero_atoms",
        "input_json",
        "pdf_file",
        "selection_mode",
    ]
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for out_i, entry in enumerate(entries, start=1):
            row = {
                "output_index": out_i,
                "original_index": entry.original_index,
                "source_index": entry.source_index,
                "label": format_label(DEFAULT_LABEL_PREFIX, entry.original_index),
                "smiles": entry.smiles,
                "can_smiles": Chem.MolToSmiles(entry.mol),
                "input_json": input_json,
                "pdf_file": pdf_name,
                "selection_mode": mode,
            }
            row.update(entry.record)
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="Render SMILES from JSON to 5x5 boxed PDF pages")
    parser.add_argument("input_json")
    parser.add_argument("--label-prefix", default=DEFAULT_LABEL_PREFIX)
    parser.add_argument("--output", default=None)
    parser.add_argument("--query-smiles", default=None)
    parser.add_argument("--similarity", type=float, default=None, help="0-100 or 0-1")
    parser.add_argument("--range", dest="range_spec", default=None, help="Examples: all, 1-25, 50-75, 1,5,10-20")
    parser.add_argument("--prompt", action="store_true", help="Deprecated; prompts are now default")
    parser.add_argument("--no-prompt", action="store_true", help="Disable interactive prompts")
    parser.add_argument("--max-pages", type=int, default=0, help="Debug: render only first N pages")
    args = parser.parse_args()

    input_json = args.input_json
    label_prefix = args.label_prefix
    out_pdf = args.output or output_pdf_name(input_json)

    print(f"Reading input JSON: {input_json}", flush=True)
    t0 = time.time()
    smiles = extract_smiles_from_json(input_json)
    print(f"Extracted {len(smiles)} candidate SMILES in JSON order", flush=True)
    if not smiles:
        raise SystemExit("No valid SMILES strings found in input JSON")

    entries = build_entries(smiles)
    print(f"Prepared {len(entries)} valid molecules", flush=True)
    if not entries:
        raise SystemExit("No renderable molecules produced from extracted SMILES")

    selected_entries, mode = resolve_selection(entries, args)
    print(f"Selection mode: {mode}", flush=True)
    print(f"Selected {len(selected_entries)} molecules", flush=True)
    if not selected_entries:
        raise SystemExit("No molecules matched the requested selection")

    page_w, page_h = PAGE_SIZE
    usable_w = page_w - 2 * MARGIN - (GRID_COLS - 1) * CELL_GAP
    usable_h = page_h - 2 * MARGIN - (GRID_ROWS - 1) * CELL_GAP
    cell_w = usable_w / GRID_COLS
    cell_h = usable_h / GRID_ROWS

    c = canvas.Canvas(out_pdf, pagesize=PAGE_SIZE)
    pages = page_selection_chunks(selected_entries)
    if args.max_pages and args.max_pages > 0:
        pages = pages[: args.max_pages]
        print(f"Debug mode: rendering first {len(pages)} page(s)", flush=True)
    total_pages = len(pages)
    print(f"Rendering {total_pages} pages", flush=True)

    for page_num, page_entries in enumerate(pages, start=1):
        print(f"Rendering page {page_num}/{total_pages} ({len(page_entries)} molecules)...", flush=True)
        highlights = page_highlight_info(page_entries)
        for idx, entry in enumerate(page_entries):
            row = idx // GRID_COLS
            col = idx % GRID_COLS
            x = MARGIN + col * (cell_w + CELL_GAP)
            y = page_h - MARGIN - (row + 1) * cell_h - row * CELL_GAP
            bond_ids, atom_colors = highlights[idx]
            draw_cell(c, entry, label_prefix, x, y, cell_w, cell_h, bond_ids, atom_colors)
        c.showPage()
    c.save()

    out_csv = output_csv_name(input_json)
    write_csv(selected_entries, input_json, out_csv, out_pdf, mode)

    print(f"Wrote {out_pdf} with {len(selected_entries)} molecules in {time.time()-t0:.1f}s", flush=True)
    print(f"Wrote {out_csv}", flush=True)


if __name__ == "__main__":
    main()
