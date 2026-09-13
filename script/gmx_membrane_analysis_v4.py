#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GROMACS membrane structural analysis

Main modules:
1) Area per lipid (APL) vs Time
2) Membrane thickness vs Time
3) Membrane SASA vs Time
4) Lipid-chain order parameter (batch gmx order)
5) Lipid / molecule MSD / diffusion coefficient (batch gmx msd)
   - residue/molecule type(s)
   - whole-membrane structural regions
   - specific AA atom / CG bead
   - selected molecule/residue in dynamic membrane regions
6) Molecular / residue density distribution (gmx density)
7) Drug MSD / diffusion in dynamically identified membrane regions
8) RMSD / structural deviation (gmx rms)
   - default GROMACS groups
   - complete residue/molecule type(s)
   - specific AA atom / CG bead subset(s)
   - individual molecule/residue instance(s)
   - same or separate fitting group
9) Membrane permeation rate (%) + molecule count
   - dynamic membrane Z range from selected lipid markers
   - molecules inside membrane / below lower leaflet / both
   - molecule-count time series and mean +/- SD
   - user-entered initial loading N0
   - cumulative permeation rate (%) vs time and final permeability (%)

Designed for approximately planar bilayers:
- membrane plane = XY
- membrane normal = Z

Definitions
-----------
APL(t) = Lx(t) * Ly(t) / N_leaflet

Thickness(t) =
    | Z_COG(upper leaflet representative head atoms)
    - Z_COG(lower leaflet representative head atoms) |

SASA(t) =
    solvent-accessible surface area of the combined [ membrane ] group

Order parameter:
    one clean order index per lipid/chain, containing only groups of
    equivalent sequential chain atoms from that lipid residue; calculated
    with gmx order -d z -szonly

Usage
-----
python3 gmx_membrane_analysis_v20.py

The program starts with a MAIN MENU. Choose only the parameter you want:
1) APL calculation
2) Membrane thickness calculation
3) Membrane SASA calculation
4) Lipid order parameter calculation
5) Lipid MSD / diffusion coefficient calculation
6) Molecular / residue density distribution calculation
7) Drug MSD / diffusion by dynamic membrane region
8) RMSD / structural deviation calculation
9) Membrane permeation rate (%) / molecule-count analysis
0) Exit

After every completed calculation you can:
1) continue calculating the current parameter
2) return to the main menu and choose another parameter
3) exit the program

If your GROMACS command is gmx_mpi:
GMX=gmx_mpi python3 gmx_membrane_analysis_v20.py
"""

from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import zipfile
from xml.sax.saxutils import escape as xml_escape
import threading
import time
import traceback
from collections import Counter
from datetime import datetime
from pathlib import Path
from statistics import mean, stdev, median


GMX = os.environ.get("GMX", "gmx")

OUTDIR = Path("Membrane_analysis")
APL_DIR = OUTDIR / "01_APL"
THICK_DIR = OUTDIR / "02_Thickness"
SASA_DIR = OUTDIR / "03_SASA"
ORDER_DIR = OUTDIR / "04_Order_Parameter"
MSD_DIR = OUTDIR / "05_MSD_Diffusion"
DENSITY_DIR = OUTDIR / "06_Density_Distribution"
DRUG_REGION_DIR = OUTDIR / "07_Drug_Dynamic_Region_Diffusion"
RMSD_DIR = OUTDIR / "08_RMSD"
PERMEATION_DIR = OUTDIR / "09_Permeation_Count"

INDEX = Path("index.ndx")
BASE_NDX = OUTDIR / "default_index.ndx"
CUSTOM_NDX = OUTDIR / "custom_membrane_index.ndx"
LIPID_CONFIG = OUTDIR / "lipid_headgroups.tsv"
SYSTEM_SUMMARY = OUTDIR / "system_summary.txt"
REPORT = OUTDIR / "analysis_report.txt"



def log_timestamp() -> str:
    """
    Timestamp used in ALL log filenames.

    Format:
        YYYYMMDD_HHMMSS_microseconds

    Microseconds avoid accidental collisions when several GROMACS commands
    finish within the same second.
    """
    return datetime.now().strftime(
        "%Y%m%d_%H%M%S_%f"
    )


def log_timestamp_text() -> str:
    """
    Human-readable local system timestamp written inside every log file.
    """
    return datetime.now().astimezone().isoformat(
        timespec="seconds"
    )


def timestamped_log_path(
    base_path: Path,
) -> Path:
    """
    Return a path with a timestamp inserted before the suffix.

    Example:
        msd_DEF.log
            ->
        msd_DEF_20260912_213045_123456.log
    """
    base_path = Path(
        base_path
    )

    suffix = (
        base_path.suffix
        if base_path.suffix
        else ".log"
    )

    stem = (
        base_path.stem
        if base_path.suffix
        else base_path.name
    )

    return base_path.with_name(
        f"{stem}_{log_timestamp()}{suffix}"
    )


def log_timestamp_header() -> str:
    """
    Standard header inserted into command/error/result log files.
    """
    return (
        "=" * 78
        + "\n"
        + f"LOG CREATED AT: {log_timestamp_text()}\n"
        + "=" * 78
        + "\n\n"
    )


def ask_positive_integer(
    prompt: str,
    default: int | None = None,
) -> int:
    """
    Positive integer input that stays at the current prompt after bad input.
    """
    while True:
        raw = ask(
            prompt,
            str(default)
            if default is not None
            else None,
        )

        try:
            value = int(
                raw
            )
        except ValueError:
            print(
                "[WARNING] Please enter a positive integer."
            )
            continue

        if value <= 0:
            print(
                "[WARNING] The value must be greater than 0."
            )
            continue

        return value


def _xlsx_col_name(number: int) -> str:
    """Convert a 1-based Excel column number to A, B, ..., AA."""
    name = ""
    while number:
        number, rem = divmod(number - 1, 26)
        name = chr(65 + rem) + name
    return name


def _xlsx_safe_text(value) -> str:
    """Remove XML-invalid characters and escape text."""
    value = "" if value is None else str(value)
    value = re.sub(r"[\x00-\x08\x0B\x0C\x0E-\x1F]", "", value)
    return xml_escape(value)


def _xlsx_is_number(value) -> bool:
    """Return True when a value should be written as a numeric Excel cell."""
    if value is None or isinstance(value, bool):
        return False
    if isinstance(value, (int, float)):
        return True
    if isinstance(value, str):
        s = value.strip()
        if not s or s.upper() in {"NA", "N/A", "NAN", "NOT_FOUND"}:
            return False
        try:
            float(s)
            return True
        except ValueError:
            return False
    return False


def write_excel_table(
    xlsx_file: Path,
    sheet_name: str,
    headers: list[str],
    rows: list[list],
):
    """
    Create a real .xlsx workbook using only Python's standard library.

    The output includes:
      - a formatted header row
      - frozen first row
      - AutoFilter
      - bounded, readable column widths
      - numeric values stored as numeric cells

    Compatible with Microsoft Excel, WPS Office and LibreOffice Calc.
    """
    xlsx_file.parent.mkdir(parents=True, exist_ok=True)

    safe_sheet = re.sub(r'[:\\/?*\[\]]', "_", sheet_name)[:31] or "Summary"

    all_rows = [headers] + rows
    ncols = len(headers)
    nrows = len(all_rows)

    widths = []
    for col_idx in range(ncols):
        longest = len(str(headers[col_idx]))
        for row in rows:
            if col_idx < len(row):
                cell = "" if row[col_idx] is None else str(row[col_idx])
                longest = max(longest, len(cell))
        widths.append(min(max(longest + 2, 10), 42))

    cols_xml = "".join(
        f'<col min="{i}" max="{i}" width="{width}" customWidth="1"/>'
        for i, width in enumerate(widths, start=1)
    )

    row_xml_parts = []
    for row_idx, row in enumerate(all_rows, start=1):
        cell_xml = []

        for col_idx in range(1, ncols + 1):
            value = row[col_idx - 1] if col_idx - 1 < len(row) else ""
            ref = f"{_xlsx_col_name(col_idx)}{row_idx}"
            style = ' s="1"' if row_idx == 1 else ""

            if row_idx > 1 and _xlsx_is_number(value):
                numeric = float(value)
                numeric_text = (
                    str(int(numeric))
                    if numeric.is_integer()
                    else repr(numeric)
                )
                cell_xml.append(
                    f'<c r="{ref}"{style}><v>{numeric_text}</v></c>'
                )
            else:
                txt = _xlsx_safe_text(value)
                cell_xml.append(
                    f'<c r="{ref}" t="inlineStr"{style}>'
                    f'<is><t xml:space="preserve">{txt}</t></is></c>'
                )

        height = ' ht="22" customHeight="1"' if row_idx == 1 else ""
        row_xml_parts.append(
            f'<row r="{row_idx}"{height}>'
            + "".join(cell_xml)
            + "</row>"
        )

    last_cell = f"{_xlsx_col_name(ncols)}{max(nrows, 1)}"
    autofilter = (
        f'<autoFilter ref="A1:{last_cell}"/>'
        if ncols > 0
        else ""
    )

    sheet_xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <sheetViews>
    <sheetView workbookViewId="0">
      <pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>
    </sheetView>
  </sheetViews>
  <sheetFormatPr defaultRowHeight="18"/>
  <cols>{cols_xml}</cols>
  <sheetData>{''.join(row_xml_parts)}</sheetData>
  {autofilter}
</worksheet>"""

    workbook_xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
          xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <sheets>
    <sheet name="{_xlsx_safe_text(safe_sheet)}" sheetId="1" r:id="rId1"/>
  </sheets>
</workbook>"""

    workbook_rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1"
    Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet"
    Target="worksheets/sheet1.xml"/>
  <Relationship Id="rId2"
    Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles"
    Target="styles.xml"/>
</Relationships>"""

    root_rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1"
    Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument"
    Target="xl/workbook.xml"/>
</Relationships>"""

    content_types = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels"
    ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/xl/workbook.xml"
    ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
  <Override PartName="/xl/worksheets/sheet1.xml"
    ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
  <Override PartName="/xl/styles.xml"
    ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
</Types>"""

    styles_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <fonts count="2">
    <font><sz val="11"/><name val="Calibri"/></font>
    <font><b/><color rgb="FFFFFFFF"/><sz val="11"/><name val="Calibri"/></font>
  </fonts>
  <fills count="3">
    <fill><patternFill patternType="none"/></fill>
    <fill><patternFill patternType="gray125"/></fill>
    <fill>
      <patternFill patternType="solid">
        <fgColor rgb="FF1F4E78"/>
        <bgColor indexed="64"/>
      </patternFill>
    </fill>
  </fills>
  <borders count="1">
    <border><left/><right/><top/><bottom/><diagonal/></border>
  </borders>
  <cellStyleXfs count="1">
    <xf numFmtId="0" fontId="0" fillId="0" borderId="0"/>
  </cellStyleXfs>
  <cellXfs count="2">
    <xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>
    <xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0"
        applyFont="1" applyFill="1" applyAlignment="1">
      <alignment horizontal="center" vertical="center"/>
    </xf>
  </cellXfs>
  <cellStyles count="1">
    <cellStyle name="Normal" xfId="0" builtinId="0"/>
  </cellStyles>
</styleSheet>"""

    with zipfile.ZipFile(
        xlsx_file,
        "w",
        compression=zipfile.ZIP_DEFLATED,
    ) as zf:
        zf.writestr("[Content_Types].xml", content_types)
        zf.writestr("_rels/.rels", root_rels)
        zf.writestr("xl/workbook.xml", workbook_xml)
        zf.writestr("xl/_rels/workbook.xml.rels", workbook_rels)
        zf.writestr("xl/styles.xml", styles_xml)
        zf.writestr("xl/worksheets/sheet1.xml", sheet_xml)

    print(f"[OK] Excel summary generated: {xlsx_file}")


def append_excel_summary(
    xlsx_file: Path,
    cache_file: Path,
    sheet_name: str,
    headers: list[str],
    row: list,
):
    """
    Append a record to a persistent Excel summary.

    A JSON cache is used internally so this script does not depend on
    openpyxl, pandas or xlsxwriter. The user-facing result is .xlsx.
    """
    records = []

    if cache_file.exists():
        try:
            loaded = json.loads(
                cache_file.read_text(encoding="utf-8")
            )
            if isinstance(loaded, list):
                records = loaded
        except Exception:
            print(
                f"[WARNING] Could not read {cache_file}; "
                "a new summary cache will be created."
            )

    records.append(row)

    cache_file.write_text(
        json.dumps(records, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    write_excel_table(
        xlsx_file,
        sheet_name,
        headers,
        records,
    )


class RecoverableAnalysisError(Exception):
    """Stop only the current calculation attempt, not the whole program."""
    pass


def die(msg: str):
    """
    Convert fatal-style helper calls into recoverable analysis errors.
    The interactive controller decides whether to retry, return to menu,
    or exit only when the user explicitly chooses to exit.
    """
    raise RecoverableAnalysisError(msg)


def ask(prompt: str, default: str | None = None) -> str:
    if default is None:
        return input(prompt).strip()
    ans = input(f"{prompt} [{default}]: ").strip()
    return ans if ans else default


def ask_float(prompt: str, default: float | None = None) -> float:
    while True:
        raw = ask(prompt, str(default) if default is not None else None)
        try:
            value = float(raw)
            if value < 0:
                raise ValueError
            return value
        except ValueError:
            print("[WARNING] Please enter a non-negative number.")



def ask_optional_float(prompt: str) -> float | None:
    """
    Optional non-negative float. Bad input stays at THIS prompt.
    """
    while True:
        raw = input(f"{prompt} [trajectory end]: ").strip()

        if not raw:
            return None

        try:
            value = float(raw)
            if value < 0:
                raise ValueError
            return value
        except ValueError:
            print(
                "[WARNING] Enter a non-negative number, "
                "or press Enter to use the trajectory end."
            )




def _animated_process(
    cmd,
    input_text: str | None = None,
    label: str = "GROMACS calculation",
):
    """
    Run a subprocess with an activity bar.

    100% DONE is shown only for return code 0. A failed command displays
    FAILED. Ctrl+C terminates the current subprocess and becomes a
    recoverable error instead of closing the analysis program.
    """
    try:
        proc = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE if input_text is not None else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError as exc:
        raise RecoverableAnalysisError(
            f"Could not start command:\\n  {' '.join(map(str, cmd))}\\n"
            f"System error: {exc}"
        ) from exc

    holder = {}

    def worker():
        stdout, stderr = proc.communicate(input=input_text)
        holder["stdout"] = stdout or ""
        holder["stderr"] = stderr or ""

    thread = threading.Thread(target=worker, daemon=True)
    thread.start()

    width = 28
    block = 6
    pos = 0
    direction = 1
    t0 = time.time()

    try:
        while thread.is_alive():
            chars = [" "] * width
            for j in range(block):
                k = pos + j
                if 0 <= k < width:
                    chars[k] = "█"

            elapsed = time.time() - t0
            print(
                f"\r[{''.join(chars)}] {label} | running {elapsed:6.1f} s",
                end="",
                flush=True,
            )

            pos += direction
            if pos + block >= width:
                direction = -1
            elif pos <= 0:
                direction = 1

            time.sleep(0.15)

    except KeyboardInterrupt:
        try:
            proc.terminate()
            proc.wait(timeout=5)
        except Exception:
            try:
                proc.kill()
            except Exception:
                pass

        print(
            f"\r[{'!' * width}] {label} | INTERRUPTED"
            + " " * 16
        )
        raise RecoverableAnalysisError(
            f"The current calculation '{label}' was interrupted. "
            "You can correct the settings and retry."
        )

    thread.join()
    elapsed = time.time() - t0

    if proc.returncode == 0:
        print(
            f"\r[{'█' * width}] {label} | 100% DONE ({elapsed:.1f} s)"
            + " " * 12
        )
    else:
        print(
            f"\r[{'!' * width}] {label} | FAILED "
            f"(exit code {proc.returncode}, {elapsed:.1f} s)"
            + " " * 12
        )

    return subprocess.CompletedProcess(
        cmd,
        proc.returncode,
        holder.get("stdout", ""),
        holder.get("stderr", ""),
    )




def run(
    cmd,
    input_text: str | None = None,
    quiet: bool = False,
    progress_label: str | None = None,
):
    """Run one external command without terminating the whole program."""
    try:
        if progress_label:
            result = _animated_process(
                cmd,
                input_text=input_text,
                label=progress_label,
            )
        elif quiet:
            result = subprocess.run(
                cmd,
                input=input_text,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        else:
            result = subprocess.run(
                cmd,
                input=input_text,
                text=True,
            )

    except RecoverableAnalysisError:
        raise
    except OSError as exc:
        raise RecoverableAnalysisError(
            f"Could not execute command:\\n  {' '.join(map(str, cmd))}\\n"
            f"System error: {exc}"
        ) from exc

    if result.returncode != 0:
        if quiet or progress_label:
            if result.stdout:
                print(result.stdout)
            if result.stderr:
                print(result.stderr, file=sys.stderr)

        raise RecoverableAnalysisError(
            "GROMACS command failed:\\n  "
            + " ".join(map(str, cmd))
            + f"\\nExit code: {result.returncode}"
        )

    return result




def check_environment() -> bool:
    """
    Check the GROMACS executable without aborting on a wrong command name.
    """
    global GMX

    while shutil.which(GMX) is None:
        print()
        print("=" * 78)
        print(" GROMACS COMMAND NOT FOUND")
        print("=" * 78)
        print(f"The command '{GMX}' is not available in PATH.")
        print("Common examples: gmx, gmx_mpi")
        print()

        candidate = ask(
            "Enter another GROMACS command, or 0 to exit",
            "gmx_mpi",
        ).strip()

        if candidate == "0":
            return False

        if not candidate:
            print("[WARNING] Please enter a command name.")
            continue

        GMX = candidate

    print(f"[OK] GROMACS command detected: {GMX}")
    return True



def parse_gro(gro_file: Path, lipid_defs: list[tuple[str, str]]):
    """
    Validate lipid residue names and representative head atom names.

    Returns:
      total_lipids
      res_counts
      head_counts
      head_records = list[(atom_index, resid, resname, atomname, z)]
      zmid
      upper_atom_indices
      lower_atom_indices
      initial thickness
      box_z
    """
    lines = gro_file.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 3:
        die("Invalid GRO file.")

    try:
        natoms = int(lines[1].strip())
    except Exception:
        die("Could not read atom count from GRO.")

    atom_lines = lines[2:2 + natoms]
    if len(atom_lines) != natoms:
        die("GRO atom count does not match file content.")

    try:
        box_values = [float(x) for x in lines[2 + natoms].split()]
    except Exception:
        die("Could not read GRO box vectors.")

    box_z = box_values[2] if len(box_values) >= 3 else 0.0

    wanted = dict(lipid_defs)
    lipid_names = set(wanted)

    res_counts = Counter()
    head_counts = Counter()
    head_records = []

    previous_res_block = None

    for atom_index, line in enumerate(atom_lines, start=1):
        if len(line) < 44:
            continue

        resnr = line[0:5].strip()
        resname = line[5:10].strip()
        atomname = line[10:15].strip()

        try:
            z = float(line[36:44])
        except ValueError:
            continue

        res_block = line[0:10]
        if res_block != previous_res_block:
            if resname in lipid_names:
                res_counts[resname] += 1
            previous_res_block = res_block

        if resname in wanted and atomname == wanted[resname]:
            head_counts[resname] += 1
            head_records.append(
                (atom_index, resnr, resname, atomname, z)
            )

    errors = []
    for resname, head in lipid_defs:
        nres = res_counts[resname]
        nhead = head_counts[resname]

        if nres == 0:
            errors.append(
                f"Residue '{resname}' was not found in GRO. "
                f"Please use the actual residue name (resname)."
            )
        elif nhead == 0:
            errors.append(
                f"Head atom/bead '{head}' was not found in residue '{resname}'."
            )
        elif nhead != nres:
            errors.append(
                f"{resname}: detected {nres} lipid residues but {nhead} "
                f"atoms/beads named '{head}'. This method expects exactly "
                f"ONE representative head atom/bead per lipid."
            )

    if errors:
        raise RecoverableAnalysisError(
            "GRO/headgroup validation failed:\n  - "
            + "\n  - ".join(errors)
        )

    if len(head_records) < 2:
        die("Too few headgroup reference atoms to identify two leaflets.")

    z_values = [x[4] for x in head_records]

    # 1D k-means, k=2
    c1, c2 = min(z_values), max(z_values)

    for _ in range(200):
        g1, g2 = [], []
        for z in z_values:
            if abs(z - c1) <= abs(z - c2):
                g1.append(z)
            else:
                g2.append(z)

        if not g1 or not g2:
            die(
                "Could not split headgroups into two leaflets. "
                "The bilayer may be broken across Z PBC."
            )

        n1 = sum(g1) / len(g1)
        n2 = sum(g2) / len(g2)

        if abs(n1 - c1) < 1e-10 and abs(n2 - c2) < 1e-10:
            c1, c2 = n1, n2
            break

        c1, c2 = n1, n2

    lower_center = min(c1, c2)
    upper_center = max(c1, c2)
    zmid = 0.5 * (lower_center + upper_center)

    upper_indices = [x[0] for x in head_records if x[4] >= zmid]
    lower_indices = [x[0] for x in head_records if x[4] < zmid]

    total_lipids = sum(res_counts[x[0]] for x in lipid_defs)
    initial_thickness = upper_center - lower_center

    return {
        "total_lipids": total_lipids,
        "res_counts": res_counts,
        "head_counts": head_counts,
        "head_records": head_records,
        "zmid": zmid,
        "upper_indices": upper_indices,
        "lower_indices": lower_indices,
        "initial_thickness": initial_thickness,
        "box_z": box_z,
        "upper_center": upper_center,
        "lower_center": lower_center,
    }


def write_ndx_group(fh, name: str, indices: list[int], per_line: int = 15):
    fh.write(f"[ {name} ]\n")
    for i in range(0, len(indices), per_line):
        chunk = indices[i:i + per_line]
        fh.write(" ".join(str(x) for x in chunk) + "\n")
    fh.write("\n")


def get_membrane_atom_indices(gro_file: Path, lipid_names: set[str]):
    lines = gro_file.read_text(encoding="utf-8", errors="replace").splitlines()
    natoms = int(lines[1].strip())
    indices = []

    for atom_index, line in enumerate(lines[2:2 + natoms], start=1):
        if len(line) < 20:
            continue
        resname = line[5:10].strip()
        if resname in lipid_names:
            indices.append(atom_index)

    return indices


def build_index(
    gro_file: Path,
    lipid_defs: list[tuple[str, str]],
    parsed: dict,
):
    print("\n---------------------------- INDEX SETUP -------------------------------")
    print(f"[INFO] Generating default GROMACS index from {gro_file} ...")

    run(
        [GMX, "make_ndx", "-f", str(gro_file), "-o", str(BASE_NDX)],
        input_text="q\n",
        quiet=True,
    )

    if not BASE_NDX.exists() or BASE_NDX.stat().st_size == 0:
        die("Default index was not generated.")

    lipid_names = {x[0] for x in lipid_defs}
    membrane_indices = get_membrane_atom_indices(gro_file, lipid_names)
    head_indices = [x[0] for x in parsed["head_records"]]
    upper_indices = parsed["upper_indices"]
    lower_indices = parsed["lower_indices"]

    with CUSTOM_NDX.open("w", encoding="utf-8") as fh:
        write_ndx_group(fh, "membrane", membrane_indices)
        write_ndx_group(fh, "membrane_heads", head_indices)
        write_ndx_group(fh, "upper_heads", upper_indices)
        write_ndx_group(fh, "lower_heads", lower_indices)

    if INDEX.exists():
        backup = INDEX.with_name(
            f"{INDEX.name}.bak"
        )
        shutil.copy2(INDEX, backup)
        print(f"[INFO] Existing {INDEX} backed up as {backup}")

    with INDEX.open("wb") as out:
        out.write(BASE_NDX.read_bytes())
        if not BASE_NDX.read_bytes().endswith(b"\n"):
            out.write(b"\n")
        out.write(CUSTOM_NDX.read_bytes())

    shutil.copy2(INDEX, OUTDIR / "index.ndx")

    print("[OK] Final index generated: index.ndx")
    print("[OK] Added groups:")
    print("     [ membrane ]")
    print("     [ membrane_heads ]")
    print("     [ upper_heads ]")
    print("     [ lower_heads ]")


def parse_xvg(path: Path, min_cols: int = 2):
    rows = []
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#") or s.startswith("@"):
                continue
            parts = s.split()
            if len(parts) < min_cols:
                continue
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    return rows


def write_stats(
    xvg_file: Path,
    out_file: Path,
    label: str,
    unit: str,
):
    rows = parse_xvg(xvg_file, min_cols=2)
    if not rows:
        die(f"No data in {xvg_file}")

    t = [r[0] for r in rows]
    values = [r[1] for r in rows]

    avg = mean(values)
    sd = stdev(values) if len(values) > 1 else 0.0

    text = (
        "============================================================\n"
        f"{label}\n"
        "============================================================\n"
        f"Time range     : {t[0]:.6f} - {t[-1]:.6f} ns\n"
        f"Data points    : {len(values)}\n"
        f"Mean           : {avg:.8f} {unit}\n"
        f"SD             : {sd:.8f} {unit}\n"
        f"Mean +/- SD    : {avg:.8f} +/- {sd:.8f} {unit}\n"
        f"Minimum        : {min(values):.8f} {unit}\n"
        f"Maximum        : {max(values):.8f} {unit}\n"
        "============================================================\n"
    )
    out_file.write_text(text, encoding="utf-8")
    print(text)


def gmx_time_args(begin_ps: float, end_ps: float | None, dt_ps: float | None = None):
    args = ["-b", str(begin_ps)]
    if end_ps is not None:
        args += ["-e", str(end_ps)]
    if dt_ps is not None and dt_ps > 0:
        args += ["-dt", str(dt_ps)]
    return args


def analyze_apl(
    edr: Path,
    nleaflet: float,
    begin_ps: float,
    end_ps: float | None,
):
    print("\n======================================================================")
    print(" [1/5] AREA PER LIPID (APL) vs Time")
    print("======================================================================")

    raw = APL_DIR / "box_xy_raw.xvg"
    out = APL_DIR / "apl_vs_time.xvg"
    csv_out = APL_DIR / "apl_vs_time.csv"
    stats_out = APL_DIR / "apl_stats.txt"

    cmd = [
        GMX, "energy",
        "-f", str(edr),
        "-o", str(raw),
        *gmx_time_args(begin_ps, end_ps),
    ]

    run(
        cmd,
        input_text="Box-X\nBox-Y\n0\n",
        quiet=True,
        progress_label="APL: extracting Box-X / Box-Y",
    )

    rows = parse_xvg(raw, min_cols=3)
    if not rows:
        die(
            "Could not extract Box-X/Box-Y. "
            f"Run '{GMX} energy -f {edr}' manually and inspect available terms."
        )

    with out.open("w", encoding="utf-8") as fh:
        fh.write('# Global projected area per lipid\n')
        fh.write('# APL(t) = Box-X(t) * Box-Y(t) / N_leaflet\n')
        fh.write(f'# N_leaflet = {nleaflet}\n')
        fh.write('@ title "Area per lipid vs Time"\n')
        fh.write('@ xaxis label "Time (ns)"\n')
        fh.write('@ yaxis label "Area per lipid (nm\\S2\\N/lipid)"\n')
        fh.write('@ s0 legend "APL"\n')

        for r in rows:
            time_ns = r[0] / 1000.0
            apl = (r[1] * r[2]) / nleaflet
            fh.write(f"{time_ns:.8f} {apl:.10f}\n")

    with csv_out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Time_ns", "APL_nm2_per_lipid"])
        for r in parse_xvg(out):
            writer.writerow([f"{r[0]:.8f}", f"{r[1]:.10f}"])

    write_stats(
        out,
        stats_out,
        "Area per lipid (APL)",
        "nm^2/lipid",
    )

    return out, csv_out, stats_out


def analyze_thickness(
    tpr: Path,
    xtc: Path,
    begin_ps: float,
    end_ps: float | None,
    dt_ps: float,
):
    print("\n======================================================================")
    print(" [2/5] MEMBRANE THICKNESS vs Time")
    print("======================================================================")

    upper_raw = THICK_DIR / "upper_heads_cog_z_raw.xvg"
    lower_raw = THICK_DIR / "lower_heads_cog_z_raw.xvg"
    out = THICK_DIR / "membrane_thickness_vs_time.xvg"
    csv_out = THICK_DIR / "membrane_thickness_vs_time.csv"
    stats_out = THICK_DIR / "membrane_thickness_stats.txt"

    time_args = gmx_time_args(begin_ps, end_ps, dt_ps)

    run(
        [
            GMX, "trajectory",
            "-f", str(xtc),
            "-s", str(tpr),
            "-n", str(INDEX),
            "-ox", str(upper_raw),
            *time_args,
            "-select", 'cog of group "upper_heads"',
            "-nox", "-noy", "-z",
        ],
        quiet=True,
        progress_label="Thickness: upper leaflet trajectory",
    )

    run(
        [
            GMX, "trajectory",
            "-f", str(xtc),
            "-s", str(tpr),
            "-n", str(INDEX),
            "-ox", str(lower_raw),
            *time_args,
            "-select", 'cog of group "lower_heads"',
            "-nox", "-noy", "-z",
        ],
        quiet=True,
        progress_label="Thickness: lower leaflet trajectory",
    )

    upper = parse_xvg(upper_raw, min_cols=2)
    lower = parse_xvg(lower_raw, min_cols=2)

    if not upper or not lower:
        die("Thickness calculation returned empty upper/lower leaflet data.")

    if len(upper) != len(lower):
        die("Upper and lower leaflet trajectories contain different frame counts.")

    with out.open("w", encoding="utf-8") as fh:
        fh.write('# Membrane thickness from upper/lower headgroup COG separation\n')
        fh.write('@ title "Membrane thickness vs Time"\n')
        fh.write('@ xaxis label "Time (ns)"\n')
        fh.write('@ yaxis label "Membrane thickness (nm)"\n')
        fh.write('@ s0 legend "Thickness"\n')

        for u, l in zip(upper, lower):
            if abs(u[0] - l[0]) > 1e-5:
                die("Upper/lower thickness time points do not match.")
            thickness = abs(u[1] - l[1])
            fh.write(f"{u[0]/1000.0:.8f} {thickness:.10f}\n")

    with csv_out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Time_ns", "Membrane_thickness_nm"])
        for r in parse_xvg(out):
            writer.writerow([f"{r[0]:.8f}", f"{r[1]:.10f}"])

    write_stats(
        out,
        stats_out,
        "Membrane thickness",
        "nm",
    )

    return out, csv_out, stats_out


def analyze_sasa(
    tpr: Path,
    xtc: Path,
    begin_ps: float,
    end_ps: float | None,
    dt_ps: float,
    probe: float,
    ndots: int,
):
    print("\n======================================================================")
    print(" [3/5] MEMBRANE SASA vs Time")
    print("======================================================================")

    raw = SASA_DIR / "membrane_sasa_raw.xvg"
    out = SASA_DIR / "membrane_sasa_vs_time.xvg"
    csv_out = SASA_DIR / "membrane_sasa_vs_time.csv"
    stats_out = SASA_DIR / "membrane_sasa_stats.txt"

    run(
        [
            GMX, "sasa",
            "-f", str(xtc),
            "-s", str(tpr),
            "-n", str(INDEX),
            "-o", str(raw),
            *gmx_time_args(begin_ps, end_ps, dt_ps),
            "-surface", 'group "membrane"',
            "-probe", str(probe),
            "-ndots", str(ndots),
        ],
        quiet=True,
        progress_label="SASA: membrane surface calculation",
    )

    rows = parse_xvg(raw, min_cols=2)
    if not rows:
        die("SASA output is empty.")

    with out.open("w", encoding="utf-8") as fh:
        fh.write('# Solvent-accessible surface area of combined membrane group\n')
        fh.write('@ title "Membrane SASA vs Time"\n')
        fh.write('@ xaxis label "Time (ns)"\n')
        fh.write('@ yaxis label "SASA (nm\\S2\\N)"\n')
        fh.write('@ s0 legend "Membrane SASA"\n')
        for r in rows:
            fh.write(f"{r[0]/1000.0:.8f} {r[1]:.10f}\n")

    with csv_out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Time_ns", "Membrane_SASA_nm2"])
        for r in parse_xvg(out):
            writer.writerow([f"{r[0]:.8f}", f"{r[1]:.10f}"])

    write_stats(
        out,
        stats_out,
        "Membrane solvent-accessible surface area (SASA)",
        "nm^2",
    )

    return out, csv_out, stats_out



def read_gro_residue_blocks(gro_file: Path):
    """Return contiguous GRO residue blocks with global 1-based atom indices."""
    lines = gro_file.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 3:
        die("Invalid GRO file.")
    try:
        natoms = int(lines[1].strip())
    except ValueError:
        die("Could not read atom count from GRO.")

    atom_lines = lines[2:2 + natoms]
    if len(atom_lines) != natoms:
        die("GRO atom count does not match file content.")

    blocks = []
    current_key = None
    current = None

    for atom_index, line in enumerate(atom_lines, start=1):
        if len(line) < 20:
            continue
        resnr = line[0:5].strip()
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        key = line[0:10]

        if key != current_key:
            current = {
                "resnr": resnr,
                "resname": resname,
                "atoms": [],
            }
            blocks.append(current)
            current_key = key

        current["atoms"].append((atom_index, atomname))

    return blocks


def lipid_residue_blocks(gro_file: Path, resname: str):
    return [b for b in read_gro_residue_blocks(gro_file) if b["resname"] == resname]


def show_first_lipid_atom_names(gro_file: Path, resname: str):
    blocks = lipid_residue_blocks(gro_file, resname)
    if not blocks:
        die(f"Residue '{resname}' was not found in {gro_file}.")

    names = [name for _, name in blocks[0]["atoms"]]
    print(f"[INFO] First {resname} residue contains {len(names)} atoms/beads in GRO order:")
    width = 8
    for i in range(0, len(names), width):
        chunk = names[i:i + width]
        print("       " + "  ".join(f"{x:>6s}" for x in chunk))
    return blocks


def collect_equivalent_atom_groups(
    gro_file: Path,
    resname: str,
    atom_names: list[str],
):
    """
    For each requested atom name, collect exactly one atom from every residue
    having the selected resname. This prevents same-named atoms from other
    lipid/molecule types entering the order index.
    """
    blocks = lipid_residue_blocks(gro_file, resname)
    if not blocks:
        die(f"Residue '{resname}' was not found in {gro_file}.")

    groups = []
    problems = []

    for atom_name in atom_names:
        indices = []
        for mol_i, block in enumerate(blocks, start=1):
            matches = [idx for idx, name in block["atoms"] if name == atom_name]
            if len(matches) != 1:
                problems.append(
                    f"{resname} molecule #{mol_i}: atom/bead name '{atom_name}' "
                    f"matched {len(matches)} times; exactly 1 is required."
                )
            else:
                indices.append(matches[0])
        groups.append((atom_name, indices))

    if problems:
        print("\n[ERROR] Order-index validation failed:", file=sys.stderr)
        for msg in problems[:20]:
            print("  - " + msg, file=sys.stderr)
        if len(problems) > 20:
            print(f"  ... and {len(problems)-20} more errors", file=sys.stderr)
        die(
            "Every selected chain atom/bead name must occur exactly once in "
            "every molecule of the selected lipid residue."
        )

    return groups, len(blocks)



def numbered_atom_bead_items(
    gro_file: Path,
    resname: str,
) -> list[str]:
    """
    Return unique atom/bead names from the first occurrence of one resname.

    The list preserves GRO order. It is used throughout membrane-related
    modules so users select by number instead of repeatedly typing names.
    """
    blocks = lipid_residue_blocks(
        gro_file,
        resname,
    )

    if not blocks:
        raise RecoverableAnalysisError(
            f"Residue/molecule type '{resname}' was not found in {gro_file}."
        )

    names = []
    seen = set()

    for _, atomname in blocks[0]["atoms"]:
        if atomname not in seen:
            seen.add(atomname)
            names.append(atomname)

    if not names:
        raise RecoverableAnalysisError(
            f"No atom/bead names could be read for '{resname}'."
        )

    return names


def show_numbered_atom_bead_items(
    gro_file: Path,
    resname: str,
    title: str | None = None,
) -> list[str]:
    """
    Show atom/bead choices in the requested compact style, for example:

        [1] CHOL-ROH   [2] CHOL-C1   [3] CHOL-C2

    This works identically for AA atoms and CG beads.
    """
    names = numbered_atom_bead_items(
        gro_file,
        resname,
    )

    print("\n" + "-" * 78)
    print(
        " "
        + (
            title
            if title is not None
            else f"ATOM / BEAD SELECTION — {resname}"
        )
    )
    print("-" * 78)

    cells = [
        f"[{i:>3d}] {resname}-{name}"
        for i, name in enumerate(
            names,
            start=1,
        )
    ]

    # Keep the terminal compact but readable.
    max_len = max(
        len(cell)
        for cell in cells
    )
    per_line = max(
        1,
        min(
            4,
            76 // max(
                max_len + 3,
                18,
            ),
        ),
    )

    for i in range(
        0,
        len(cells),
        per_line,
    ):
        print(
            "  "
            + "   ".join(
                f"{cell:<{max_len}s}"
                for cell in cells[i:i + per_line]
            )
        )

    print("-" * 78)
    return names


def prompt_numbered_atom_bead_selection(
    gro_file: Path,
    resname: str,
    *,
    title: str,
    instruction: str,
    allow_multiple: bool,
    min_count: int = 1,
    tpr: Path | None = None,
    require_once_per_residue: bool = True,
) -> list[str]:
    """
    Generic number-based atom/bead selector.

    Accepted input:
        1
        1 3 5
        1,3,5
        2-6

    When require_once_per_residue=True, each selected name must occur exactly
    once in every molecule/residue block of the selected resname. This is
    appropriate for lipid headgroup markers, chain positions, and most
    residue-specific membrane analyses.
    """
    while True:
        names = show_numbered_atom_bead_items(
            gro_file,
            resname,
            title=title,
        )

        print()
        print(instruction)

        if allow_multiple:
            print(
                "Input example: 1 3 5   or   1,3,5"
            )
        else:
            print(
                "Enter ONE number."
            )

        raw = ask(
            "Select atom/bead number"
            + (
                "(s)"
                if allow_multiple
                else ""
            )
        )

        numbers, invalid = _parse_numbered_selection(
            raw,
            len(names),
        )

        if invalid or not numbers:
            if invalid:
                print(
                    "[WARNING] Invalid atom/bead selection(s): "
                    + ", ".join(invalid)
                )
            print(
                f"[WARNING] Valid numbers are 1-{len(names)}."
            )
            continue

        if (
            not allow_multiple
            and len(numbers) != 1
        ):
            print(
                "[WARNING] Please select exactly ONE atom/bead number."
            )
            continue

        if len(numbers) < min_count:
            print(
                f"[WARNING] Please select at least {min_count} atom/bead "
                "position(s)."
            )
            continue

        selected = [
            names[number - 1]
            for number in numbers
        ]

        errors = []

        if require_once_per_residue:
            blocks = lipid_residue_blocks(
                gro_file,
                resname,
            )

            for atomname in selected:
                bad = []

                for mol_i, block in enumerate(
                    blocks,
                    start=1,
                ):
                    matches = [
                        idx
                        for idx, name in block["atoms"]
                        if name == atomname
                    ]

                    if len(matches) != 1:
                        bad.append(
                            f"molecule #{mol_i}: {len(matches)} match(es)"
                        )

                if bad:
                    errors.append(
                        f"{resname}-{atomname} is not present exactly once "
                        f"in every residue/molecule block "
                        f"({'; '.join(bad[:6])})."
                    )

        if tpr is not None:
            errors.extend(
                validate_atom_names_in_tpr(
                    tpr,
                    resname,
                    selected,
                )
            )

        if errors:
            print(
                "\n[WARNING] Selected atom/bead item(s) failed validation:"
            )

            for error in errors:
                print(
                    f"  - {error}"
                )

            print(
                "[ACTION] Please select again at THIS step."
            )
            continue

        print("\n[SELECTED]")
        for number in numbers:
            print(
                f"  [{number:>3d}] {resname}-{names[number - 1]}"
            )

        if yes_no(
            "Confirm this atom/bead selection?",
            True,
        ):
            return selected


def ask_chain_atom_names(
    gro_file: Path,
    resname: str,
    chain_label: str,
    tpr: Path | None = None,
):
    """
    Number-based sequential chain definition.

    Instead of repeatedly typing atom/bead names, the user sees entries such as:
        [1] HSPC-C1A  [2] HSPC-C2A  [3] HSPC-C3A

    The entered NUMBER ORDER is preserved and therefore defines the chain
    sequence from headgroup-side toward the tail end.
    """
    while True:
        selected = prompt_numbered_atom_bead_selection(
            gro_file,
            resname,
            title=f"ORDER PARAMETER — {resname} / {chain_label}",
            instruction=(
                f"Select at least 3 sequential {chain_label} atom/bead "
                "positions from the headgroup-side toward the tail end. "
                "IMPORTANT: enter the numbers IN CHAIN ORDER."
            ),
            allow_multiple=True,
            min_count=3,
            tpr=tpr,
            require_once_per_residue=True,
        )

        try:
            collect_equivalent_atom_groups(
                gro_file,
                resname,
                selected,
            )
        except RecoverableAnalysisError as exc:
            print(
                f"[WARNING] {exc}"
            )
            print(
                "[ACTION] Please select the chain positions again."
            )
            continue

        print(
            f"[OK] {chain_label} sequence: "
            + " > ".join(
                f"{resname}-{name}"
                for name in selected
            )
        )

        return selected




def write_order_index(
    gro_file: Path,
    resname: str,
    atom_names: list[str],
    ndx_file: Path,
    chain_label: str,
):
    groups, nmol = collect_equivalent_atom_groups(gro_file, resname, atom_names)

    # This file is intentionally created from scratch. gmx order documentation
    # requires ONLY calculation groups and explicitly warns against generic
    # groups such as System/Protein.
    with ndx_file.open("w", encoding="utf-8") as fh:
        for pos, (atom_name, indices) in enumerate(groups, start=1):
            group_name = f"{resname}_{chain_label}_P{pos:02d}_{atom_name}"
            write_ndx_group(fh, group_name, indices)

    print(f"[OK] Clean order index created: {ndx_file}")
    print(f"     Lipid molecules included : {nmol}")
    print(f"     Sequential index groups  : {len(groups)}")
    print("     Generic groups retained  : 0")
    return nmol


def run_order_command(
    xtc: Path,
    tpr: Path,
    ndx_file: Path,
    xvg_file: Path,
    begin_ps: float,
    end_ps: float | None,
    dt_ps: float,
):
    cmd = [
        GMX, "order",
        "-f", str(xtc),
        "-s", str(tpr),
        "-n", str(ndx_file),
        "-o", str(xvg_file),
        *gmx_time_args(begin_ps, end_ps, dt_ps),
        "-d", "z",
        "-szonly",
    ]

    log_file = timestamped_log_path(xvg_file.with_suffix(".log"))
    result = _animated_process(
        cmd,
        label=f"Order parameter: {ndx_file.stem}",
    )
    log_file.write_text(
        log_timestamp_header() + "COMMAND\n" + " ".join(cmd) + "\n\nSTDOUT\n" +
        (result.stdout or "") + "\nSTDERR\n" + (result.stderr or ""),
        encoding="utf-8",
    )

    if result.returncode != 0:
        print(result.stdout or "")
        print(result.stderr or "", file=sys.stderr)
        die(
            f"gmx order failed for {ndx_file.name}. Full log: {log_file}"
        )

    if not xvg_file.exists() or xvg_file.stat().st_size == 0:
        die(f"gmx order did not generate {xvg_file}")

    print(f"[OK] Z-axis order parameter generated: {xvg_file}")
    print(f"[INFO] Command log: {log_file}")
    return log_file


def analyze_order_parameters(
    gro: Path,
    tpr: Path,
    xtc: Path,
    begin_ps: float,
    end_ps: float | None,
    dt_ps: float,
):
    print("\n======================================================================")
    print(" [4/5] LIPID ORDER PARAMETER — BATCH gmx order")
    print("======================================================================")
    print(
        "The dedicated order index for each chain is created FROM SCRATCH and "
        "contains only equivalent sequential chain-atom groups."
    )
    print(
        "Calculation command: gmx order -f XTC -s TPR -n ORDER.ndx "
        "-o ORDER.xvg -d z -szonly"
    )
    print()
    print("*** GROMACS MANUAL LIMITATION ***")
    print(
        "Current GROMACS documentation states that gmx order only works for "
        "saturated carbons and united-atom force fields. For other models "
        "(including typical coarse-grained Martini or many all-atom models), "
        "the numeric output should not be treated as a validated lipid order "
        "parameter without an appropriate alternative method."
    )

    choice = ask("Run the order-parameter module? (Y/n)", "Y")
    if choice.lower().startswith("n"):
        print("[INFO] Order-parameter module skipped by user.")
        return []

    ORDER_DIR.mkdir(parents=True, exist_ok=True)
    results = []

    while True:
        print("\n------------------------- NEW LIPID FOR ORDER --------------------------")
        resname = ask("Lipid residue name (resname) to analyze")
        if not resname or " " in resname:
            print("[WARNING] Invalid residue name; try again.")
            continue

        blocks = lipid_residue_blocks(gro, resname)
        if not blocks:
            print(f"[WARNING] Residue '{resname}' was not found in {gro}.")
            again = ask("Try another residue name? (Y/n)", "Y")
            if again.lower().startswith("n"):
                break
            continue

        print(f"[OK] Found {len(blocks)} molecules/residue blocks named {resname}.")

        while True:
            chain_raw = ask("Does this lipid have one chain or two chains? Enter 1 or 2")
            if chain_raw in {"1", "2"}:
                nchain = int(chain_raw)
                break
            print("[WARNING] Please enter 1 or 2.")

        if nchain == 1:
            atom_names = ask_chain_atom_names(gro, resname, "single-chain")
            base = f"{resname}-order"
            ndx = ORDER_DIR / f"{base}.ndx"
            xvg = ORDER_DIR / f"{base}.xvg"
            nmol = write_order_index(
                gro, resname, atom_names, ndx, "chain"
            )
            log = run_order_command(
                xtc, tpr, ndx, xvg, begin_ps, end_ps, dt_ps
            )
            results.append({
                "lipid": resname,
                "chain": "single",
                "molecules": nmol,
                "atoms": atom_names,
                "ndx": ndx,
                "xvg": xvg,
                "log": log,
            })

        else:
            for chain_no in (1, 2):
                print(
                    f"\n================ {resname}: BUILD CHAIN {chain_no} ================"
                )
                atom_names = ask_chain_atom_names(
                    gro, resname, f"chain{chain_no}"
                )
                base = f"{resname}-order-chain{chain_no}"
                ndx = ORDER_DIR / f"{base}.ndx"
                xvg = ORDER_DIR / f"{base}.xvg"
                nmol = write_order_index(
                    gro, resname, atom_names, ndx, f"chain{chain_no}"
                )
                log = run_order_command(
                    xtc, tpr, ndx, xvg, begin_ps, end_ps, dt_ps
                )
                results.append({
                    "lipid": resname,
                    "chain": f"chain{chain_no}",
                    "molecules": nmol,
                    "atoms": atom_names,
                    "ndx": ndx,
                    "xvg": xvg,
                    "log": log,
                })

            print(f"[OK] Both chains of {resname} have been completed.")

        more = ask("Continue calculating another lipid molecule? (y/N)", "N")
        if not more.lower().startswith("y"):
            break

    if results:
        summary_file = ORDER_DIR / "order_parameter_summary.xlsx"
        headers = [
            "lipid",
            "chain",
            "molecules",
            "sequential_atoms",
            "index_file",
            "order_file",
        ]
        rows = [
            [
                r["lipid"],
                r["chain"],
                r["molecules"],
                " > ".join(r["atoms"]),
                str(r["ndx"]),
                str(r["xvg"]),
            ]
            for r in results
        ]
        write_excel_table(
            summary_file,
            "Order Parameter",
            headers,
            rows,
        )
        print(f"[OK] Batch order Excel summary: {summary_file}")

    return results


def make_plots(apl_csv: Path, thick_csv: Path, sasa_csv: Path):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "[INFO] matplotlib is not installed; PNG plotting skipped.\n"
            "       You can use xmgrace on the XVG files."
        )
        return

    jobs = [
        (
            apl_csv,
            APL_DIR / "apl_vs_time.png",
            "Time_ns",
            "APL_nm2_per_lipid",
            "Area per lipid vs Time",
            "Time (ns)",
            r"Area per lipid (nm$^2$/lipid)",
        ),
        (
            thick_csv,
            THICK_DIR / "membrane_thickness_vs_time.png",
            "Time_ns",
            "Membrane_thickness_nm",
            "Membrane thickness vs Time",
            "Time (ns)",
            "Membrane thickness (nm)",
        ),
        (
            sasa_csv,
            SASA_DIR / "membrane_sasa_vs_time.png",
            "Time_ns",
            "Membrane_SASA_nm2",
            "Membrane SASA vs Time",
            "Time (ns)",
            r"SASA (nm$^2$)",
        ),
    ]

    for (
        csv_file,
        png_file,
        xkey,
        ykey,
        title,
        xlabel,
        ylabel,
    ) in jobs:
        x, y = [], []

        with csv_file.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                x.append(float(row[xkey]))
                y.append(float(row[ykey]))

        fig, ax = plt.subplots(figsize=(7.2, 5.0))
        ax.plot(x, y, linewidth=1.2)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        fig.tight_layout()
        fig.savefig(png_file, dpi=300)
        plt.close(fig)

    print("[OK] Three 300-dpi PNG plots generated.")



# ======================================================================
# MENU-DRIVEN INTERACTIVE LAYER (v5)
# ======================================================================

def _files_with_suffix(suffix: str) -> list[Path]:
    """Detect candidate files in the current working directory."""
    suffix = suffix.lower().lstrip(".")
    cwd = Path.cwd()
    files = [
        p for p in cwd.iterdir()
        if p.is_file() and p.suffix.lower() == f".{suffix}"
    ]

    preferred_words = (
        "production", "prod", "md", "step7", "npt", "nvt", "trjconv",
    )

    def sort_key(p: Path):
        low = p.name.lower()
        priority = min(
            (i for i, word in enumerate(preferred_words) if word in low),
            default=len(preferred_words),
        )
        try:
            newer_first = -p.stat().st_mtime
        except OSError:
            newer_first = 0
        return (priority, newer_first, low)

    return sorted(files, key=sort_key)


def _choose_detected_file(
    candidates: list[Path],
    label: str,
    previous: Path | None = None,
) -> Path:
    """Select among several detected candidates using a number."""
    if len(candidates) == 1:
        return candidates[0]

    default_number = 1
    print(f"\nMultiple {label} files were detected:")
    for i, path in enumerate(candidates, start=1):
        mark = ""
        if previous is not None:
            try:
                if path.resolve() == previous.resolve():
                    default_number = i
                    mark = "  <- previously used"
            except OSError:
                pass
        print(f"  [{i}] {path.name}{mark}")

    while True:
        raw = ask(f"Choose {label} file number", str(default_number))
        try:
            number = int(raw)
            if 1 <= number <= len(candidates):
                return candidates[number - 1]
        except ValueError:
            pass
        print(f"[WARNING] Enter a number from 1 to {len(candidates)}.")



def validate_selected_input_file(
    path: Path,
    extension: str,
) -> tuple[bool, str]:
    """
    Validate an input file at the file-selection step.

    GRO is structurally parsed. Binary files are checked for existence and
    non-zero size here; their detailed content is validated by GROMACS when
    the corresponding analysis runs.
    """
    if not path.exists() or not path.is_file():
        return False, f"File does not exist: {path}"

    if path.stat().st_size <= 0:
        return False, f"File is empty: {path}"

    extension = extension.lower().lstrip(".")

    if extension == "gro":
        try:
            lines = path.read_text(
                encoding="utf-8",
                errors="replace",
            ).splitlines()
            if len(lines) < 3:
                return False, "GRO file has too few lines."

            natoms = int(lines[1].strip())
            if natoms <= 0:
                return False, "GRO file reports zero atoms."

            if len(lines) < natoms + 3:
                return False, (
                    "GRO atom count does not match the file length."
                )
        except Exception as exc:
            return False, f"GRO file cannot be parsed: {exc}"

    return True, "OK"



def resolve_analysis_files(
    requirements: list[tuple[str, str, str]],
    state: dict,
) -> dict[str, Path]:
    """
    Detect the file types required by the selected calculation in the CURRENT
    working directory.

    The user confirms the detected set once. If several files of one type are
    present, the user chooses by number instead of typing full file names.
    Missing types alone fall back to manual path entry.
    """
    print("\n---------------------- AUTOMATIC FILE DETECTION ------------------------")
    print(f"Current working directory: {Path.cwd()}")

    detected = {}
    all_present = True

    for key, extension, label in requirements:
        candidates = _files_with_suffix(extension)
        detected[key] = candidates

        if not candidates:
            all_present = False
            print(f"  {label:<18}: NOT FOUND (*.{extension})")
        elif len(candidates) == 1:
            print(f"  {label:<18}: {candidates[0].name}")
        else:
            print(f"  {label:<18}: {len(candidates)} candidates")
            for i, path in enumerate(candidates, start=1):
                print(f"      [{i}] {path.name}")

    if all_present:
        use_detected = yes_no(
            "Use detected file(s) in the current folder for this analysis?",
            True,
        )
    else:
        print(
            "[INFO] Some required file types are missing. "
            "Detected files can still be reused."
        )
        use_detected = yes_no(
            "Use detected files where available and enter only missing ones?",
            True,
        )

    chosen = {}

    for key, extension, label in requirements:
        previous = state.get(key)
        previous_path = Path(previous) if previous is not None else None

        while True:
            if use_detected and detected[key]:
                candidate = _choose_detected_file(
                    detected[key],
                    label,
                    previous_path,
                )
            else:
                default = (
                    str(previous)
                    if previous is not None
                    else f"md.{extension}"
                )
                candidate = ask_existing_path(
                    label,
                    default,
                )

            ok_file, file_detail = validate_selected_input_file(
                candidate,
                extension,
            )

            if ok_file:
                chosen[key] = candidate
                break

            print(
                f"[WARNING] The selected {label} is not usable: "
                f"{file_detail}"
            )
            print(
                f"[ACTION] Please select/enter the {label} again "
                "at THIS step."
            )

            # If an auto-detected candidate is unusable, keep the user at
            # THIS file-input step but allow manual replacement immediately.
            if detected[key]:
                if yes_no(
                    f"Choose another detected {label} file?",
                    len(detected[key]) > 1,
                ):
                    continue

                use_detected = False
                print(
                    f"[ACTION] Please enter the correct {label} path "
                    "at THIS step."
                )
                continue

    state.update(chosen)

    print("\n[INPUT FILES CONFIRMED]")
    for key, extension, label in requirements:
        print(f"  {label:<18}: {chosen[key]}")

    return chosen


def ask_existing_path(prompt: str, default: str | None = None) -> Path:
    """Prompt until an existing file is supplied."""
    while True:
        raw = ask(prompt, default) if default is not None else ask(prompt)
        p = Path(raw)
        if p.exists() and p.is_file():
            return p
        print(f"[WARNING] File not found: {p}")


def yes_no(prompt: str, default: bool = True) -> bool:
    default_text = "Y" if default else "N"
    while True:
        ans = ask(prompt, default_text).strip().lower()
        if ans in {"y", "yes"}:
            return True
        if ans in {"n", "no"}:
            return False
        print("[WARNING] Please enter Y or N.")



def show_completion_banner(parameter_name: str):
    """
    Show a success banner after one parameter calculation finishes normally.
    """
    print("\n")
    print("=" * 94)
    print(f"[COMPLETED] {parameter_name} calculation completed successfully.")
    print("=" * 94)
    print()
    print("                         CONGRATULATIONS")
    print()
    print(r"   ██████╗ ██████╗ ███╗   ██╗ ██████╗ ██████╗  █████╗ ████████╗███████╗")
    print(r"  ██╔════╝██╔═══██╗████╗  ██║██╔════╝ ██╔══██╗██╔══██╗╚══██╔══╝██╔════╝")
    print(r"  ██║     ██║   ██║██╔██╗ ██║██║  ███╗██████╔╝███████║   ██║   ███████╗")
    print(r"  ██║     ██║   ██║██║╚██╗██║██║   ██║██╔══██╗██╔══██║   ██║   ╚════██║")
    print(r"  ╚██████╗╚██████╔╝██║ ╚████║╚██████╔╝██║  ██║██║  ██║   ██║   ███████║")
    print(r"   ╚═════╝ ╚═════╝ ╚═╝  ╚═══╝ ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝   ╚═╝   ╚══════╝")
    print()
    print("      You have completed one parameter calculation and are one step closer to scientific research.")
    print("      If you have any problems using this tool, please do not hesitate to contact minfangfeng@bucm.edu.cn.")
    print("=" * 94)



def post_calculation_menu(parameter_name: str, continue_text: str) -> str:
    """
    Unified action menu after one parameter calculation.

    Returns:
      "continue" -> calculate current parameter again
      "menu"     -> return to main parameter menu
      "exit"     -> terminate program
    """
    show_completion_banner(parameter_name)

    while True:
        print("\n======================================================================")
        print(f" {parameter_name} CALCULATION COMPLETED")
        print("======================================================================")
        print(f" [1] {continue_text}")
        print(" [2] Return to MAIN MENU and calculate another parameter")
        print(" [3] Exit program")
        print("======================================================================")
        choice = ask("Choose next action", "2")
        if choice == "1":
            return "continue"
        if choice == "2":
            return "menu"
        if choice == "3":
            return "exit"
        print("[WARNING] Please choose 1, 2, or 3.")


def prompt_time_window(require_dt: bool = False):
    print("\n--------------------------- ANALYSIS WINDOW ----------------------------")
    print("Time unit: ps (1000 ps = 1 ns).")
    begin_ps = ask_float("Analysis start time", 0.0)
    end_ps = ask_optional_float("Analysis end time")

    if require_dt:
        dt_ps = ask_float(
            "Sampling interval; 0 means every stored trajectory frame",
            0.0,
        )
    else:
        dt_ps = 0.0

    return begin_ps, end_ps, dt_ps


def residue_counts_from_gro(gro: Path) -> Counter:
    counts = Counter()
    for block in read_gro_residue_blocks(gro):
        counts[block["resname"]] += 1
    return counts


def show_residue_names(gro: Path):
    counts = residue_counts_from_gro(gro)
    print(f"\n[INFO] Residue names detected in {gro}:")
    items = sorted(counts.items(), key=lambda x: x[0])
    for i in range(0, len(items), 4):
        chunk = items[i:i + 4]
        print(
            "       " +
            "   ".join(f"{name}({n})" for name, n in chunk)
        )




def _count_ndx_atoms(ndx_file: Path) -> int:
    """Count atom indices in a temporary one-group NDX file."""
    if not ndx_file.exists():
        return 0

    count = 0
    with ndx_file.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("["):
                continue
            for token in s.split():
                try:
                    int(token)
                    count += 1
                except ValueError:
                    pass
    return count


def tpr_selection_count(tpr: Path | None, selection: str) -> int | None:
    """
    Evaluate a static GROMACS selection against TPR.

    Returns selected atom count. Returns None if this validation mechanism
    itself is unavailable, in which case GRO validation is still used.
    """
    if tpr is None:
        return None

    OUTDIR.mkdir(parents=True, exist_ok=True)
    tmp = OUTDIR / f".tpr_check_{abs(hash((str(tpr), selection))) % 100000000}.ndx"

    try:
        if tmp.exists():
            tmp.unlink()

        result = subprocess.run(
            [
                GMX, "select",
                "-s", str(tpr),
                "-on", str(tmp),
                "-select", selection,
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        if result.returncode != 0:
            return None

        return _count_ndx_atoms(tmp)

    except Exception:
        return None

    finally:
        try:
            if tmp.exists():
                tmp.unlink()
        except OSError:
            pass


def validate_resname_in_tpr(
    tpr: Path | None,
    resname: str,
) -> tuple[bool, str]:
    """Check whether a residue/molecule type exists in the selected TPR."""
    if tpr is None:
        return True, "TPR validation skipped."

    count = tpr_selection_count(
        tpr,
        f'resname "{resname}"',
    )

    if count is None:
        return True, "TPR validation unavailable; GRO validation retained."

    if count == 0:
        return False, (
            f"resname '{resname}' exists in GRO but was not found in "
            f"TPR '{tpr}'. Please check that GRO and TPR belong to the same system."
        )

    return True, f"TPR matched {count} atoms."


def validate_atom_in_tpr(
    tpr: Path | None,
    resname: str,
    atomname: str,
) -> tuple[bool, str]:
    """Check one atom/bead name under one resname in TPR."""
    if tpr is None:
        return True, "TPR validation skipped."

    count = tpr_selection_count(
        tpr,
        f'resname "{resname}" and name "{atomname}"',
    )

    if count is None:
        return True, "TPR validation unavailable; GRO validation retained."

    if count == 0:
        return False, (
            f"atom/bead '{atomname}' under resname '{resname}' "
            f"was not found in TPR '{tpr}'."
        )

    return True, f"TPR matched {count} atoms."


def validate_atom_names_in_tpr(
    tpr: Path | None,
    resname: str,
    atom_names: list[str],
) -> list[str]:
    """Return all TPR validation errors for a set of atom/bead names."""
    errors = []
    for atomname in atom_names:
        ok, detail = validate_atom_in_tpr(
            tpr,
            resname,
            atomname,
        )
        if not ok:
            errors.append(detail)
    return errors



def numbered_residue_items(gro: Path) -> list[tuple[str, int]]:
    """
    Return all residue/molecule types in first-appearance order together
    with their number of residue blocks in the GRO file.

    No automatic lipid/non-lipid classification is attempted because custom
    GROMACS/Martini residue names can be arbitrary.
    """
    blocks = read_gro_residue_blocks(gro)
    counts = Counter(block["resname"] for block in blocks)

    ordered = []
    seen = set()
    for block in blocks:
        name = block["resname"]
        if name not in seen:
            seen.add(name)
            ordered.append(name)

    return [(name, counts[name]) for name in ordered]


def show_numbered_residue_items(
    gro: Path,
    title: str = "RESIDUE / MOLECULE TYPE SELECTION",
) -> list[tuple[str, int]]:
    """
    Print every residue/molecule type as:
        [1] HSPC  286
        [2] CHOL  188
        ...
    """
    items = numbered_residue_items(gro)
    if not items:
        die(f"No residue/molecule types could be read from {gro}.")

    print("\n" + "-" * 78)
    print(f" {title}")
    print("-" * 78)
    print(f"Source GRO: {gro}")
    print()
    print(f"  {'No.':>4s}  {'resname':<18s} {'count':>10s}")
    print("  " + "-" * 40)

    for number, (resname, count) in enumerate(items, start=1):
        print(f"  [{number:>2d}]  {resname:<18s} {count:>10d}")

    return items


def _parse_numbered_selection(
    raw: str,
    nitems: int,
) -> tuple[list[int], list[str]]:
    """
    Parse 1-based numeric selections separated by whitespace, commas,
    semicolons, or tabs. Preserve order and remove duplicates.
    """
    tokens = [x for x in re.split(r"[\s,;]+", raw.strip()) if x]
    valid = []
    invalid = []

    for token in tokens:
        try:
            number = int(token)
        except ValueError:
            invalid.append(token)
            continue

        if not (1 <= number <= nitems):
            invalid.append(token)
            continue

        if number not in valid:
            valid.append(number)

    return valid, invalid


def prompt_numbered_residue_selection(
    gro: Path,
    *,
    title: str,
    instruction: str,
    allow_multiple: bool,
    selection_name: str,
    tpr: Path | None = None,
) -> list[str]:
    """
    Generic numbered residue/molecule selector used by Thickness, SASA,
    Order Parameter and MSD.

    Multiple-selection input supports:
      - one line: 1 2 3
      - one line: 1,2,3
      - line-by-line:
            1 <Enter>
            2 <Enter>
            3 <Enter>
            <blank Enter to finish>
    """
    items = show_numbered_residue_items(gro, title=title)
    nitems = len(items)

    print()
    print(instruction)

    if allow_multiple:
        print()
        print("Input examples:")
        print("  Same line : 1 2 3")
        print("  Same line : 1,2,3")
        print("  New lines : enter 1, then 2, then 3; blank Enter finishes")
    else:
        print("Enter ONE number from the table above.")

    selected_numbers = []

    while True:
        raw = input(
            f"Select {selection_name} number"
            f"{'(s)' if allow_multiple else ''}: "
        ).strip()

        if not raw:
            print("[WARNING] Please make a selection.")
            continue

        numbers, invalid = _parse_numbered_selection(raw, nitems)

        if invalid:
            print(
                "[WARNING] Invalid selection(s): " + ", ".join(invalid)
            )
            print(f"          Valid numbers are 1-{nitems}.")
            continue

        if not numbers:
            print("[WARNING] No valid number was selected.")
            continue

        if not allow_multiple and len(numbers) != 1:
            print("[WARNING] Please select exactly ONE number.")
            continue

        for number in numbers:
            if number not in selected_numbers:
                selected_numbers.append(number)

        if not allow_multiple:
            break

        token_count = len(
            [x for x in re.split(r"[\s,;]+", raw) if x]
        )

        # If several numbers were provided in one line, treat that input
        # as the completed batch selection.
        if token_count > 1:
            break

        # If only one number was entered, allow line-by-line continuation.
        while True:
            raw_more = input(
                f"Add another {selection_name} number "
                "(or press Enter to finish): "
            ).strip()

            if not raw_more:
                break

            more_numbers, invalid = _parse_numbered_selection(
                raw_more,
                nitems,
            )

            if invalid:
                print(
                    "[WARNING] Invalid selection(s): "
                    + ", ".join(invalid)
                )
                print(f"          Valid numbers are 1-{nitems}.")
                continue

            for number in more_numbers:
                if number not in selected_numbers:
                    selected_numbers.append(number)

        break

    selected = [items[i - 1][0] for i in selected_numbers]

    if tpr is not None:
        tpr_errors = []
        for resname in selected:
            ok, detail = validate_resname_in_tpr(
                tpr,
                resname,
            )
            if not ok:
                tpr_errors.append(detail)

        if tpr_errors:
            print()
            print(
                "[WARNING] The selected residue/molecule type is not "
                "present in the current TPR:"
            )
            for error in tpr_errors:
                print(f"  - {error}")
            print(
                "[ACTION] Please make the selection again at THIS step."
            )
            return prompt_numbered_residue_selection(
                gro,
                title=title,
                instruction=instruction,
                allow_multiple=allow_multiple,
                selection_name=selection_name,
                tpr=tpr,
            )

    print(f"\n[SELECTED {selection_name.upper()}]")
    for number in selected_numbers:
        resname, count = items[number - 1]
        print(
            f"  [{number:>2d}] {resname:<18s} "
            f"count = {count}"
        )

    if not yes_no(
        f"Confirm the selected {selection_name}"
        f"{'s' if len(selected) != 1 else ''}?",
        True,
    ):
        return prompt_numbered_residue_selection(
            gro,
            title=title,
            instruction=instruction,
            allow_multiple=allow_multiple,
            selection_name=selection_name,
            tpr=tpr,
        )

    return selected



def prompt_membrane_lipid_names(
    gro: Path,
    tpr: Path | None = None,
) -> list[str]:
    """
    Select all lipid residue/molecule types that constitute the membrane
    using numbered choices rather than typing resnames manually.

    Used by:
      - membrane SASA
      - whole-membrane regional MSD
    """
    return prompt_numbered_residue_selection(
        gro,
        title="MEMBRANE LIPID SELECTION",
        instruction=(
            "Select ALL residue/molecule types that are membrane lipids. "
            "Do not select water, ions, protein residues, free ligand, "
            "or other non-membrane species."
        ),
        allow_multiple=True,
        selection_name="membrane lipid",
        tpr=tpr,
    )




def prompt_thickness_lipid_definitions(
    gro: Path,
    tpr: Path | None = None,
) -> list[tuple[str, str]]:
    """
    Membrane-thickness workflow with fully numbered atom/bead choices:
      1) select all membrane lipid species by NUMBER;
      2) for each lipid, select ONE representative headgroup atom/bead by
         NUMBER, e.g. [1] CHOL-ROH [2] CHOL-C1 ...

    The selected marker must occur exactly once in every molecule of the lipid.
    """
    lipid_names = prompt_numbered_residue_selection(
        gro,
        title="MEMBRANE THICKNESS — LIPID SELECTION",
        instruction=(
            "Select ALL lipid residue/molecule types that form the membrane. "
            "The next step will show the available atoms/beads for each lipid "
            "as numbered choices."
        ),
        allow_multiple=True,
        selection_name="membrane lipid",
        tpr=tpr,
    )

    lipid_defs = []

    print(
        "\n---------------- THICKNESS HEADGROUP REFERENCE SETUP ----------------"
    )
    print(
        "For each selected lipid, choose exactly ONE representative "
        "headgroup atom/bead by NUMBER."
    )

    for i, resname in enumerate(
        lipid_names,
        start=1,
    ):
        print(
            f"\nLipid {i}/{len(lipid_names)}: {resname}"
        )

        head = prompt_numbered_atom_bead_selection(
            gro,
            resname,
            title=f"THICKNESS HEADGROUP MARKER — {resname}",
            instruction=(
                "Choose ONE representative headgroup atom/bead. "
                "For example, CHOL usually uses ROH; phospholipids often "
                "use a phosphorus/phosphate marker when appropriate."
            ),
            allow_multiple=False,
            min_count=1,
            tpr=tpr,
            require_once_per_residue=True,
        )[0]

        lipid_defs.append(
            (
                resname,
                head,
            )
        )

    return lipid_defs




def build_membrane_only_index(gro_file: Path, lipid_names: list[str]):
    """
    Generate the normal GROMACS default index and append a clean
    user-defined [ membrane ] group containing ALL atoms of all specified
    lipid residues.
    """
    print("\n---------------------------- INDEX SETUP -------------------------------")
    print(f"[INFO] Generating default GROMACS index from {gro_file} ...")

    run(
        [GMX, "make_ndx", "-f", str(gro_file), "-o", str(BASE_NDX)],
        input_text="q\n",
        quiet=True,
    )

    if not BASE_NDX.exists() or BASE_NDX.stat().st_size == 0:
        die("Default GROMACS index was not generated.")

    membrane_indices = get_membrane_atom_indices(
        gro_file,
        set(lipid_names),
    )
    if not membrane_indices:
        die("No membrane atoms were selected.")

    with CUSTOM_NDX.open("w", encoding="utf-8") as fh:
        write_ndx_group(fh, "membrane", membrane_indices)

    if INDEX.exists():
        backup = INDEX.with_name(f"{INDEX.name}.bak")
        shutil.copy2(INDEX, backup)
        print(f"[INFO] Existing {INDEX} backed up as {backup}")

    base_bytes = BASE_NDX.read_bytes()
    with INDEX.open("wb") as out:
        out.write(base_bytes)
        if not base_bytes.endswith(b"\n"):
            out.write(b"\n")
        out.write(CUSTOM_NDX.read_bytes())

    shutil.copy2(INDEX, OUTDIR / "index.ndx")

    print("[OK] Final index generated: index.ndx")
    print(
        "[OK] [ membrane ] contains all atoms from: " +
        ", ".join(lipid_names)
    )


def make_single_plot(
    csv_file: Path,
    png_file: Path,
    xkey: str,
    ykey: str,
    title: str,
    xlabel: str,
    ylabel: str,
):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "[INFO] matplotlib is not installed; PNG plotting skipped. "
            "The XVG file can be opened with xmgrace."
        )
        return

    x, y = [], []
    with csv_file.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            x.append(float(row[xkey]))
            y.append(float(row[ykey]))

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(x, y, linewidth=1.2)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(png_file, dpi=300)
    plt.close(fig)
    print(f"[OK] Plot generated: {png_file}")



def apl_numbered_residue_table(gro: Path) -> list[tuple[str, int]]:
    """
    Return every unique GRO residue name and its residue/molecule count,
    preserving the order in which residue types first appear in the GRO.

    The list is deliberately NOT auto-classified as lipid/non-lipid because
    custom force fields and Martini systems may use arbitrary residue names.
    """
    blocks = read_gro_residue_blocks(gro)
    counts = Counter(block["resname"] for block in blocks)

    ordered_names = []
    seen = set()
    for block in blocks:
        name = block["resname"]
        if name not in seen:
            seen.add(name)
            ordered_names.append(name)

    return [(name, counts[name]) for name in ordered_names]


def show_apl_numbered_residue_table(gro: Path) -> list[tuple[str, int]]:
    """
    Print all residue/molecule types with 1, 2, 3... labels specifically
    for APL membrane-lipid selection.
    """
    items = apl_numbered_residue_table(gro)

    if not items:
        die(f"No residues/molecules could be read from {gro}.")

    print("\n-------------------- APL MEMBRANE LIPID SELECTION ---------------------")
    print(f"All residue/molecule types detected in: {gro}")
    print()
    print(f"  {'No.':>4s}  {'resname':<16s} {'count':>10s}")
    print("  " + "-" * 36)

    for number, (resname, count) in enumerate(items, start=1):
        print(f"  [{number:>2d}]  {resname:<16s} {count:>10d}")

    print()
    print(
        "Select ALL lipid species that actually constitute the membrane."
    )
    print(
        "Do NOT select water, ions, protein residues, free ligand, or other "
        "non-membrane molecules."
    )
    print(
        "IMPORTANT: this produces one GLOBAL average APL for the selected "
        "membrane lipids. Selecting only one species from a mixed membrane "
        "does NOT produce a species-specific Voronoi APL."
    )

    return items


def _parse_apl_selection_numbers(
    raw: str,
    nitems: int,
) -> tuple[list[int], list[str]]:
    """
    Parse numeric choices separated by spaces, commas, semicolons or tabs.
    Returns valid unique 1-based numbers and invalid tokens.
    """
    tokens = [x for x in re.split(r"[\s,;]+", raw.strip()) if x]
    valid = []
    invalid = []

    for token in tokens:
        try:
            number = int(token)
        except ValueError:
            invalid.append(token)
            continue

        if not (1 <= number <= nitems):
            invalid.append(token)
            continue

        if number not in valid:
            valid.append(number)

    return valid, invalid


def prompt_apl_lipid_selection(
    gro: Path,
) -> tuple[list[str], dict[str, int], int]:
    """
    Let the user select membrane lipid species by NUMBER.

    Supported input styles:
      - one-line: 1 2 3
      - one-line: 1,2,3
      - line-by-line:
            1 <Enter>
            2 <Enter>
            3 <Enter>
            <blank Enter to finish>

    Returns:
      selected resnames
      selected count dictionary
      total selected lipid molecules
    """
    items = show_apl_numbered_residue_table(gro)
    nitems = len(items)

    print()
    print("Selection input examples:")
    print("  Same line : 1 2 3")
    print("  Same line : 1,2,3")
    print("  New lines : enter 1, then 2, then 3; press blank Enter to finish")
    print()

    selected_numbers = []

    # First entry. If the user enters multiple numbers at once, treat that
    # as a complete one-line selection. A single number enters line-by-line
    # mode so more can be added.
    while True:
        raw = input("Select membrane lipid number(s): ").strip()

        if not raw:
            print("[WARNING] Please select at least one membrane lipid.")
            continue

        numbers, invalid = _parse_apl_selection_numbers(raw, nitems)

        if invalid:
            print(
                "[WARNING] Invalid selection(s): " + ", ".join(invalid)
            )
            print(f"          Valid numbers are 1-{nitems}.")
            continue

        if not numbers:
            print("[WARNING] No valid selection was entered.")
            continue

        for number in numbers:
            if number not in selected_numbers:
                selected_numbers.append(number)

        # More than one token means the user used the space/comma form,
        # so no extra blank line is required.
        token_count = len(
            [x for x in re.split(r"[\s,;]+", raw) if x]
        )
        if token_count > 1:
            break

        # Single-number form: allow additional selections on new lines.
        while True:
            raw_more = input(
                "Add another lipid number (or press Enter to finish): "
            ).strip()

            if not raw_more:
                break

            more_numbers, invalid = _parse_apl_selection_numbers(
                raw_more,
                nitems,
            )

            if invalid:
                print(
                    "[WARNING] Invalid selection(s): " + ", ".join(invalid)
                )
                print(f"          Valid numbers are 1-{nitems}.")
                continue

            for number in more_numbers:
                if number not in selected_numbers:
                    selected_numbers.append(number)

        break

    selected_resnames = [items[i - 1][0] for i in selected_numbers]
    selected_counts = {
        items[i - 1][0]: items[i - 1][1]
        for i in selected_numbers
    }
    total_selected = sum(selected_counts.values())

    print("\n[APL SELECTED MEMBRANE LIPIDS]")
    for number in selected_numbers:
        resname, count = items[number - 1]
        print(f"  [{number:>2d}] {resname:<16s} count = {count}")

    print(f"\n  Selected lipid species : {len(selected_resnames)}")
    print(f"  Total selected lipids  : {total_selected}")

    if not yes_no("Confirm this membrane-lipid selection?", True):
        return prompt_apl_lipid_selection(gro)

    return selected_resnames, selected_counts, total_selected


def write_apl_selection_summary(
    gro: Path,
    selected_resnames: list[str],
    selected_counts: dict[str, int],
    total_selected: int,
    nleaflet: float,
):
    """Record exactly which molecules were used as the APL denominator."""
    APL_DIR.mkdir(parents=True, exist_ok=True)
    out = APL_DIR / "apl_membrane_lipid_selection.txt"

    lines = [
        "======================================================================",
        "APL MEMBRANE LIPID SELECTION",
        "======================================================================",
        f"GRO file                 : {gro}",
        "",
        "Selected membrane lipid species:",
    ]

    for resname in selected_resnames:
        lines.append(
            f"  {resname:<16s} : {selected_counts[resname]} molecules"
        )

    lines += [
        "",
        f"Total selected lipids    : {total_selected}",
        f"N_leaflet used for APL   : {nleaflet:g}",
        "",
        "APL definition:",
        "  APL(t) = Box-X(t) * Box-Y(t) / N_leaflet",
        "",
        "Note:",
        "  The selected species are treated together as the membrane lipid set.",
        "  This is a GLOBAL projected average APL, not a lipid-type-specific",
        "  Voronoi area per molecule.",
        "======================================================================",
        "",
    ]

    out.write_text("\n".join(lines), encoding="utf-8")
    return out




def prompt_apl_nleaflet(
    gro: Path,
) -> tuple[float, list[str], dict[str, int], int]:
    """
    APL denominator is built only from USER-SELECTED membrane lipid species.

    Workflow:
      1) enumerate every resname in GRO with a number;
      2) user selects all membrane lipid species by number;
      3) count only those selected lipid residues/molecules;
      4) propose N_leaflet = N_selected_total / 2 for a symmetric bilayer;
      5) allow manual override for asymmetric membranes.
    """
    (
        selected_resnames,
        selected_counts,
        total_selected,
    ) = prompt_apl_lipid_selection(gro)

    print("\n--------------------------- APL LIPID NUMBER ---------------------------")
    print(
        "Only the membrane lipid species selected above are included in the "
        "APL denominator."
    )
    print(f"Selected total membrane lipids = {total_selected}")

    proposed = total_selected / 2.0

    if total_selected % 2 == 0:
        print(
            f"Assuming a symmetric bilayer: "
            f"N_leaflet = {total_selected}/2 = {proposed:g}"
        )
    else:
        print(
            "[WARNING] The selected total lipid count is odd, so the two "
            "leaflets cannot contain exactly equal lipid numbers."
        )
        print(
            f"N_total/2 = {proposed:g} can only be an average denominator."
        )

    if yes_no(
        f"Use N_leaflet = {proposed:g} for GLOBAL APL?",
        True,
    ):
        nleaflet = proposed
    else:
        nleaflet = ask_float(
            "Enter the actual N_leaflet divisor for APL"
        )
        if nleaflet <= 0:
            die("N_leaflet must be > 0.")

    return (
        nleaflet,
        selected_resnames,
        selected_counts,
        total_selected,
    )




def run_apl_menu(state: dict) -> bool:
    """
    APL workflow:
      - automatically detect EDR and GRO from the current directory;
      - enumerate all GRO residue/molecule types with numeric labels;
      - user selects all membrane lipids by number;
      - count ONLY selected lipid molecules;
      - calculate GLOBAL APL with selected N_leaflet.
    """
    while True:
        print("\n======================================================================")
        print(" APL CALCULATION")
        print("======================================================================")
        print(
            "APL will use only the lipid species that YOU identify as membrane "
            "components from the numbered GRO residue list."
        )

        # APL requires EDR for Box-X/Box-Y and GRO for membrane-lipid counting.
        files = resolve_analysis_files(
            [
                ("edr", "edr", "EDR file"),
                ("gro", "gro", "GRO file"),
            ],
            state,
        )
        edr = files["edr"]
        gro = files["gro"]

        (
            nleaflet,
            selected_resnames,
            selected_counts,
            total_selected,
        ) = prompt_apl_nleaflet(gro)

        selection_file = write_apl_selection_summary(
            gro,
            selected_resnames,
            selected_counts,
            total_selected,
            nleaflet,
        )

        print("\n[APL DENOMINATOR CONFIRMED]")
        print(
            "  Selected lipids : "
            + ", ".join(
                f"{name}({selected_counts[name]})"
                for name in selected_resnames
            )
        )
        print(f"  N_total selected: {total_selected}")
        print(f"  N_leaflet used  : {nleaflet:g}")
        print(f"  Selection record: {selection_file}")

        begin_ps, end_ps, _ = prompt_time_window(require_dt=False)

        apl_xvg, apl_csv, apl_stats = analyze_apl(
            edr,
            nleaflet,
            begin_ps,
            end_ps,
        )

        # Add the membrane composition to the XVG header without altering
        # the numerical APL calculation.
        try:
            raw = apl_xvg.read_text(encoding="utf-8")
            composition = ", ".join(
                f"{name}:{selected_counts[name]}"
                for name in selected_resnames
            )
            prefix = (
                f"# Selected membrane lipid resnames: "
                f"{', '.join(selected_resnames)}\n"
                f"# Selected lipid counts: {composition}\n"
                f"# Total selected membrane lipids: {total_selected}\n"
            )
            apl_xvg.write_text(prefix + raw, encoding="utf-8")
        except Exception as exc:
            print(
                f"[WARNING] Could not add APL selection metadata to XVG: {exc}"
            )

        # Append selection metadata to the statistics report.
        try:
            old_stats = apl_stats.read_text(encoding="utf-8")
            selection_text = (
                "\nAPL membrane selection\n"
                "------------------------------------------------------------\n"
                f"Selected resnames : {', '.join(selected_resnames)}\n"
                f"Selected counts   : "
                + ", ".join(
                    f"{name}={selected_counts[name]}"
                    for name in selected_resnames
                )
                + "\n"
                f"Selected N_total  : {total_selected}\n"
                f"N_leaflet used    : {nleaflet:g}\n"
            )
            apl_stats.write_text(
                old_stats.rstrip() + "\n" + selection_text,
                encoding="utf-8",
            )
        except Exception as exc:
            print(
                f"[WARNING] Could not append APL selection metadata to stats: {exc}"
            )

        make_single_plot(
            apl_csv,
            APL_DIR / "apl_vs_time.png",
            "Time_ns",
            "APL_nm2_per_lipid",
            "Area per lipid vs Time",
            "Time (ns)",
            r"Area per lipid (nm$^2$/lipid)",
        )

        print(f"[RESULT] {apl_xvg}")
        print(f"[STATS ] {apl_stats}")
        print(f"[SELECT] {selection_file}")

        action = post_calculation_menu(
            "APL",
            "Calculate APL again (e.g. another membrane/system)",
        )
        if action == "continue":
            continue
        if action == "menu":
            return False
        return True



def run_thickness_menu(state: dict) -> bool:
    while True:
        print("\n======================================================================")
        print(" MEMBRANE THICKNESS CALCULATION")
        print("======================================================================")

        files = resolve_analysis_files(
            [
                ("gro", "gro", "GRO file"),
                ("tpr", "tpr", "TPR file"),
                ("xtc", "xtc", "XTC trajectory"),
            ],
            state,
        )
        gro, tpr, xtc = files["gro"], files["tpr"], files["xtc"]

        begin_ps, end_ps, dt_ps = prompt_time_window(require_dt=True)
        lipid_defs = prompt_thickness_lipid_definitions(gro, tpr)

        parsed = parse_gro(gro, lipid_defs)

        print("\n------------------------- STRUCTURE VALIDATION -------------------------")
        print(f"[OK] Total membrane lipids: {parsed['total_lipids']}")
        for resname, head in lipid_defs:
            print(
                f"     {resname}: {parsed['res_counts'][resname]} lipids; "
                f"head '{head}': {parsed['head_counts'][resname]} matches"
            )
        print(f"[OK] Upper leaflet heads: {len(parsed['upper_indices'])}")
        print(f"[OK] Lower leaflet heads: {len(parsed['lower_indices'])}")
        print(
            f"[OK] Initial headgroup separation: "
            f"{parsed['initial_thickness']:.4f} nm"
        )

        if (
            parsed["box_z"] > 0
            and parsed["initial_thickness"] > 0.5 * parsed["box_z"]
        ):
            print(
                "[WARNING] Apparent leaflet separation exceeds half Box-Z. "
                "The membrane may be split across Z PBC."
            )

        with LIPID_CONFIG.open("w", encoding="utf-8") as fh:
            fh.write("resname\thead_atom\n")
            for resname, head in lipid_defs:
                fh.write(f"{resname}\t{head}\n")

        build_index(gro, lipid_defs, parsed)

        thick_xvg, thick_csv, thick_stats = analyze_thickness(
            tpr,
            xtc,
            begin_ps,
            end_ps,
            dt_ps,
        )

        make_single_plot(
            thick_csv,
            THICK_DIR / "membrane_thickness_vs_time.png",
            "Time_ns",
            "Membrane_thickness_nm",
            "Membrane thickness vs Time",
            "Time (ns)",
            "Membrane thickness (nm)",
        )

        print(f"[RESULT] {thick_xvg}")
        print(f"[STATS ] {thick_stats}")

        action = post_calculation_menu(
            "MEMBRANE THICKNESS",
            "Calculate membrane thickness again (e.g. another membrane/system)",
        )
        if action == "continue":
            continue
        if action == "menu":
            return False
        return True


def run_sasa_menu(state: dict) -> bool:
    while True:
        print("\n======================================================================")
        print(" MEMBRANE SASA CALCULATION")
        print("======================================================================")

        files = resolve_analysis_files(
            [
                ("gro", "gro", "GRO file"),
                ("tpr", "tpr", "TPR file"),
                ("xtc", "xtc", "XTC trajectory"),
            ],
            state,
        )
        gro, tpr, xtc = files["gro"], files["tpr"], files["xtc"]

        begin_ps, end_ps, dt_ps = prompt_time_window(require_dt=True)

        print(
            "\nSelect ALL membrane lipid residue/molecule types by number. "
            "They will be merged into one [ membrane ] index group."
        )
        lipid_names = prompt_membrane_lipid_names(gro, tpr)
        build_membrane_only_index(gro, lipid_names)

        probe = ask_float("SASA probe radius in nm", 0.14)

        while True:
            try:
                ndots = int(ask("SASA dots per sphere", "24"))
                if ndots <= 0:
                    raise ValueError
                break
            except ValueError:
                print("[WARNING] Please enter a positive integer.")

        sasa_xvg, sasa_csv, sasa_stats = analyze_sasa(
            tpr,
            xtc,
            begin_ps,
            end_ps,
            dt_ps,
            probe,
            ndots,
        )

        make_single_plot(
            sasa_csv,
            SASA_DIR / "membrane_sasa_vs_time.png",
            "Time_ns",
            "Membrane_SASA_nm2",
            "Membrane SASA vs Time",
            "Time (ns)",
            r"SASA (nm$^2$)",
        )

        print(f"[RESULT] {sasa_xvg}")
        print(f"[STATS ] {sasa_stats}")

        action = post_calculation_menu(
            "MEMBRANE SASA",
            "Calculate membrane SASA again (e.g. another membrane/system)",
        )
        if action == "continue":
            continue
        if action == "menu":
            return False
        return True


def analyze_one_order_lipid(
    gro: Path,
    tpr: Path,
    xtc: Path,
    begin_ps: float,
    end_ps: float | None,
    dt_ps: float,
):
    """
    Calculate the order parameter for exactly ONE lipid residue.
    If the lipid is double-chain, chain1 and chain2 are completed
    continuously before this function returns.
    """
    print("\n------------------------- LIPID ORDER SETUP ----------------------------")

    selected = prompt_numbered_residue_selection(
        gro,
        title="ORDER PARAMETER — LIPID SELECTION",
        instruction=(
            "Select ONE lipid residue/molecule type whose chain order "
            "parameter will be calculated. After it finishes, the program "
            "can return here to calculate another lipid."
        ),
        allow_multiple=False,
        selection_name="lipid",
        tpr=tpr,
    )
    resname = selected[0]
    blocks = lipid_residue_blocks(gro, resname)

    print(
        f"[OK] Selected lipid: {resname} "
        f"({len(blocks)} molecule/residue blocks)."
    )

    while True:
        chain_raw = ask("Does this lipid have one chain or two chains? Enter 1 or 2")
        if chain_raw in {"1", "2"}:
            nchain = int(chain_raw)
            break
        print("[WARNING] Please enter 1 or 2.")

    ORDER_DIR.mkdir(parents=True, exist_ok=True)
    results = []

    if nchain == 1:
        atom_names = ask_chain_atom_names(gro, resname, "single-chain", tpr)
        base = f"{resname}-order"
        ndx = ORDER_DIR / f"{base}.ndx"
        xvg = ORDER_DIR / f"{base}.xvg"

        nmol = write_order_index(
            gro,
            resname,
            atom_names,
            ndx,
            "chain",
        )
        log = run_order_command(
            xtc,
            tpr,
            ndx,
            xvg,
            begin_ps,
            end_ps,
            dt_ps,
        )
        results.append(
            {
                "lipid": resname,
                "chain": "single",
                "molecules": nmol,
                "atoms": atom_names,
                "ndx": ndx,
                "xvg": xvg,
                "log": log,
            }
        )
    else:
        for chain_no in (1, 2):
            print(
                f"\n================ {resname}: BUILD CHAIN {chain_no} ================"
            )
            atom_names = ask_chain_atom_names(
                gro,
                resname,
                f"chain{chain_no}",
                tpr,
            )
            base = f"{resname}-order-chain{chain_no}"
            ndx = ORDER_DIR / f"{base}.ndx"
            xvg = ORDER_DIR / f"{base}.xvg"

            nmol = write_order_index(
                gro,
                resname,
                atom_names,
                ndx,
                f"chain{chain_no}",
            )
            log = run_order_command(
                xtc,
                tpr,
                ndx,
                xvg,
                begin_ps,
                end_ps,
                dt_ps,
            )
            results.append(
                {
                    "lipid": resname,
                    "chain": f"chain{chain_no}",
                    "molecules": nmol,
                    "atoms": atom_names,
                    "ndx": ndx,
                    "xvg": xvg,
                    "log": log,
                }
            )

        print(f"[OK] Both chains of {resname} have been completed.")

    summary_file = ORDER_DIR / "order_parameter_summary.xlsx"
    summary_cache = ORDER_DIR / ".order_parameter_summary_data.json"

    headers = [
        "lipid",
        "chain",
        "molecules",
        "sequential_atoms",
        "index_file",
        "order_file",
    ]

    for r in results:
        append_excel_summary(
            summary_file,
            summary_cache,
            "Order Parameter",
            headers,
            [
                r["lipid"],
                r["chain"],
                r["molecules"],
                " > ".join(r["atoms"]),
                str(r["ndx"]),
                str(r["xvg"]),
            ],
        )

    print(
        f"[OK] Order-parameter Excel summary updated: "
        f"{summary_file}"
    )
    return results


def run_order_menu(state: dict) -> bool:
    """
    File/time settings are entered before the first lipid.
    If the user chooses "calculate another lipid molecule", the current
    GRO/TPR/XTC/time settings can be reused without interruption.
    """
    print("\n======================================================================")
    print(" LIPID ORDER PARAMETER CALCULATION")
    print("======================================================================")
    print(
        "Method: gmx order -f XTC -s TPR -n lipid-order.ndx "
        "-o lipid-order.xvg -d z -szonly"
    )
    print()
    print("*** IMPORTANT GROMACS LIMITATION ***")
    print(
        "gmx order is documented for saturated carbons / united-atom style "
        "use. Do not automatically interpret its output as a validated "
        "Martini coarse-grained order parameter."
    )

    files = resolve_analysis_files(
        [
            ("gro", "gro", "GRO file"),
            ("tpr", "tpr", "TPR file"),
            ("xtc", "xtc", "XTC trajectory"),
        ],
        state,
    )
    gro, tpr, xtc = files["gro"], files["tpr"], files["xtc"]

    begin_ps, end_ps, dt_ps = prompt_time_window(require_dt=True)

    while True:
        results = analyze_one_order_lipid(
            gro,
            tpr,
            xtc,
            begin_ps,
            end_ps,
            dt_ps,
        )

        print("\n[ORDER RESULTS]")
        for r in results:
            print(f"  Index : {r['ndx']}")
            print(f"  Order : {r['xvg']}")
            print(f"  Log   : {r['log']}")

        action = post_calculation_menu(
            "LIPID ORDER PARAMETER",
            "Calculate the order parameter of another lipid molecule",
        )

        if action == "menu":
            return False
        if action == "exit":
            return True

        # Continue current parameter: normally reuse same system.
        reuse = yes_no(
            "Reuse the current GRO/TPR/XTC files and analysis time window?",
            True,
        )
        if reuse:
            continue

        files = resolve_analysis_files(
            [
                ("gro", "gro", "GRO file"),
                ("tpr", "tpr", "TPR file"),
                ("xtc", "xtc", "XTC trajectory"),
            ],
            state,
        )
        gro, tpr, xtc = files["gro"], files["tpr"], files["xtc"]
        begin_ps, end_ps, dt_ps = prompt_time_window(require_dt=True)



# ======================================================================
# 5) LIPID MSD / DIFFUSION COEFFICIENT
# ======================================================================

def sanitize_filename(text: str) -> str:
    """Create a filesystem-safe, human-readable label."""
    clean = re.sub(r"[^A-Za-z0-9_.-]+", "_", text.strip())
    clean = clean.strip("._")
    return clean or "selection"


def split_atom_names(raw: str) -> list[str]:
    """
    Accept atom/bead names separated by spaces, commas, semicolons or tabs.
    Preserve entry order and remove duplicates.
    """
    tokens = [x for x in re.split(r"[\s,;]+", raw.strip()) if x]
    unique = []
    seen = set()
    for token in tokens:
        if token not in seen:
            unique.append(token)
            seen.add(token)
    return unique


def atoms_for_resname(gro_file: Path, resname: str) -> list[int]:
    """Return all global 1-based atom indices belonging to one resname."""
    indices = []
    for block in read_gro_residue_blocks(gro_file):
        if block["resname"] == resname:
            indices.extend(idx for idx, _ in block["atoms"])
    return indices


def validate_region_atom_names(
    gro_file: Path,
    resname: str,
    atom_names: list[str],
    region_label: str,
) -> list[int]:
    """
    Select region atoms strictly as:
        selected lipid resname + selected atom/bead name

    Every listed name must occur exactly once in every residue/molecule of the
    selected lipid. This prevents same-named atoms from other lipids entering
    the membrane-region MSD groups.
    """
    blocks = lipid_residue_blocks(gro_file, resname)
    if not blocks:
        die(f"Residue '{resname}' was not found in {gro_file}.")

    if not atom_names:
        die(
            f"No atom/bead names were supplied for {resname} "
            f"{region_label} region."
        )

    selected = []
    problems = []

    for atom_name in atom_names:
        for mol_i, block in enumerate(blocks, start=1):
            matches = [idx for idx, name in block["atoms"] if name == atom_name]
            if len(matches) != 1:
                problems.append(
                    f"{resname} molecule #{mol_i}: '{atom_name}' matched "
                    f"{len(matches)} times in {region_label}; exactly 1 is required."
                )
            else:
                selected.append(matches[0])

    if problems:
        print(
            f"\n[ERROR] Validation failed for {resname} / {region_label}:",
            file=sys.stderr,
        )
        for problem in problems[:20]:
            print("  - " + problem, file=sys.stderr)
        if len(problems) > 20:
            print(
                f"  ... and {len(problems) - 20} additional errors",
                file=sys.stderr,
            )
        die(
            "Please correct the lipid-specific atom/bead names. "
            "Atom names are always restricted by the entered resname."
        )

    return selected


def prompt_region_names_for_lipid(
    gro_file: Path,
    resname: str,
    tpr: Path | None = None,
) -> dict[str, list[str]]:
    """
    Number-based lipid region definition.

    The user no longer types atom/bead names. Each region displays choices as:
        [1] CHOL-ROH  [2] CHOL-C1 ...
    """
    print("\n" + "-" * 78)
    print(
        f" MSD REGION DEFINITION FOR LIPID: {resname}"
    )
    print("-" * 78)

    result = {}

    region_info = [
        (
            "HEADGROUP",
            (
                "Select the polar/headgroup atoms or beads for this lipid."
            ),
        ),
        (
            "BARRIER",
            (
                "Select the interfacial/glycerol/sterol-ring/barrier-region "
                "atoms or beads for this lipid."
            ),
        ),
        (
            "TAIL",
            (
                "Select the hydrophobic tail/core atoms or beads for this lipid."
            ),
        ),
    ]

    for region, instruction in region_info:
        result[region] = prompt_numbered_atom_bead_selection(
            gro_file,
            resname,
            title=f"{region} REGION — {resname}",
            instruction=instruction,
            allow_multiple=True,
            min_count=1,
            tpr=tpr,
            require_once_per_residue=True,
        )

    return result



def write_clean_msd_index(
    ndx_file: Path,
    groups: list[tuple[str, list[int]]],
):
    """
    Create a dedicated MSD index from scratch.
    Only MSD calculation groups are written; this avoids ambiguity with
    default System/Protein/etc. groups and keeps output reproducible.
    """
    ndx_file.parent.mkdir(parents=True, exist_ok=True)

    with ndx_file.open("w", encoding="utf-8") as fh:
        for group_name, indices in groups:
            if not indices:
                die(f"MSD index group '{group_name}' is empty.")
            write_ndx_group(fh, group_name, indices)

    print(f"[OK] MSD index generated: {ndx_file}")
    for group_name, indices in groups:
        print(f"     [ {group_name} ] : {len(indices)} atoms/beads")


def get_gmx_msd_capabilities() -> dict[str, bool]:
    """
    Detect options at runtime for compatibility with both legacy GROMACS
    (e.g. 2020.6) and current selection-based gmx msd versions.
    """
    result = subprocess.run(
        [GMX, "msd", "-h"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    help_text = result.stdout or ""

    return {
        "sel": "-sel " in help_text or "\n-sel" in help_text,
        "dt": "-dt " in help_text or "\n-dt" in help_text,
        "maxtau": "-maxtau " in help_text or "\n-maxtau" in help_text,
    }


def prompt_msd_dimension_mode() -> tuple[str, list[str]]:
    """
    Return (human-readable description, gmx msd option list).

    Default 3D matches the user's example:
        gmx msd -f ... -s ... -n ... -o ...
    For membrane lateral diffusion, -lateral z computes motion in the XY plane.
    """
    print("\n------------------------- MSD DIMENSION MODE ---------------------------")
    print(" [1] 3D MSD / diffusion coefficient")
    print("     Matches ordinary 'gmx msd' without a direction flag.")
    print(" [2] Lateral XY diffusion of membrane lipids")
    print("     Uses: -lateral z")
    print(" [3] X-direction only")
    print(" [4] Y-direction only")
    print(" [5] Z-direction only")

    while True:
        choice = ask("Choose MSD dimensional mode", "1")

        if choice == "1":
            return "3D", []
        if choice == "2":
            return "Lateral_XY_perpendicular_to_Z", ["-lateral", "z"]
        if choice == "3":
            return "X_only", ["-type", "x"]
        if choice == "4":
            return "Y_only", ["-type", "y"]
        if choice == "5":
            return "Z_only", ["-type", "z"]

        print("[WARNING] Please choose 1, 2, 3, 4, or 5.")


def ask_fit_value(prompt: str, default: float = -1.0) -> float:
    """Allow -1, which is GROMACS' automatic 10%-90% fitting setting."""
    while True:
        raw = ask(prompt, str(default))
        try:
            value = float(raw)
            if value < 0 and value != -1:
                raise ValueError
            return value
        except ValueError:
            print("[WARNING] Enter -1 or a non-negative time value in ps.")


def prompt_msd_settings() -> dict:
    """
    Ask only the settings that materially affect MSD/D fitting.
    Simple users can retain GROMACS defaults.
    """
    begin_ps, end_ps, dt_ps = prompt_time_window(require_dt=True)
    dimension_label, dimension_args = prompt_msd_dimension_mode()

    print("\n-------------------------- MSD FIT SETTINGS ----------------------------")
    print(
        "The diffusion coefficient is obtained from a linear fit to MSD(t). "
        "GROMACS uses 10%-90% of the MSD time range when beginfit/endfit = -1."
    )

    use_defaults = yes_no(
        "Use default MSD fitting settings (trestart=10 ps, beginfit=-1, endfit=-1)?",
        True,
    )

    if use_defaults:
        trestart = 10.0
        beginfit = -1.0
        endfit = -1.0
        maxtau = None
    else:
        trestart = ask_float("Time between MSD restart/reference points (-trestart), ps", 10.0)
        beginfit = ask_fit_value("Diffusion fit start (-beginfit), ps; -1 = automatic", -1.0)
        endfit = ask_fit_value("Diffusion fit end (-endfit), ps; -1 = automatic", -1.0)

        while True:
            raw = ask(
                "Maximum MSD time lag (-maxtau), ps; press Enter for no limit",
                "",
            )

            if not raw.strip():
                maxtau = None
                break

            try:
                maxtau = float(raw)
                if maxtau <= 0:
                    raise ValueError
                break
            except ValueError:
                print(
                    "[WARNING] -maxtau must be greater than 0 ps, "
                    "or press Enter for no limit."
                )

    return {
        "begin_ps": begin_ps,
        "end_ps": end_ps,
        "dt_ps": dt_ps,
        "dimension_label": dimension_label,
        "dimension_args": dimension_args,
        "trestart": trestart,
        "beginfit": beginfit,
        "endfit": endfit,
        "maxtau": maxtau,
    }


def diffusion_pattern_matches(text: str):
    """
    Extract GROMACS diffusion legend/output lines such as:
      D[       OST] = 0.8069 (+/- 0.0496) (1e-5 cm^2/s)
    Also tolerates versions that omit '=' or the outer unit parentheses.
    """
    number = r"(?:[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?|[-+]?nan)"
    pattern = re.compile(
        rf"D\[\s*(?P<label>[^\]]+?)\s*\]\s*=?\s*"
        rf"(?P<D>{number})\s*"
        rf"\(\+/-\s*(?P<err>{number})\s*\)\s*"
        rf"\(?\s*(?P<unit>1e-5\s+cm\^2/s)\s*\)?",
        re.IGNORECASE,
    )
    return list(pattern.finditer(text))


def extract_diffusion_coefficient(
    xvg_file: Path,
    log_file: Path,
    requested_group: str,
) -> dict:
    """
    Search both the XVG legend and command log for the diffusion coefficient.
    """
    texts = []

    if xvg_file.exists():
        texts.append(xvg_file.read_text(encoding="utf-8", errors="replace"))
    if log_file.exists():
        texts.append(log_file.read_text(encoding="utf-8", errors="replace"))

    matches = []
    for source_text in texts:
        matches.extend(diffusion_pattern_matches(source_text))

    if not matches:
        return {
            "group": requested_group,
            "D": None,
            "error": None,
            "unit": "1e-5 cm^2/s",
            "raw": "NOT_FOUND",
        }

    # Prefer a label matching the requested group; otherwise use the last match.
    chosen = None
    for m in matches:
        label = m.group("label").strip()
        if label == requested_group:
            chosen = m
    if chosen is None:
        chosen = matches[-1]

    return {
        "group": chosen.group("label").strip(),
        "D": chosen.group("D"),
        "error": chosen.group("err"),
        "unit": re.sub(r"\s+", " ", chosen.group("unit")).strip(),
        "raw": chosen.group(0),
    }


def run_msd_command(
    xtc: Path,
    tpr: Path,
    ndx_file: Path,
    group_name: str,
    xvg_file: Path,
    settings: dict,
) -> tuple[Path, dict]:
    """
    Run one gmx msd calculation.

    Compatibility strategy:
    - current GROMACS: use -sel 'group "NAME"'
    - older GROMACS (e.g. 2020.6): select NAME interactively via stdin
    - pass -dt / -maxtau only when the installed gmx msd supports them
    """
    caps = get_gmx_msd_capabilities()

    cmd = [
        GMX,
        "msd",
        "-f",
        str(xtc),
        "-s",
        str(tpr),
        "-n",
        str(ndx_file),
        "-o",
        str(xvg_file),
        "-b",
        str(settings["begin_ps"]),
    ]

    if settings["end_ps"] is not None:
        cmd += ["-e", str(settings["end_ps"])]

    if settings["dt_ps"] > 0:
        if caps["dt"]:
            cmd += ["-dt", str(settings["dt_ps"])]
        else:
            print(
                "[WARNING] Installed gmx msd does not support -dt; "
                "the requested frame interval will not be passed."
            )

    cmd += [
        "-trestart",
        str(settings["trestart"]),
        "-beginfit",
        str(settings["beginfit"]),
        "-endfit",
        str(settings["endfit"]),
    ]

    if settings["maxtau"] is not None:
        if caps["maxtau"]:
            cmd += ["-maxtau", str(settings["maxtau"])]
        else:
            print(
                "[WARNING] Installed gmx msd does not support -maxtau; "
                "this setting will be ignored."
            )

    cmd += settings["dimension_args"]

    if caps["sel"]:
        cmd += ["-sel", f'group "{group_name}"']
        input_text = None
    else:
        # Legacy index-group selection, compatible with gmx msd 2020.x.
        input_text = group_name + "\n"

    log_file = timestamped_log_path(xvg_file.with_suffix(".log"))

    result = _animated_process(
        cmd,
        input_text=input_text,
        label=f"MSD / diffusion: {group_name}",
    )

    log_file.write_text(
        log_timestamp_header() + "COMMAND\n"
        + " ".join(cmd)
        + "\n\nSELECTION\n"
        + group_name
        + "\n\nSTDOUT\n"
        + (result.stdout or "")
        + "\nSTDERR\n"
        + (result.stderr or ""),
        encoding="utf-8",
    )

    if result.returncode != 0:
        print(result.stdout or "")
        print(result.stderr or "", file=sys.stderr)
        die(
            f"gmx msd failed for group '{group_name}'. "
            f"Full command log: {log_file}"
        )

    if not xvg_file.exists() or xvg_file.stat().st_size == 0:
        die(
            f"gmx msd did not generate the requested MSD curve: {xvg_file}"
        )

    diffusion = extract_diffusion_coefficient(
        xvg_file,
        log_file,
        group_name,
    )

    print(f"[OK] MSD curve generated: {xvg_file}")
    if diffusion["D"] is not None:
        print(
            f"[OK] D[{diffusion['group']}] = {diffusion['D']} "
            f"(+/- {diffusion['error']}) ({diffusion['unit']})"
        )
    else:
        print(
            "[WARNING] Could not automatically extract D from the XVG/log. "
            f"Inspect: {log_file}"
        )

    return log_file, diffusion


def write_diffusion_result_file(
    result_file: Path,
    object_type: str,
    object_name: str,
    region: str,
    group_name: str,
    diffusion: dict,
    ndx_file: Path,
    xvg_file: Path,
    log_file: Path,
    settings: dict,
):
    """Write one clear human-readable diffusion result file."""
    lines = [
        "======================================================================",
        "GROMACS MSD / DIFFUSION COEFFICIENT RESULT",
        "======================================================================",
        f"Calculation target type : {object_type}",
        f"Calculation target      : {object_name}",
        f"Region             : {region}",
        f"Index group        : {group_name}",
        f"MSD dimension      : {settings['dimension_label']}",
        "",
        f"Trajectory begin   : {settings['begin_ps']} ps",
        f"Trajectory end     : "
        + (
            f"{settings['end_ps']} ps"
            if settings["end_ps"] is not None
            else "trajectory end"
        ),
        f"Sampling interval  : {settings['dt_ps']} ps",
        f"trestart           : {settings['trestart']} ps",
        f"beginfit           : {settings['beginfit']} ps",
        f"endfit             : {settings['endfit']} ps",
        f"maxtau             : "
        + (
            f"{settings['maxtau']} ps"
            if settings["maxtau"] is not None
            else "not set"
        ),
        "",
    ]

    if diffusion["D"] is not None:
        lines += [
            f"Diffusion D         : {diffusion['D']} ({diffusion['unit']})",
            f"Estimated error     : +/- {diffusion['error']} ({diffusion['unit']})",
            f"GROMACS legend      : {diffusion['raw']}",
        ]
    else:
        lines += [
            "Diffusion D         : NOT AUTOMATICALLY EXTRACTED",
            f"GROMACS legend      : {diffusion['raw']}",
        ]

    lines += [
        "",
        f"Index file          : {ndx_file}",
        f"MSD curve           : {xvg_file}",
        f"Command log         : {log_file}",
        "======================================================================",
        "",
    ]

    result_file.write_text("\n".join(lines), encoding="utf-8")



def append_global_msd_summary(
    object_type: str,
    object_name: str,
    region: str,
    group_name: str,
    diffusion: dict,
    ndx_file: Path,
    xvg_file: Path,
    result_file: Path,
    settings: dict,
):
    """Append one MSD/D record to an Excel summary workbook."""
    summary = MSD_DIR / "MSD_diffusion_summary.xlsx"
    cache = MSD_DIR / ".MSD_diffusion_summary_data.json"

    headers = [
        "calculation_target_type",
        "calculation_target",
        "region",
        "group",
        "D_1e-5_cm2_s",
        "error_1e-5_cm2_s",
        "dimension",
        "trestart_ps",
        "beginfit_ps",
        "endfit_ps",
        "index_file",
        "msd_curve",
        "result_file",
    ]

    row = [
        object_type,
        object_name,
        region,
        group_name,
        diffusion["D"]
        if diffusion["D"] is not None
        else "NA",
        diffusion["error"]
        if diffusion["error"] is not None
        else "NA",
        settings["dimension_label"],
        settings["trestart"],
        settings["beginfit"],
        settings["endfit"],
        str(ndx_file),
        str(xvg_file),
        str(result_file),
    ]

    append_excel_summary(
        summary,
        cache,
        "MSD Diffusion",
        headers,
        row,
    )


def make_msd_plot(xvg_file: Path, png_file: Path, title: str):
    """Optional PNG copy of the GROMACS MSD curve."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return

    rows = parse_xvg(xvg_file, min_cols=2)
    if not rows:
        return

    x = [r[0] for r in rows]
    y = [r[1] for r in rows]

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(x, y, linewidth=1.2)
    ax.set_title(title)
    ax.set_xlabel("MSD lag time (ps)")
    ax.set_ylabel(r"MSD (nm$^2$)")
    fig.tight_layout()
    fig.savefig(png_file, dpi=300)
    plt.close(fig)
    print(f"[OK] MSD plot generated: {png_file}")




def prompt_multiple_resnames_for_msd(
    gro: Path,
    tpr: Path | None = None,
) -> list[str]:
    """
    Select one or multiple residue/molecule TYPES for separate MSD/D
    calculations using numbered choices.

    This does NOT mean one individual molecule instance. If CHOL is selected,
    the CHOL group contains all GRO residues/molecules whose resname is CHOL.
    """
    return prompt_numbered_residue_selection(
        gro,
        title="MSD / DIFFUSION — RESIDUE/MOLECULE TYPE SELECTION",
        instruction=(
            "Select one or more residue/molecule TYPES to calculate. "
            "Each selected type will be calculated separately and will "
            "receive its own MSD curve and diffusion coefficient. "
            "For example, selecting CHOL and HSPC calculates one CHOL MSD/D "
            "and one HSPC MSD/D."
        ),
        allow_multiple=True,
        selection_name="MSD target type",
        tpr=tpr,
    )



def calculate_single_species_msd(
    gro: Path,
    tpr: Path,
    xtc: Path,
    settings: dict,
):
    """
    MSD mode 1: one or multiple molecular species selected by resname.

    All selected resnames are entered in one step. One clean index is built
    containing one independent group per resname. For compatibility with both
    legacy and current GROMACS versions, gmx msd is then executed sequentially
    for each selected group using the same files and fitting settings.
    """
    print("\n======================================================================")
    print(" MSD MODE 1 — SELECTED RESIDUE/MOLECULE TYPES")
    print("======================================================================")
    print(
        "Select one or more residue/molecule TYPES by number. "
        "Each selected type is calculated separately and receives its own "
        "MSD curve and diffusion coefficient."
    )

    resnames = prompt_multiple_resnames_for_msd(gro, tpr)

    batch_dir = MSD_DIR / "Single_Resnames_Batch"
    batch_dir.mkdir(parents=True, exist_ok=True)

    groups = []
    group_names = {}

    for resname in resnames:
        safe = sanitize_filename(resname)
        group_name = f"MSD_{safe}"
        indices = atoms_for_resname(gro, resname)
        if not indices:
            die(f"No atoms/beads found for resname '{resname}'.")
        groups.append((group_name, indices))
        group_names[resname] = group_name

    ndx_file = batch_dir / "MSD_single_resnames_batch_index.ndx"
    write_clean_msd_index(ndx_file, groups)

    batch_summary = batch_dir / "MSD_single_resnames_batch_diffusion_summary.xlsx"
    results = []

    print("\n======================================================================")
    print(f" STARTING MSD BATCH: {len(resnames)} RESNAME(S)")
    print("======================================================================")

    for i, resname in enumerate(resnames, start=1):
        safe = sanitize_filename(resname)
        group_name = group_names[resname]

        print(f"\n---------------- MSD {i}/{len(resnames)} : {resname} ----------------")

        xvg_file = batch_dir / f"MSD_single_{safe}_curve.xvg"
        result_file = batch_dir / f"MSD_single_{safe}_diffusion_result.txt"
        png_file = batch_dir / f"MSD_single_{safe}_curve.png"

        log_file, diffusion = run_msd_command(
            xtc,
            tpr,
            ndx_file,
            group_name,
            xvg_file,
            settings,
        )

        write_diffusion_result_file(
            result_file,
            object_type="molecular species selected by resname (batch mode)",
            object_name=resname,
            region="whole selected molecular species",
            group_name=group_name,
            diffusion=diffusion,
            ndx_file=ndx_file,
            xvg_file=xvg_file,
            log_file=log_file,
            settings=settings,
        )

        append_global_msd_summary(
            "batch_resname",
            resname,
            "whole_species",
            group_name,
            diffusion,
            ndx_file,
            xvg_file,
            result_file,
            settings,
        )

        make_msd_plot(
            xvg_file,
            png_file,
            f"MSD — {resname}",
        )

        results.append(
            {
                "resname": resname,
                "group": group_name,
                "D": diffusion["D"],
                "error": diffusion["error"],
                "unit": diffusion["unit"],
                "curve": xvg_file,
                "result": result_file,
                "log": log_file,
            }
        )

        print(
            f"[BATCH] Completed {i}/{len(resnames)} selected resname(s): {resname}"
        )

    batch_headers = [
        "resname",
        "group",
        "D_1e-5_cm2_s",
        "error_1e-5_cm2_s",
        "dimension",
        "msd_curve",
        "result_file",
        "log_file",
    ]
    batch_rows = [
        [
            row["resname"],
            row["group"],
            row["D"] if row["D"] is not None else "NA",
            row["error"] if row["error"] is not None else "NA",
            settings["dimension_label"],
            str(row["curve"]),
            str(row["result"]),
            str(row["log"]),
        ]
        for row in results
    ]
    write_excel_table(
        batch_summary,
        "MSD Batch Summary",
        batch_headers,
        batch_rows,
    )

    print("\n======================================================================")
    print(" MULTI-RESNAME MSD / DIFFUSION SUMMARY")
    print("======================================================================")
    for row in results:
        if row["D"] is not None:
            print(
                f" {row['resname']:<14s}: "
                f"D = {row['D']} (+/- {row['error']}) ({row['unit']})"
            )
        else:
            print(
                f" {row['resname']:<14s}: "
                "D = NOT AUTOMATICALLY EXTRACTED"
            )

    print()
    print(f"  Batch index   : {ndx_file}")
    print(f"  Batch summary : {batch_summary}")
    print(f"  Output folder : {batch_dir}")

    return results


def selected_lipid_atom_catalog(
    gro_file: Path,
    lipid_names: list[str],
) -> dict[str, list[str]]:
    """Collect unique atom/bead names for each selected lipid species."""
    catalog = {}

    for resname in lipid_names:
        blocks = lipid_residue_blocks(gro_file, resname)
        if not blocks:
            die(f"Residue '{resname}' was not found in {gro_file}.")

        names = []
        seen = set()
        for _, atomname in blocks[0]["atoms"]:
            if atomname not in seen:
                seen.add(atomname)
                names.append(atomname)

        catalog[resname] = names

    return catalog


def show_selected_lipid_atom_catalog(
    gro_file: Path,
    lipid_names: list[str],
):
    """
    Show all selected membrane lipids and their atoms/beads with NUMBER labels.

    Example:
        [1] CHOL-ROH  [2] CHOL-C1 ...
    """
    counts = residue_counts_from_gro(
        gro_file
    )

    print("\n" + "=" * 78)
    print(" SELECTED MEMBRANE LIPIDS — NUMBERED ATOM/BEAD CATALOG")
    print("=" * 78)

    catalog = {}

    for lipid_no, resname in enumerate(
        lipid_names,
        start=1,
    ):
        print(
            f"\nLipid [{lipid_no}] {resname} "
            f"({counts[resname]} molecule/residue blocks)"
        )

        catalog[resname] = show_numbered_atom_bead_items(
            gro_file,
            resname,
            title=f"{resname} ATOM / BEAD CATALOG",
        )

    return catalog



def _parse_region_batch_entry(
    raw: str,
    lipid_names: list[str],
) -> tuple[dict[str, list[str]], list[str]]:
    """
    Parse one-line region definitions covering all selected membrane lipids.

    Supported examples:
      1:NC3,PO4; 2:ROH; 3:NC3,PO4
      HSPC:NC3 PO4 | CHOL:ROH | DSPE:NC3 PO4
    """
    definitions = {}
    errors = []

    entries = [
        item.strip()
        for item in re.split(r"[;|]+", raw.strip())
        if item.strip()
    ]

    if not entries:
        return {}, ["No lipid-region definitions were entered."]

    for entry in entries:
        if ":" not in entry:
            errors.append(
                f"'{entry}' has no ':' separator. "
                "Use a format such as 1:NC3,PO4."
            )
            continue

        key, atom_text = entry.split(":", 1)
        key = key.strip()
        atom_text = atom_text.strip()

        if not key:
            errors.append(
                f"Missing lipid number/resname in '{entry}'."
            )
            continue

        resname = None

        try:
            number = int(key)
            if 1 <= number <= len(lipid_names):
                resname = lipid_names[number - 1]
            else:
                errors.append(
                    f"Lipid number '{key}' is outside 1-{len(lipid_names)}."
                )
                continue
        except ValueError:
            if key in lipid_names:
                resname = key
            else:
                errors.append(
                    f"'{key}' is neither a valid lipid number nor a "
                    "selected resname."
                )
                continue

        if resname in definitions:
            errors.append(
                f"Lipid '{resname}' was defined more than once."
            )
            continue

        atom_names = split_atom_names(atom_text)

        if not atom_names:
            errors.append(
                f"No atom/bead names were provided for '{resname}'."
            )
            continue

        definitions[resname] = atom_names

    missing = [
        resname
        for resname in lipid_names
        if resname not in definitions
    ]

    if missing:
        errors.append(
            "Missing definition(s) for selected lipid(s): "
            + ", ".join(missing)
        )

    return definitions, errors


def prompt_whole_membrane_region_batch(
    gro_file: Path,
    lipid_names: list[str],
    region: str,
    tpr: Path | None = None,
) -> dict[str, list[str]]:
    """
    Define ONE structural region across all selected membrane lipids using
    numbered atom/bead choices only.

    Previous manual input such as:
        HSPC:PO4,NC3; CHOL:ROH
    is replaced by direct numbered selections for each lipid.
    """
    descriptions = {
        "HEADGROUP": "polar / hydrophilic headgroup atoms or beads",
        "BARRIER": (
            "interfacial / glycerol / sterol-ring / barrier-region "
            "atoms or beads"
        ),
        "TAIL": "hydrophobic tail-chain/core atoms or beads",
    }

    print("\n" + "-" * 78)
    print(
        f" DEFINE {region} REGION FOR ALL SELECTED MEMBRANE LIPIDS"
    )
    print("-" * 78)
    print(
        f"Region meaning: {descriptions.get(region, region)}."
    )
    print(
        "For each lipid, select the corresponding atoms/beads by NUMBER."
    )

    definitions = {}

    for number, resname in enumerate(
        lipid_names,
        start=1,
    ):
        print(
            f"\nLipid {number}/{len(lipid_names)}: {resname}"
        )

        definitions[resname] = prompt_numbered_atom_bead_selection(
            gro_file,
            resname,
            title=f"{region} — {resname}",
            instruction=(
                f"Select one or more {region} atoms/beads for {resname}."
            ),
            allow_multiple=True,
            min_count=1,
            tpr=tpr,
            require_once_per_residue=True,
        )

    print(
        f"\n[{region} REGION CONFIRMED]"
    )

    for number, resname in enumerate(
        lipid_names,
        start=1,
    ):
        print(
            f"  [{number}] {resname:<12s}: "
            + ", ".join(
                f"{resname}-{name}"
                for name in definitions[resname]
            )
        )

    return definitions



def prompt_all_membrane_regions_batch(
    gro_file: Path,
    lipid_names: list[str],
    tpr: Path | None = None,
) -> dict[str, dict[str, list[str]]]:
    """
    Define all whole-membrane MSD regions in only THREE prompts:
      1) HEADGROUP for all selected lipids
      2) BARRIER   for all selected lipids
      3) TAIL      for all selected lipids
    """
    show_selected_lipid_atom_catalog(
        gro_file,
        lipid_names,
    )

    region_map = {}

    for region in (
        "HEADGROUP",
        "BARRIER",
        "TAIL",
    ):
        region_map[region] = prompt_whole_membrane_region_batch(
            gro_file,
            lipid_names,
            region,
            tpr,
        )

    definitions = {
        resname: {}
        for resname in lipid_names
    }

    for region in (
        "HEADGROUP",
        "BARRIER",
        "TAIL",
    ):
        for resname in lipid_names:
            definitions[resname][region] = (
                region_map[region][resname]
            )

    print("\n" + "=" * 78)
    print(" WHOLE-MEMBRANE REGION DEFINITIONS COMPLETED")
    print("=" * 78)

    for region in (
        "HEADGROUP",
        "BARRIER",
        "TAIL",
    ):
        print(f"\n {region}:")
        for number, resname in enumerate(
            lipid_names,
            start=1,
        ):
            print(
                f"   [{number}] {resname:<12s}: "
                + ", ".join(definitions[resname][region])
            )

    print("=" * 78)

    return definitions




def calculate_membrane_region_msd(
    gro: Path,
    tpr: Path,
    xtc: Path,
    settings: dict,
):
    """
    MSD mode 2 — whole membrane divided into structural regions.

    Interaction:
      1) select all membrane lipid types by number;
      2) show all selected lipid atom/bead catalogs once;
      3) enter HEADGROUP mappings for all lipids in one line;
      4) enter BARRIER mappings for all lipids in one line;
      5) enter TAIL mappings for all lipids in one line;
      6) combine all selected atoms/beads into three membrane-wide groups.
    """
    print("\n======================================================================")
    print(" MSD MODE 2 — WHOLE MEMBRANE DIVIDED INTO STRUCTURAL REGIONS")
    print("======================================================================")
    print(
        "The whole membrane will be analyzed as three independent MSD/D "
        "regions:"
    )
    print("  [ MSD_MEMBRANE_HEADGROUP ]")
    print("  [ MSD_MEMBRANE_BARRIER ]")
    print("  [ MSD_MEMBRANE_TAIL ]")
    print()
    print(
        "After selecting membrane lipid types, each region is configured "
        "for ALL selected lipids in one line."
    )

    lipid_names = prompt_membrane_lipid_names(gro, tpr)

    # New interaction: only three region-level prompts.
    definitions = prompt_all_membrane_regions_batch(
        gro,
        lipid_names,
        tpr,
    )

    region_indices = {
        "HEADGROUP": [],
        "BARRIER": [],
        "TAIL": [],
    }

    for region in (
        "HEADGROUP",
        "BARRIER",
        "TAIL",
    ):
        for resname in lipid_names:
            indices = validate_region_atom_names(
                gro,
                resname,
                definitions[resname][region],
                region,
            )
            region_indices[region].extend(indices)

    subdir = MSD_DIR / "Whole_Membrane_Regions"
    subdir.mkdir(parents=True, exist_ok=True)

    ndx_file = subdir / "MSD_membrane_regions_index.ndx"

    groups = [
        (
            "MSD_MEMBRANE_HEADGROUP",
            region_indices["HEADGROUP"],
        ),
        (
            "MSD_MEMBRANE_BARRIER",
            region_indices["BARRIER"],
        ),
        (
            "MSD_MEMBRANE_TAIL",
            region_indices["TAIL"],
        ),
    ]

    write_clean_msd_index(
        ndx_file,
        groups,
    )

    definition_file = (
        subdir
        / "MSD_membrane_region_definitions.xlsx"
    )

    definition_headers = [
        "lipid",
        "region",
        "atom_or_bead_names",
    ]

    definition_rows = []

    for region in (
        "HEADGROUP",
        "BARRIER",
        "TAIL",
    ):
        for resname in lipid_names:
            definition_rows.append(
                [
                    resname,
                    region,
                    " ".join(
                        definitions[resname][region]
                    ),
                ]
            )

    write_excel_table(
        definition_file,
        "Region Definitions",
        definition_headers,
        definition_rows,
    )

    summary_rows = []

    for region, group_name in [
        (
            "HEADGROUP",
            "MSD_MEMBRANE_HEADGROUP",
        ),
        (
            "BARRIER",
            "MSD_MEMBRANE_BARRIER",
        ),
        (
            "TAIL",
            "MSD_MEMBRANE_TAIL",
        ),
    ]:
        region_lower = region.lower()

        xvg_file = (
            subdir
            / f"MSD_membrane_{region_lower}_curve.xvg"
        )
        result_file = (
            subdir
            / f"MSD_membrane_{region_lower}_diffusion_result.txt"
        )
        png_file = (
            subdir
            / f"MSD_membrane_{region_lower}_curve.png"
        )

        log_file, diffusion = run_msd_command(
            xtc,
            tpr,
            ndx_file,
            group_name,
            xvg_file,
            settings,
        )

        write_diffusion_result_file(
            result_file,
            object_type="whole membrane structural region",
            object_name="all specified membrane lipids",
            region=region,
            group_name=group_name,
            diffusion=diffusion,
            ndx_file=ndx_file,
            xvg_file=xvg_file,
            log_file=log_file,
            settings=settings,
        )

        append_global_msd_summary(
            "whole_membrane_region",
            "all_specified_membrane_lipids",
            region,
            group_name,
            diffusion,
            ndx_file,
            xvg_file,
            result_file,
            settings,
        )

        make_msd_plot(
            xvg_file,
            png_file,
            f"Membrane MSD — {region}",
        )

        summary_rows.append(
            {
                "region": region,
                "group": group_name,
                "D": diffusion["D"],
                "error": diffusion["error"],
                "unit": diffusion["unit"],
                "curve": xvg_file,
                "result": result_file,
            }
        )

    region_summary = (
        subdir
        / "MSD_membrane_regions_diffusion_summary.xlsx"
    )

    region_headers = [
        "region",
        "group",
        "D_1e-5_cm2_s",
        "error_1e-5_cm2_s",
        "dimension",
        "msd_curve",
        "result_file",
    ]

    region_rows = [
        [
            row["region"],
            row["group"],
            (
                row["D"]
                if row["D"] is not None
                else "NA"
            ),
            (
                row["error"]
                if row["error"] is not None
                else "NA"
            ),
            settings["dimension_label"],
            str(row["curve"]),
            str(row["result"]),
        ]
        for row in summary_rows
    ]

    write_excel_table(
        region_summary,
        "Membrane Regions",
        region_headers,
        region_rows,
    )

    print("\n======================================================================")
    print(" WHOLE-MEMBRANE MSD / DIFFUSION SUMMARY")
    print("======================================================================")

    for row in summary_rows:
        if row["D"] is not None:
            print(
                f" {row['region']:<10s}: "
                f"D = {row['D']} "
                f"(+/- {row['error']}) "
                f"({row['unit']})"
            )
        else:
            print(
                f" {row['region']:<10s}: "
                "D = NOT EXTRACTED"
            )

    print()
    print(
        f"  Region definitions : {definition_file}"
    )
    print(
        f"  Region index       : {ndx_file}"
    )
    print(
        f"  Summary workbook   : {region_summary}"
    )





# ======================================================================
# MSD EXTENSION — SPECIFIC ATOM/BEAD + TARGET MOLECULE IN MEMBRANE REGIONS
# ======================================================================
#
# IMPORTANT DESIGN RULE
# ---------------------
# This section ADDS functions to the original MSD module. It does not remove
# or replace the existing calculations:
#
#   Existing MSD mode 1:
#       selected residue/molecule TYPES
#
#   Existing MSD mode 2:
#       whole membrane divided into HEADGROUP / BARRIER / TAIL regions
#
# New MSD mode 3:
#       a specific atom (AA) or bead (CG) inside a selected residue/molecule
#
# New MSD mode 4:
#       a selected molecule/residue type in dynamically identified membrane
#       HEADGROUP / BARRIER / TAIL regions
#
# Main-menu module 7 (dynamic drug-region diffusion) is also retained exactly
# as an independent advanced workflow.
# ======================================================================


def show_numbered_atom_bead_items_for_msd(
    gro: Path,
    resname: str,
) -> list[str]:
    """
    Show the atom/bead names of one selected residue/molecule as a numbered
    list, so users do not need to type names manually.

    The first residue block is used as the topology/name template. The
    selected names are later validated against all residue blocks and TPR.
    """
    blocks = lipid_residue_blocks(
        gro,
        resname,
    )

    if not blocks:
        raise RecoverableAnalysisError(
            f"Residue/molecule type '{resname}' was not found in {gro}."
        )

    names = []
    seen = set()

    for _, atomname in blocks[0]["atoms"]:
        if atomname not in seen:
            seen.add(atomname)
            names.append(atomname)

    if not names:
        raise RecoverableAnalysisError(
            f"No atom/bead names could be read for '{resname}'."
        )

    print("\n" + "-" * 78)
    print(f" MSD SPECIFIC ATOM / BEAD SELECTION — {resname}")
    print("-" * 78)
    print(
        "Choose by NUMBER. The same interface is used for all-atom atoms "
        "and coarse-grained beads."
    )
    print()

    # More compact than printing one item per line while still remaining
    # readable for large AA molecules.
    per_line = 4
    cells = []

    for number, name in enumerate(
        names,
        start=1,
    ):
        cells.append(
            f"[{number:>3d}] {name:<12s}"
        )

    for i in range(
        0,
        len(cells),
        per_line,
    ):
        print(
            "  " + "   ".join(
                cells[i:i + per_line]
            )
        )

    print("-" * 78)

    return names


def prompt_specific_atom_bead_names_for_msd(
    gro: Path,
    tpr: Path,
    resname: str,
) -> list[str]:
    """
    Select one or multiple atom/bead names by number.

    Examples:
        AA : P, O3, C12 ...
        CG : PO4, ROH, C1A ...

    Every selected name must occur exactly once in every residue/molecule
    block of the selected resname. This makes a combined site-MSD physically
    interpretable and prevents accidental duplicate-name selections.
    """
    while True:
        names = show_numbered_atom_bead_items_for_msd(
            gro,
            resname,
        )

        print()
        print("Selection examples:")
        print("  One site       : 5")
        print("  Several sites  : 5 8 11")
        print("  Comma form     : 5,8,11")
        print()

        raw = ask(
            "Select atom/bead number(s)"
        )

        numbers, invalid = _parse_numbered_selection(
            raw,
            len(names),
        )

        if invalid or not numbers:
            if invalid:
                print(
                    "[WARNING] Invalid atom/bead number(s): "
                    + ", ".join(invalid)
                )
            print(
                f"[ACTION] Enter one or more numbers from 1-{len(names)}."
            )
            continue

        selected_names = [
            names[number - 1]
            for number in numbers
        ]

        validation_errors = []

        for atomname in selected_names:
            try:
                validate_region_atom_names(
                    gro,
                    resname,
                    [atomname],
                    "SPECIFIC_SITE_MSD",
                )
            except RecoverableAnalysisError as exc:
                validation_errors.append(
                    str(exc)
                )

            ok_tpr, detail_tpr = validate_atom_in_tpr(
                tpr,
                resname,
                atomname,
            )

            if not ok_tpr:
                validation_errors.append(
                    detail_tpr
                )

        if validation_errors:
            print()
            print(
                "[WARNING] One or more selected atom/bead names failed "
                "GRO/TPR validation:"
            )
            for error in validation_errors:
                print(
                    "  - " + error
                )
            print(
                "[ACTION] Please make the atom/bead selection again."
            )
            continue

        print("\n[SELECTED SPECIFIC ATOM/BEAD TARGETS]")
        for number in numbers:
            print(
                f"  [{number:>3d}] {names[number - 1]}"
            )

        if yes_no(
            "Confirm these atom/bead target(s)?",
            True,
        ):
            return selected_names


def prompt_specific_site_scope_for_msd(
    gro: Path,
    resname: str,
) -> tuple[str, list[int]]:
    """
    Decide whether a selected atom/bead is analyzed:
      1) as one combined group across ALL molecules of this resname,
      2) in selected individual molecule/residue instances separately,
      3) in ALL individual instances separately.

    Returned occurrence numbers are 1-based within the selected resname.
    """
    blocks = lipid_residue_blocks(
        gro,
        resname,
    )

    if not blocks:
        raise RecoverableAnalysisError(
            f"No residue/molecule blocks named '{resname}' were found."
        )

    print("\n" + "-" * 78)
    print(" SPECIFIC ATOM / BEAD — MOLECULE SCOPE")
    print("-" * 78)
    print(
        f"Detected {len(blocks)} individual '{resname}' "
        "molecule/residue occurrence(s)."
    )
    print()
    print(" [1] ALL molecules combined")
    print(
        "     The selected atom/bead from every molecule is combined into one "
        "MSD group. Output = ensemble-average site diffusion."
    )
    print()
    print(" [2] SELECTED individual molecule/residue instance(s)")
    print(
        "     Choose occurrence numbers. Each selected molecule receives its "
        "own independent MSD curve and diffusion coefficient."
    )
    print()
    print(" [3] ALL individual molecule/residue instances separately")
    print(
        "     Every molecule receives its own MSD curve and D value. "
        "Use with care for very large molecule counts."
    )

    while True:
        mode = ask(
            "Choose molecule scope",
            "1",
        )

        if mode in {
            "1",
            "2",
            "3",
        }:
            break

        print(
            "[WARNING] Please choose 1, 2, or 3."
        )

    if mode == "1":
        return (
            "all_combined",
            list(
                range(
                    1,
                    len(blocks) + 1,
                )
            ),
        )

    if mode == "3":
        count = len(blocks)

        if count > 100:
            print(
                f"[WARNING] This will create {count} individual molecule "
                "targets for EACH selected atom/bead."
            )

            if not yes_no(
                "Continue with all individual instances?",
                False,
            ):
                return prompt_specific_site_scope_for_msd(
                    gro,
                    resname,
                )

        return (
            "all_individual",
            list(
                range(
                    1,
                    len(blocks) + 1,
                )
            ),
        )

    _show_instance_preview(
        blocks,
        resname,
    )

    while True:
        raw = ask(
            "Select individual occurrence number(s), e.g. 1 2 5"
        )

        numbers, invalid = _parse_numbered_selection(
            raw,
            len(blocks),
        )

        if invalid or not numbers:
            if invalid:
                print(
                    "[WARNING] Invalid occurrence number(s): "
                    + ", ".join(invalid)
                )

            print(
                f"[ACTION] Enter one or more numbers from "
                f"1-{len(blocks)}."
            )
            continue

        print("\n[SELECTED MOLECULE/RESIDUE INSTANCES]")
        for number in numbers:
            block = blocks[number - 1]

            print(
                f"  [{number:>5d}] {resname}; "
                f"GRO resnr={block['resnr']}; "
                f"atoms/beads={len(block['atoms'])}"
            )

        if yes_no(
            "Confirm these individual molecule/residue instances?",
            True,
        ):
            return (
                "selected_individual",
                numbers,
            )


def _specific_site_index_in_block(
    block: dict,
    atomname: str,
    resname: str,
    occurrence_number: int,
) -> int:
    """
    Return exactly one global atom index for one atom/bead name in one
    individual residue/molecule block.
    """
    matches = [
        idx
        for idx, name in block["atoms"]
        if name == atomname
    ]

    if len(matches) != 1:
        raise RecoverableAnalysisError(
            f"{resname} occurrence #{occurrence_number}: atom/bead "
            f"'{atomname}' matched {len(matches)} time(s); exactly one "
            "is required."
        )

    return matches[0]


def build_specific_site_msd_targets(
    gro: Path,
    resname: str,
    atom_names: list[str],
    scope_mode: str,
    occurrence_numbers: list[int],
) -> list[dict]:
    """
    Construct clean independent MSD targets for the new specific-site mode.
    """
    blocks = lipid_residue_blocks(
        gro,
        resname,
    )

    targets = []

    if scope_mode == "all_combined":
        for atomname in atom_names:
            indices = []

            for number, block in enumerate(
                blocks,
                start=1,
            ):
                indices.append(
                    _specific_site_index_in_block(
                        block,
                        atomname,
                        resname,
                        number,
                    )
                )

            safe_res = sanitize_filename(
                resname
            )
            safe_atom = sanitize_filename(
                atomname
            )

            targets.append(
                {
                    "label": (
                        f"{resname}:{atomname} — "
                        f"all {len(blocks)} molecules combined"
                    ),
                    "group": (
                        f"MSD_SITE_{safe_res}_{safe_atom}_ALL"
                    ),
                    "indices": indices,
                    "resname": resname,
                    "atomname": atomname,
                    "scope": "all molecules combined",
                    "occurrence": "ALL",
                    "gro_resnr": "ALL",
                }
            )

        return targets

    for occurrence_number in occurrence_numbers:
        block = blocks[
            occurrence_number - 1
        ]

        for atomname in atom_names:
            atom_index = _specific_site_index_in_block(
                block,
                atomname,
                resname,
                occurrence_number,
            )

            safe_res = sanitize_filename(
                resname
            )
            safe_atom = sanitize_filename(
                atomname
            )

            targets.append(
                {
                    "label": (
                        f"{resname}:{atomname} — molecule "
                        f"#{occurrence_number} / GRO resnr {block['resnr']}"
                    ),
                    "group": (
                        f"MSD_SITE_{safe_res}_{safe_atom}_"
                        f"MOL{occurrence_number:05d}"
                    ),
                    "indices": [
                        atom_index
                    ],
                    "resname": resname,
                    "atomname": atomname,
                    "scope": "individual molecule/residue",
                    "occurrence": occurrence_number,
                    "gro_resnr": block["resnr"],
                }
            )

    return targets


def calculate_specific_atom_bead_msd(
    gro: Path,
    tpr: Path,
    xtc: Path,
    settings: dict,
):
    """
    NEW MSD MODE 3.

    Calculate diffusion coefficients for a specific atom (AA) or bead (CG)
    belonging to a selected residue/molecule type.

    The user can calculate:
      - ensemble-average diffusion of that site across all molecules,
      - selected individual molecules separately,
      - every individual molecule separately.

    Existing MSD modes remain unchanged.
    """
    print("\n" + "=" * 78)
    print(" MSD MODE 3 — SPECIFIC ATOM / BEAD DIFFUSION")
    print("=" * 78)
    print(
        "This mode calculates MSD/D for a chosen atom in an all-atom model "
        "or a chosen bead in a coarse-grained model."
    )
    print()
    print(
        "The residue/molecule type, atom/bead, and optional individual "
        "molecule instances are selected by NUMBER."
    )

    selected = prompt_numbered_residue_selection(
        gro,
        title="SPECIFIC ATOM/BEAD MSD — RESIDUE/MOLECULE TYPE",
        instruction=(
            "Select ONE residue/molecule TYPE that contains the atom/bead "
            "whose diffusion coefficient you want to calculate."
        ),
        allow_multiple=False,
        selection_name="site-MSD residue/molecule type",
        tpr=tpr,
    )

    resname = selected[0]

    atom_names = prompt_specific_atom_bead_names_for_msd(
        gro,
        tpr,
        resname,
    )

    (
        scope_mode,
        occurrence_numbers,
    ) = prompt_specific_site_scope_for_msd(
        gro,
        resname,
    )

    targets = build_specific_site_msd_targets(
        gro,
        resname,
        atom_names,
        scope_mode,
        occurrence_numbers,
    )

    if not targets:
        raise RecoverableAnalysisError(
            "No valid specific atom/bead MSD target was generated."
        )

    if len(targets) > 100:
        print(
            f"[WARNING] This batch contains {len(targets)} independent "
            "specific-site MSD calculations."
        )

        if not yes_no(
            "Continue with this large MSD batch?",
            False,
        ):
            raise RecoverableAnalysisError(
                "Large specific-site MSD batch cancelled by user."
            )

    subdir = (
        MSD_DIR
        / "Specific_Atom_Bead_Diffusion"
        / sanitize_filename(
            resname
        )
    )
    subdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    ndx_file = (
        subdir
        / "MSD_specific_atom_bead_index.ndx"
    )

    groups = [
        (
            target["group"],
            target["indices"],
        )
        for target in targets
    ]

    write_clean_msd_index(
        ndx_file,
        groups,
    )

    definition_file = (
        subdir
        / "MSD_specific_atom_bead_targets.xlsx"
    )

    definition_headers = [
        "resname",
        "atom_or_bead",
        "scope",
        "molecule_occurrence",
        "GRO_resnr",
        "index_group",
        "selected_particle_count",
    ]

    definition_rows = [
        [
            target["resname"],
            target["atomname"],
            target["scope"],
            target["occurrence"],
            target["gro_resnr"],
            target["group"],
            len(
                target["indices"]
            ),
        ]
        for target in targets
    ]

    write_excel_table(
        definition_file,
        "Specific Site Targets",
        definition_headers,
        definition_rows,
    )

    results = []

    print("\n" + "=" * 78)
    print(
        f" STARTING SPECIFIC ATOM/BEAD MSD BATCH: "
        f"{len(targets)} TARGET(S)"
    )
    print("=" * 78)

    for number, target in enumerate(
        targets,
        start=1,
    ):
        safe = sanitize_filename(
            (
                f"{target['resname']}_"
                f"{target['atomname']}_"
                f"{target['occurrence']}"
            )
        )

        print(
            f"\n---------------- SITE MSD "
            f"{number}/{len(targets)}: "
            f"{target['label']} ----------------"
        )

        xvg_file = (
            subdir
            / f"MSD_site_{safe}_curve.xvg"
        )

        result_file = (
            subdir
            / f"MSD_site_{safe}_diffusion_result.txt"
        )

        png_file = (
            subdir
            / f"MSD_site_{safe}_curve.png"
        )

        log_file, diffusion = run_msd_command(
            xtc,
            tpr,
            ndx_file,
            target["group"],
            xvg_file,
            settings,
        )

        write_diffusion_result_file(
            result_file,
            object_type=(
                "specific atom/bead inside selected "
                "residue/molecule"
            ),
            object_name=target["label"],
            region=(
                "whole trajectory; no membrane-region "
                "classification"
            ),
            group_name=target["group"],
            diffusion=diffusion,
            ndx_file=ndx_file,
            xvg_file=xvg_file,
            log_file=log_file,
            settings=settings,
        )

        append_global_msd_summary(
            "specific_atom_bead",
            target["label"],
            "whole_trajectory",
            target["group"],
            diffusion,
            ndx_file,
            xvg_file,
            result_file,
            settings,
        )

        make_msd_plot(
            xvg_file,
            png_file,
            f"MSD — {target['label']}",
        )

        results.append(
            {
                **target,
                "D": diffusion["D"],
                "error": diffusion["error"],
                "unit": diffusion["unit"],
                "curve": xvg_file,
                "result": result_file,
                "log": log_file,
            }
        )

        print(
            f"[BATCH] Completed {number}/{len(targets)} "
            "specific atom/bead MSD target(s)."
        )

    summary_file = (
        subdir
        / "MSD_specific_atom_bead_diffusion_summary.xlsx"
    )

    summary_headers = [
        "resname",
        "atom_or_bead",
        "scope",
        "molecule_occurrence",
        "GRO_resnr",
        "D_1e-5_cm2_s",
        "error_1e-5_cm2_s",
        "dimension",
        "MSD_curve",
        "result_file",
        "log_file",
    ]

    summary_rows = [
        [
            row["resname"],
            row["atomname"],
            row["scope"],
            row["occurrence"],
            row["gro_resnr"],
            (
                row["D"]
                if row["D"] is not None
                else "NA"
            ),
            (
                row["error"]
                if row["error"] is not None
                else "NA"
            ),
            settings["dimension_label"],
            str(
                row["curve"]
            ),
            str(
                row["result"]
            ),
            str(
                row["log"]
            ),
        ]
        for row in results
    ]

    write_excel_table(
        summary_file,
        "Specific Site Diffusion",
        summary_headers,
        summary_rows,
    )

    print("\n" + "=" * 78)
    print(" SPECIFIC ATOM/BEAD MSD / DIFFUSION SUMMARY")
    print("=" * 78)

    for row in results:
        if row["D"] is not None:
            print(
                f" {row['label']}: "
                f"D = {row['D']} "
                f"(+/- {row['error']}) "
                f"({row['unit']})"
            )
        else:
            print(
                f" {row['label']}: "
                "D = NOT AUTOMATICALLY EXTRACTED"
            )

    print()
    print(
        f"  Target definitions : {definition_file}"
    )
    print(
        f"  MSD index          : {ndx_file}"
    )
    print(
        f"  Summary workbook   : {summary_file}"
    )

    return results


def calculate_selected_molecule_dynamic_region_msd(
    gro: Path,
    tpr: Path,
    xtc: Path,
):
    """
    NEW MSD MODE 4.

    Reuse the existing, detailed dynamic-region engine without deleting or
    replacing Main Menu module 7.

    The existing engine is already general in practice: its target selector
    lists every residue/molecule type from the GRO, so the selected target
    does not have to be a pharmaceutical "drug". It can be any molecule or
    residue that should be tracked through membrane regions.

    The method:
      - identifies membrane HEADGROUP/BARRIER/TAIL boundaries from density;
      - tracks each selected molecule/residue occurrence frame by frame;
      - applies hysteresis;
      - detects continuous residence segments;
      - calculates region-specific MSD/D only within continuous residence
        segments, so separate visits are never concatenated.
    """
    print("\n" + "=" * 78)
    print(
        " MSD MODE 4 — SELECTED MOLECULE / RESIDUE "
        "IN DYNAMIC MEMBRANE REGIONS"
    )
    print("=" * 78)
    print(
        "This mode reuses the complete dynamic-region engine already present "
        "in the original script."
    )
    print()
    print(
        "You may select ANY residue/molecule type from the GRO. "
        "The target is not restricted to a drug."
    )
    print()
    print(
        "For each selected individual target molecule, the trajectory is "
        "classified dynamically into HEADGROUP / BARRIER / TAIL / WATER."
    )
    print()
    print(
        "Region-specific diffusion coefficients are then calculated from "
        "continuous residence segments only."
    )

    settings = prompt_dynamic_region_settings()

    analyze_dynamic_drug_region_diffusion(
        gro,
        tpr,
        xtc,
        settings,
    )

    # Add a lightweight cross-reference inside the standard MSD folder so a
    # user entering through MSD mode 4 can easily locate the existing advanced
    # dynamic-region output directory without duplicating hundreds of files.
    pointer_dir = (
        MSD_DIR
        / "Selected_Molecule_Dynamic_Region_Diffusion"
    )
    pointer_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    pointer_file = (
        pointer_dir
        / "OUTPUT_LOCATION.txt"
    )

    pointer_file.write_text(
        "This MSD mode uses the script's original dynamic-region diffusion "
        "engine without deleting or duplicating it.\n\n"
        f"Full dynamic-region results are stored in:\n"
        f"  {DRUG_REGION_DIR}\n\n"
        "The target selector accepts any residue/molecule type from the GRO; "
        "it is not restricted to a pharmaceutical drug.\n",
        encoding="utf-8",
    )

    print(
        f"[INFO] MSD mode-4 output pointer: {pointer_file}"
    )




def run_msd_menu(state: dict) -> bool:
    """
    Main interactive controller for MSD / diffusion coefficient.

    ORIGINAL FUNCTIONS ARE PRESERVED:
      [1] selected residue/molecule type(s)
      [2] whole membrane divided into structural regions

    NEW FUNCTIONS ARE ADDED:
      [3] specific atom (AA) / bead (CG) inside selected molecule/residue
      [4] selected molecule/residue dynamically classified into membrane
          HEADGROUP / BARRIER / TAIL regions

    Main-menu module 7 remains available independently.
    """
    print("\n======================================================================")
    print(" MSD / DIFFUSION COEFFICIENT CALCULATION")
    print("======================================================================")
    print(
        "GROMACS method: gmx msd calculates mean-square displacement and "
        "obtains D from linear regression of MSD(t)."
    )
    print()
    print("The original MSD functions are retained, with two additional modes.")
    print()
    print(" [1] One or more selected residue/molecule TYPES")
    print(
        "     ORIGINAL MODE. Example: CHOL, HSPC, CER. "
        "Each selected type receives its own MSD curve and D value."
    )
    print()
    print(" [2] The WHOLE MEMBRANE divided into structural regions")
    print(
        "     ORIGINAL MODE. Combine selected membrane lipids and calculate "
        "HEADGROUP / BARRIER / TAIL MSD and D."
    )
    print()
    print(" [3] A SPECIFIC ATOM (AA) or BEAD (CG) inside a molecule/residue")
    print(
        "     NEW. Select residue type -> select atom/bead by NUMBER -> "
        "analyze all molecules together or selected individual molecules."
    )
    print()
    print(
        " [4] A selected MOLECULE / RESIDUE in different dynamic "
        "membrane regions"
    )
    print(
        "     NEW. Track selected molecules through HEADGROUP / BARRIER / "
        "TAIL / WATER and calculate D from continuous residence segments."
    )

    files = resolve_analysis_files(
        [
            ("gro", "gro", "GRO file"),
            ("tpr", "tpr", "TPR file"),
            ("xtc", "xtc", "XTC trajectory"),
        ],
        state,
    )

    gro = files["gro"]
    tpr = files["tpr"]
    xtc = files["xtc"]

    while True:
        print("\n" + "-" * 78)
        print(" MSD CALCULATION TARGET MODES")
        print("-" * 78)
        print("  [1] Selected residue/molecule type(s)")
        print("  [2] Whole membrane HEADGROUP/BARRIER/TAIL")
        print("  [3] Specific atom/bead inside selected molecule/residue")
        print("  [4] Selected molecule/residue in dynamic membrane regions")
        print("-" * 78)

        mode = ask(
            "Choose MSD mode",
            "1",
        )

        if mode not in {
            "1",
            "2",
            "3",
            "4",
        }:
            print(
                "[WARNING] Please choose 1, 2, 3, or 4."
            )
            continue

        # Modes 1/2/3 use the standard gmx msd fitting settings.
        if mode in {
            "1",
            "2",
            "3",
        }:
            settings = prompt_msd_settings()

            if mode == "1":
                calculate_single_species_msd(
                    gro,
                    tpr,
                    xtc,
                    settings,
                )

            elif mode == "2":
                calculate_membrane_region_msd(
                    gro,
                    tpr,
                    xtc,
                    settings,
                )

            else:
                calculate_specific_atom_bead_msd(
                    gro,
                    tpr,
                    xtc,
                    settings,
                )

        # Mode 4 uses the existing continuous-residence dynamic-region engine,
        # which has its own scientifically relevant settings.
        else:
            calculate_selected_molecule_dynamic_region_msd(
                gro,
                tpr,
                xtc,
            )

        action = post_calculation_menu(
            "MSD / DIFFUSION",
            "Calculate another MSD/D mode or target",
        )

        if action == "menu":
            return False

        if action == "exit":
            return True

        reuse = yes_no(
            "Reuse current GRO/TPR/XTC files?",
            True,
        )

        if reuse:
            continue

        files = resolve_analysis_files(
            [
                ("gro", "gro", "GRO file"),
                ("tpr", "tpr", "TPR file"),
                ("xtc", "xtc", "XTC trajectory"),
            ],
            state,
        )

        gro = files["gro"]
        tpr = files["tpr"]
        xtc = files["xtc"]


# ======================================================================
# 6) MOLECULAR / RESIDUE DENSITY DISTRIBUTION — gmx density
# ======================================================================

def prompt_density_target_mode() -> str:
    """
    Clarify exactly WHAT is being used as the density-distribution target.
    """
    print("\n" + "=" * 78)
    print(" DENSITY DISTRIBUTION — CALCULATION TARGET")
    print("=" * 78)
    print("Choose what particles should contribute to the density profile:")
    print()
    print(" [1] Whole residue/molecule TYPE(S) selected by resname")
    print("     Example: CHOL, HSPC, SOL, NA, LIG.")
    print("     All atoms/beads belonging to each selected resname are included.")
    print("     Multiple residue/molecule types can be selected in one batch.")
    print()
    print(" [2] Specific atom/bead subset(s) inside selected resname(s)")
    print("     Example: CHOL:ROH or HSPC:PO4,NC3.")
    print("     Useful for headgroup, marker-atom, ion, or bead density profiles.")
    print()
    print(" [3] Specific individual residue/molecule INSTANCE(S)")
    print("     First select one resname, then select occurrence numbers such as")
    print("     molecule 1, molecule 2, molecule 5, etc.")
    print("     Each selected individual instance receives its own density curve.")
    print("=" * 78)

    while True:
        choice = ask(
            "Choose density target mode: "
            "1 = whole resname type(s); "
            "2 = atom/bead subset(s); "
            "3 = individual molecule/residue instance(s)",
            "1",
        )

        if choice in {"1", "2", "3"}:
            return choice

        print("[WARNING] Please choose 1, 2, or 3.")


def _density_safe_group_name(label: str) -> str:
    clean = sanitize_filename(label)
    clean = clean[:60]
    return f"DENS_{clean}"


def density_target_whole_resnames(
    gro: Path,
    tpr: Path,
) -> list[dict]:
    """
    One or more complete residue/molecule TYPES.
    """
    resnames = prompt_numbered_residue_selection(
        gro,
        title="DENSITY — WHOLE RESIDUE/MOLECULE TYPE SELECTION",
        instruction=(
            "Select one or more residue/molecule TYPES. "
            "Each selected resname will be calculated separately."
        ),
        allow_multiple=True,
        selection_name="density target type",
        tpr=tpr,
    )

    targets = []

    for resname in resnames:
        indices = atoms_for_resname(gro, resname)

        if not indices:
            print(
                f"[WARNING] No atoms/beads were found for resname '{resname}'."
            )
            continue

        targets.append(
            {
                "label": resname,
                "group": _density_safe_group_name(resname),
                "indices": indices,
                "mode": "whole_resname_type",
                "resname": resname,
                "selection": f"all atoms/beads of resname {resname}",
            }
        )

    if not targets:
        raise RecoverableAnalysisError(
            "No valid whole-resname density target was created."
        )

    return targets


def show_density_atom_catalog(
    gro: Path,
    selected_resnames: list[str],
):
    """
    Show the atom/bead names of every selected residue/molecule type ONCE.
    """
    counts = residue_counts_from_gro(gro)

    print("\n" + "=" * 78)
    print(" DENSITY — SELECTED RESNAME ATOM/BEAD CATALOG")
    print("=" * 78)

    for number, resname in enumerate(selected_resnames, start=1):
        blocks = lipid_residue_blocks(gro, resname)
        if not blocks:
            continue

        names = []
        seen = set()
        for _, atomname in blocks[0]["atoms"]:
            if atomname not in seen:
                seen.add(atomname)
                names.append(atomname)

        print()
        print(
            f" [{number}] {resname} "
            f"({counts[resname]} residue/molecule blocks)"
        )

        width = 10
        for i in range(0, len(names), width):
            print(
                "      "
                + "  ".join(
                    f"{name:>6s}"
                    for name in names[i:i + width]
                )
            )

    print("=" * 78)


def prompt_density_atom_mapping(
    gro: Path,
    tpr: Path,
    selected_resnames: list[str],
) -> dict[str, list[str]]:
    """
    Number-based atom/bead subsets for density calculations.

    Every selected resname is handled independently and displayed as:
        [1] CHOL-ROH [2] CHOL-C1 ...
    """
    definitions = {}

    print("\n" + "=" * 78)
    print(" DENSITY — NUMBERED ATOM / BEAD SUBSET SELECTION")
    print("=" * 78)

    for number, resname in enumerate(
        selected_resnames,
        start=1,
    ):
        print(
            f"\nTarget {number}/{len(selected_resnames)}: {resname}"
        )

        definitions[resname] = prompt_numbered_atom_bead_selection(
            gro,
            resname,
            title=f"DENSITY ATOMS/BEADS — {resname}",
            instruction=(
                "Select the atom/bead subset that should contribute to the "
                "density profile."
            ),
            allow_multiple=True,
            min_count=1,
            tpr=tpr,
            require_once_per_residue=True,
        )

    print("\n[DENSITY ATOM/BEAD SUBSETS]")
    for number, resname in enumerate(
        selected_resnames,
        start=1,
    ):
        print(
            f"  [{number}] {resname:<14s}: "
            + ", ".join(
                f"{resname}-{name}"
                for name in definitions[resname]
            )
        )

    return definitions



def density_target_atom_subsets(
    gro: Path,
    tpr: Path,
) -> list[dict]:
    """
    Build one density target per selected residue/molecule type, restricted
    to user-selected atom/bead names.
    """
    resnames = prompt_numbered_residue_selection(
        gro,
        title="DENSITY — RESNAME SELECTION FOR ATOM/BEAD SUBSETS",
        instruction=(
            "First select one or more residue/molecule TYPES. "
            "The next step will ask which atom/bead names inside each type "
            "should contribute to the density profile."
        ),
        allow_multiple=True,
        selection_name="resname",
        tpr=tpr,
    )

    definitions = prompt_density_atom_mapping(
        gro,
        tpr,
        resnames,
    )

    targets = []

    for resname in resnames:
        atom_names = definitions[resname]
        indices = validate_region_atom_names(
            gro,
            resname,
            atom_names,
            "DENSITY",
        )

        atom_tag = "_".join(atom_names)
        label = f"{resname}_{atom_tag}"

        targets.append(
            {
                "label": label,
                "group": _density_safe_group_name(label),
                "indices": indices,
                "mode": "atom_bead_subset",
                "resname": resname,
                "selection": (
                    f"resname {resname}; atom/bead names: "
                    + ", ".join(atom_names)
                ),
            }
        )

    return targets


def _show_instance_preview(
    blocks: list[dict],
    resname: str,
):
    """
    Avoid flooding the terminal when a residue type has thousands of
    instances. Show all if small, otherwise first 25 and last 5.
    """
    n = len(blocks)

    print("\n" + "-" * 78)
    print(
        f" INDIVIDUAL INSTANCE SELECTION — {resname} "
        f"({n} instances)"
    )
    print("-" * 78)
    print(
        "Occurrence No. is the sequential occurrence of this resname "
        "inside the GRO file."
    )
    print()

    if n <= 40:
        positions = list(range(1, n + 1))
    else:
        positions = list(range(1, 26)) + list(range(n - 4, n + 1))

    previous = None
    for pos in positions:
        if previous is not None and pos > previous + 1:
            print("       ...")
        block = blocks[pos - 1]
        print(
            f"  [{pos:>5d}] GRO resnr={block['resnr']:<8s} "
            f"atoms/beads={len(block['atoms'])}"
        )
        previous = pos

    print()
    print(f"Valid occurrence-number range: 1-{n}")
    print("-" * 78)


def density_target_individual_instances(
    gro: Path,
    tpr: Path,
) -> list[dict]:
    """
    Select one resname, then one or multiple individual residue/molecule
    occurrences of that resname.
    """
    selected = prompt_numbered_residue_selection(
        gro,
        title="DENSITY — INDIVIDUAL MOLECULE/RESIDUE TYPE",
        instruction=(
            "Select ONE residue/molecule TYPE first. "
            "You will then choose specific individual occurrences "
            "of this type."
        ),
        allow_multiple=False,
        selection_name="residue/molecule type",
        tpr=tpr,
    )

    resname = selected[0]
    blocks = lipid_residue_blocks(gro, resname)

    if not blocks:
        raise RecoverableAnalysisError(
            f"No residue/molecule blocks named '{resname}' were found."
        )

    _show_instance_preview(
        blocks,
        resname,
    )

    while True:
        raw = ask(
            "Select individual occurrence number(s), e.g. 1 2 5"
        )

        numbers, invalid = _parse_numbered_selection(
            raw,
            len(blocks),
        )

        if invalid or not numbers:
            if invalid:
                print(
                    "[WARNING] Invalid occurrence selection(s): "
                    + ", ".join(invalid)
                )
            print(
                f"[ACTION] Enter one or more numbers from 1-{len(blocks)}."
            )
            continue

        print("\n[SELECTED INDIVIDUAL INSTANCES]")
        for number in numbers:
            block = blocks[number - 1]
            print(
                f"  [{number}] {resname}; "
                f"GRO resnr={block['resnr']}; "
                f"atoms/beads={len(block['atoms'])}"
            )

        if yes_no(
            "Confirm these individual density targets?",
            True,
        ):
            break

    targets = []

    for number in numbers:
        block = blocks[number - 1]
        indices = [idx for idx, _ in block["atoms"]]
        label = (
            f"{resname}_instance{number}_resnr{block['resnr']}"
        )

        targets.append(
            {
                "label": label,
                "group": _density_safe_group_name(label),
                "indices": indices,
                "mode": "individual_instance",
                "resname": resname,
                "selection": (
                    f"{resname} occurrence {number}; "
                    f"GRO resnr {block['resnr']}; "
                    "whole individual residue/molecule block"
                ),
            }
        )

    return targets


def prompt_density_targets(
    gro: Path,
    tpr: Path,
) -> tuple[str, list[dict]]:
    """
    Return density target mode and one or more calculation targets.
    """
    mode = prompt_density_target_mode()

    if mode == "1":
        return mode, density_target_whole_resnames(
            gro,
            tpr,
        )

    if mode == "2":
        return mode, density_target_atom_subsets(
            gro,
            tpr,
        )

    return mode, density_target_individual_instances(
        gro,
        tpr,
    )


def prompt_density_axis() -> tuple[str, list[str]]:
    """
    GROMACS default is Z. If Z is selected, do not add -d.
    """
    print("\n---------------------- DENSITY PROFILE DIRECTION -----------------------")
    print("GROMACS default direction is Z.")
    print("For a planar membrane whose normal is Z, keep the default.")
    print(" [1] Z direction — GROMACS default")
    print(" [2] X direction")
    print(" [3] Y direction")

    while True:
        choice = ask(
            "Choose density-profile direction",
            "1",
        )

        if choice == "1":
            return "Z", []
        if choice == "2":
            return "X", ["-d", "X"]
        if choice == "3":
            return "Y", ["-d", "Y"]

        print("[WARNING] Please choose 1, 2, or 3.")


def prompt_density_slices() -> tuple[int, list[str]]:
    """
    GROMACS default is 50 slices. If default is accepted, omit -sl.
    """
    print("\n------------------------- DENSITY SLICES -------------------------------")
    print("GROMACS default: 50 slices.")

    if yes_no(
        "Use the GROMACS default number of slices (50)?",
        True,
    ):
        return 50, []

    while True:
        raw = ask(
            "Enter number of density slices (-sl)"
        )
        try:
            value = int(raw)
            if value <= 0:
                raise ValueError
            return value, ["-sl", str(value)]
        except ValueError:
            print("[WARNING] Please enter a positive integer.")


def prompt_density_options() -> dict:
    """
    User-requested density options:
      - optional number density
      - optional center
      - optional symm
    Defaults are implemented by NOT adding command flags.
    """
    begin_ps, end_ps, dt_ps = prompt_time_window(
        require_dt=True,
    )

    axis, axis_args = prompt_density_axis()
    slices, slice_args = prompt_density_slices()

    print("\n------------------------- DENSITY TYPE --------------------------------")
    print("GROMACS default is MASS density.")
    print(
        "If you choose No below, '-dens' is NOT added and the "
        "GROMACS default mass density is used."
    )

    number_density = yes_no(
        "Calculate NUMBER density instead of the default MASS density?",
        False,
    )

    if number_density:
        density_type = "number"
        density_args = ["-dens", "number"]
        unit = "1/nm^3 (selected atoms/beads)"
        print(
            "[INFO] Using: -dens number"
        )
        print(
            "[NOTE] Number density counts the selected particles/atoms/beads. "
            "For a molecular number-density profile, selecting ONE "
            "representative atom/bead per molecule is usually the cleanest "
            "definition."
        )
    else:
        density_type = "mass"
        density_args = []
        unit = "kg/m^3"
        print(
            "[INFO] Using GROMACS default MASS density; no -dens flag added."
        )

    print("\n------------------------- CENTER / SYMM -------------------------------")
    print(
        "If No is selected, the corresponding flag is NOT added and "
        "GROMACS default behavior is retained."
    )

    use_center = yes_no(
        "Use -center (bin relative to a changing center reference)?",
        False,
    )

    use_symm = yes_no(
        "Use -symm (symmetrize profile around the center)?",
        False,
    )

    if use_symm:
        print(
            "[INFO] GROMACS -symm automatically turns on centering."
        )

    return {
        "begin_ps": begin_ps,
        "end_ps": end_ps,
        "dt_ps": dt_ps,
        "axis": axis,
        "axis_args": axis_args,
        "slices": slices,
        "slice_args": slice_args,
        "density_type": density_type,
        "density_args": density_args,
        "unit": unit,
        "use_center": use_center,
        "use_symm": use_symm,
        "center_effective": use_center or use_symm,
    }


def prompt_density_center_reference(
    gro: Path,
    tpr: Path,
    targets: list[dict],
) -> tuple[str | None, list[int], str]:
    """
    If center/symm is active, choose ONE common center reference for all
    target curves so the profiles stay mutually comparable.

    Default: combined atoms/beads of all selected density targets.
    """
    print("\n--------------------- DENSITY CENTER REFERENCE ------------------------")
    print(
        "All density curves in this batch should use the SAME centering "
        "reference for direct comparison."
    )

    if yes_no(
        "Use ALL selected density targets combined as the center reference?",
        True,
    ):
        indices = []
        seen = set()

        for target in targets:
            for idx in target["indices"]:
                if idx not in seen:
                    seen.add(idx)
                    indices.append(idx)

        return (
            "DENSITY_CENTER_ALL_TARGETS",
            indices,
            "all selected density targets combined",
        )

    print()
    print(
        "Select one or more whole residue/molecule TYPES to build the "
        "common center reference."
    )

    resnames = prompt_numbered_residue_selection(
        gro,
        title="DENSITY — CUSTOM CENTER REFERENCE",
        instruction=(
            "Select the residue/molecule type(s) whose combined atoms/beads "
            "define the center. For a membrane-centered profile, select the "
            "membrane components that should define the bilayer center."
        ),
        allow_multiple=True,
        selection_name="center-reference type",
        tpr=tpr,
    )

    indices = []
    seen = set()

    for resname in resnames:
        for idx in atoms_for_resname(
            gro,
            resname,
        ):
            if idx not in seen:
                seen.add(idx)
                indices.append(idx)

    if not indices:
        raise RecoverableAnalysisError(
            "The custom density center reference is empty."
        )

    return (
        "DENSITY_CENTER_REFERENCE",
        indices,
        "combined resnames: " + ", ".join(resnames),
    )


def write_clean_density_index(
    ndx_file: Path,
    targets: list[dict],
    center_group: str | None = None,
    center_indices: list[int] | None = None,
):
    """
    Create a density-specific NDX from scratch.
    """
    ndx_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with ndx_file.open(
        "w",
        encoding="utf-8",
    ) as fh:

        if center_group is not None:
            if not center_indices:
                raise RecoverableAnalysisError(
                    "Density center-reference group is empty."
                )
            write_ndx_group(
                fh,
                center_group,
                center_indices,
            )

        for target in targets:
            if not target["indices"]:
                raise RecoverableAnalysisError(
                    f"Density target '{target['label']}' is empty."
                )

            write_ndx_group(
                fh,
                target["group"],
                target["indices"],
            )

    print(f"[OK] Density index generated: {ndx_file}")

    if center_group is not None:
        print(
            f"     [ {center_group} ] : "
            f"{len(center_indices)} atoms/beads"
        )

    for target in targets:
        print(
            f"     [ {target['group']} ] : "
            f"{len(target['indices'])} atoms/beads"
        )


def run_density_command(
    xtc: Path,
    tpr: Path,
    ndx_file: Path,
    target: dict,
    xvg_file: Path,
    settings: dict,
    center_group: str | None,
) -> Path:
    """
    Run one gmx density profile.

    One target is run at a time to give each target an independent XVG,
    CSV, PNG and clearly named result file.
    """
    cmd = [
        GMX,
        "density",
        "-f",
        str(xtc),
        "-s",
        str(tpr),
        "-n",
        str(ndx_file),
        "-o",
        str(xvg_file),
        *gmx_time_args(
            settings["begin_ps"],
            settings["end_ps"],
            settings["dt_ps"],
        ),
        *settings["axis_args"],
        *settings["slice_args"],
        *settings["density_args"],
    ]

    if settings["use_center"]:
        cmd.append("-center")

    if settings["use_symm"]:
        cmd.append("-symm")

    if settings["center_effective"]:
        if not center_group:
            raise RecoverableAnalysisError(
                "Centering is active but no center-reference group exists."
            )
        input_text = (
            center_group
            + "\n"
            + target["group"]
            + "\n"
        )
    else:
        input_text = target["group"] + "\n"

    log_file = timestamped_log_path(xvg_file.with_suffix(".log"))

    result = _animated_process(
        cmd,
        input_text=input_text,
        label=f"Density: {target['label']}",
    )

    log_file.write_text(
        log_timestamp_header() + "COMMAND\n"
        + " ".join(map(str, cmd))
        + "\n\nCENTER GROUP\n"
        + (
            center_group
            if settings["center_effective"]
            else "GROMACS default / no center flag"
        )
        + "\n\nDENSITY TARGET GROUP\n"
        + target["group"]
        + "\n\nSTDOUT\n"
        + (result.stdout or "")
        + "\nSTDERR\n"
        + (result.stderr or ""),
        encoding="utf-8",
    )

    if result.returncode != 0:
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(
                result.stderr,
                file=sys.stderr,
            )
        raise RecoverableAnalysisError(
            f"gmx density failed for target '{target['label']}'. "
            f"Full command log: {log_file}"
        )

    if (
        not xvg_file.exists()
        or xvg_file.stat().st_size == 0
    ):
        raise RecoverableAnalysisError(
            f"gmx density did not generate: {xvg_file}"
        )

    print(
        f"[OK] Density profile generated: {xvg_file}"
    )
    print(
        f"[INFO] Command log: {log_file}"
    )

    return log_file


def density_profile_statistics(
    xvg_file: Path,
) -> dict:
    """
    Basic profile statistics for summary/reporting.
    """
    rows = parse_xvg(
        xvg_file,
        min_cols=2,
    )

    if not rows:
        raise RecoverableAnalysisError(
            f"No density-profile data could be read from {xvg_file}."
        )

    x = [r[0] for r in rows]
    y = [r[1] for r in rows]

    max_idx = max(
        range(len(y)),
        key=lambda i: y[i],
    )
    min_idx = min(
        range(len(y)),
        key=lambda i: y[i],
    )

    return {
        "points": len(y),
        "coord_min": min(x),
        "coord_max": max(x),
        "density_mean": mean(y),
        "density_min": y[min_idx],
        "density_min_position": x[min_idx],
        "density_max": y[max_idx],
        "density_max_position": x[max_idx],
    }


def write_density_csv(
    xvg_file: Path,
    csv_file: Path,
    axis: str,
):
    rows = parse_xvg(
        xvg_file,
        min_cols=2,
    )

    with csv_file.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                f"{axis}_coordinate_nm",
                "density",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    f"{row[0]:.10f}",
                    f"{row[1]:.10f}",
                ]
            )


def make_density_plot(
    xvg_file: Path,
    png_file: Path,
    title: str,
    axis: str,
    unit: str,
):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "[INFO] matplotlib is not installed; density PNG skipped. "
            "The XVG/CSV files remain available."
        )
        return

    rows = parse_xvg(
        xvg_file,
        min_cols=2,
    )

    if not rows:
        return

    x = [r[0] for r in rows]
    y = [r[1] for r in rows]

    fig, ax = plt.subplots(
        figsize=(7.2, 5.0)
    )
    ax.plot(
        x,
        y,
        linewidth=1.2,
    )
    ax.set_title(title)
    ax.set_xlabel(
        f"{axis} coordinate (nm)"
    )
    ax.set_ylabel(
        f"Density ({unit})"
    )
    fig.tight_layout()
    fig.savefig(
        png_file,
        dpi=300,
    )
    plt.close(fig)

    print(
        f"[OK] Density plot generated: {png_file}"
    )


def write_density_result_file(
    result_file: Path,
    target: dict,
    stats: dict,
    settings: dict,
    center_description: str,
    ndx_file: Path,
    xvg_file: Path,
    csv_file: Path,
    log_file: Path,
):
    """
    Write a human-readable density-profile result summary.
    """
    lines = [
        "=" * 78,
        "GROMACS DENSITY DISTRIBUTION RESULT",
        "=" * 78,
        f"Calculation target mode : {target['mode']}",
        f"Calculation target      : {target['label']}",
        f"Source resname          : {target['resname']}",
        f"Selection definition    : {target['selection']}",
        f"Selected atoms/beads    : {len(target['indices'])}",
        "",
        f"Density type            : {settings['density_type']}",
        f"Density unit            : {settings['unit']}",
        f"Profile direction       : {settings['axis']}",
        f"Number of slices        : {settings['slices']}",
        f"Use -center             : {settings['use_center']}",
        f"Use -symm               : {settings['use_symm']}",
        f"Effective centering     : {settings['center_effective']}",
        f"Center reference        : {center_description}",
        "",
        f"Trajectory begin        : {settings['begin_ps']} ps",
        f"Trajectory end          : "
        + (
            f"{settings['end_ps']} ps"
            if settings["end_ps"] is not None
            else "trajectory end"
        ),
        f"Sampling interval       : {settings['dt_ps']} ps",
        "",
        f"Profile data points     : {stats['points']}",
        f"Coordinate range        : "
        f"{stats['coord_min']:.8f} to {stats['coord_max']:.8f} nm",
        f"Mean profile density    : "
        f"{stats['density_mean']:.8f} {settings['unit']}",
        f"Minimum density         : "
        f"{stats['density_min']:.8f} {settings['unit']} "
        f"at {stats['density_min_position']:.8f} nm",
        f"Maximum density         : "
        f"{stats['density_max']:.8f} {settings['unit']} "
        f"at {stats['density_max_position']:.8f} nm",
        "",
        f"Index file              : {ndx_file}",
        f"Density XVG             : {xvg_file}",
        f"Density CSV             : {csv_file}",
        f"Command log             : {log_file}",
        "=" * 78,
        "",
    ]

    result_file.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def write_density_target_definition_excel(
    path: Path,
    targets: list[dict],
    center_description: str,
):
    headers = [
        "target_label",
        "target_mode",
        "source_resname",
        "selection_definition",
        "selected_atom_bead_count",
        "index_group",
        "center_reference",
    ]

    rows = [
        [
            target["label"],
            target["mode"],
            target["resname"],
            target["selection"],
            len(target["indices"]),
            target["group"],
            center_description,
        ]
        for target in targets
    ]

    write_excel_table(
        path,
        "Density Targets",
        headers,
        rows,
    )


def append_density_summary_excel(
    target: dict,
    stats: dict,
    settings: dict,
    center_description: str,
    ndx_file: Path,
    xvg_file: Path,
    csv_file: Path,
    result_file: Path,
    log_file: Path,
):
    """
    Persistent Excel summary across density calculations.
    """
    summary = (
        DENSITY_DIR
        / "Density_distribution_summary.xlsx"
    )
    cache = (
        DENSITY_DIR
        / ".Density_distribution_summary_data.json"
    )

    headers = [
        "target_mode",
        "target",
        "source_resname",
        "selection_definition",
        "selected_atom_bead_count",
        "density_type",
        "density_unit",
        "profile_axis",
        "slices",
        "use_center",
        "use_symm",
        "center_reference",
        "mean_profile_density",
        "minimum_density",
        "minimum_position_nm",
        "maximum_density",
        "maximum_position_nm",
        "index_file",
        "density_xvg",
        "density_csv",
        "result_file",
        "log_file",
    ]

    row = [
        target["mode"],
        target["label"],
        target["resname"],
        target["selection"],
        len(target["indices"]),
        settings["density_type"],
        settings["unit"],
        settings["axis"],
        settings["slices"],
        str(settings["use_center"]),
        str(settings["use_symm"]),
        center_description,
        stats["density_mean"],
        stats["density_min"],
        stats["density_min_position"],
        stats["density_max"],
        stats["density_max_position"],
        str(ndx_file),
        str(xvg_file),
        str(csv_file),
        str(result_file),
        str(log_file),
    ]

    append_excel_summary(
        summary,
        cache,
        "Density Summary",
        headers,
        row,
    )


def run_density_batch(
    gro: Path,
    tpr: Path,
    xtc: Path,
    mode: str,
    targets: list[dict],
    settings: dict,
):
    """
    Create the index once and automatically calculate every requested
    density target.
    """
    mode_names = {
        "1": "Whole_Resname_Types",
        "2": "Atom_Bead_Subsets",
        "3": "Individual_Instances",
    }

    subdir = (
        DENSITY_DIR
        / mode_names.get(
            mode,
            "Density_Targets",
        )
    )
    subdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if settings["center_effective"]:
        (
            center_group,
            center_indices,
            center_description,
        ) = prompt_density_center_reference(
            gro,
            tpr,
            targets,
        )
    else:
        center_group = None
        center_indices = []
        center_description = (
            "GROMACS default; no -center/-symm flag added"
        )

    ndx_file = (
        subdir
        / "Density_targets_index.ndx"
    )

    write_clean_density_index(
        ndx_file,
        targets,
        center_group=center_group,
        center_indices=center_indices,
    )

    definition_file = (
        subdir
        / "Density_target_definitions.xlsx"
    )

    write_density_target_definition_excel(
        definition_file,
        targets,
        center_description,
    )

    print("\n" + "=" * 78)
    print(
        f" STARTING DENSITY BATCH: {len(targets)} TARGET(S)"
    )
    print("=" * 78)

    results = []

    density_tag = settings["density_type"]
    axis_tag = settings["axis"]

    for number, target in enumerate(
        targets,
        start=1,
    ):
        safe = sanitize_filename(
            target["label"]
        )

        print(
            f"\n---------------- DENSITY {number}/{len(targets)}: "
            f"{target['label']} ----------------"
        )

        base = (
            f"Density_{safe}_{density_tag}_{axis_tag}"
        )

        xvg_file = subdir / f"{base}.xvg"
        csv_file = subdir / f"{base}.csv"
        png_file = subdir / f"{base}.png"
        result_file = subdir / f"{base}_result.txt"

        log_file = run_density_command(
            xtc,
            tpr,
            ndx_file,
            target,
            xvg_file,
            settings,
            center_group,
        )

        write_density_csv(
            xvg_file,
            csv_file,
            settings["axis"],
        )

        stats = density_profile_statistics(
            xvg_file,
        )

        write_density_result_file(
            result_file,
            target,
            stats,
            settings,
            center_description,
            ndx_file,
            xvg_file,
            csv_file,
            log_file,
        )

        make_density_plot(
            xvg_file,
            png_file,
            f"Density distribution — {target['label']}",
            settings["axis"],
            settings["unit"],
        )

        append_density_summary_excel(
            target,
            stats,
            settings,
            center_description,
            ndx_file,
            xvg_file,
            csv_file,
            result_file,
            log_file,
        )

        results.append(
            {
                "target": target["label"],
                "xvg": xvg_file,
                "csv": csv_file,
                "png": png_file,
                "result": result_file,
                "log": log_file,
                "stats": stats,
            }
        )

        print(
            f"[BATCH] Completed {number}/{len(targets)} "
            f"density target(s): {target['label']}"
        )

    print("\n" + "=" * 78)
    print(" DENSITY DISTRIBUTION SUMMARY")
    print("=" * 78)

    for result in results:
        stats = result["stats"]
        print(
            f" {result['target']:<24s}: "
            f"max = {stats['density_max']:.6g} "
            f"{settings['unit']} at "
            f"{stats['density_max_position']:.6g} nm"
        )

    print()
    print(f"  Target definitions : {definition_file}")
    print(f"  Density index      : {ndx_file}")
    print(
        f"  Global Excel summary: "
        f"{DENSITY_DIR / 'Density_distribution_summary.xlsx'}"
    )

    return results


def run_density_menu(state: dict) -> bool:
    """
    Interactive density-distribution module.
    """
    print("\n" + "=" * 78)
    print(" MOLECULAR / RESIDUE DENSITY DISTRIBUTION")
    print("=" * 78)
    print(
        "Method: gmx density computes one-dimensional partial-density "
        "profiles across the simulation box."
    )
    print(
        "For planar membrane systems, the usual profile direction is Z."
    )
    print()

    files = resolve_analysis_files(
        [
            ("gro", "gro", "GRO file"),
            ("tpr", "tpr", "TPR file"),
            ("xtc", "xtc", "XTC trajectory"),
        ],
        state,
    )

    gro = files["gro"]
    tpr = files["tpr"]
    xtc = files["xtc"]

    while True:
        mode, targets = prompt_density_targets(
            gro,
            tpr,
        )

        settings = prompt_density_options()

        run_density_batch(
            gro,
            tpr,
            xtc,
            mode,
            targets,
            settings,
        )

        action = post_calculation_menu(
            "DENSITY DISTRIBUTION",
            "Calculate another density distribution "
            "(another target/setting/system)",
        )

        if action == "menu":
            return False

        if action == "exit":
            return True

        reuse = yes_no(
            "Reuse the current GRO/TPR/XTC files for another density calculation?",
            True,
        )

        if reuse:
            continue

        files = resolve_analysis_files(
            [
                ("gro", "gro", "GRO file"),
                ("tpr", "tpr", "TPR file"),
                ("xtc", "xtc", "XTC trajectory"),
            ],
            state,
        )

        gro = files["gro"]
        tpr = files["tpr"]
        xtc = files["xtc"]




# ======================================================================
# 7) DRUG MSD / DIFFUSION BY DYNAMIC MEMBRANE REGION
# ======================================================================

def prompt_dynamic_drug_type(
    gro: Path,
    tpr: Path,
) -> tuple[str, list[dict]]:
    """
    Select ONE drug/guest residue type, then choose all or selected individual
    molecule/residue occurrences. Every occurrence is tracked independently.
    """
    selected = prompt_numbered_residue_selection(
        gro,
        title="DYNAMIC REGION DIFFUSION — DRUG / GUEST TYPE",
        instruction=(
            "Select ONE drug/guest residue or molecule TYPE. "
            "Each individual occurrence of this resname can then be tracked "
            "independently through HEADGROUP / BARRIER / TAIL / WATER regions."
        ),
        allow_multiple=False,
        selection_name="drug/guest type",
        tpr=tpr,
    )

    resname = selected[0]
    blocks = lipid_residue_blocks(
        gro,
        resname,
    )

    if not blocks:
        raise RecoverableAnalysisError(
            f"No molecules/residue blocks named '{resname}' were found."
        )

    print()
    print(
        f"[INFO] Detected {len(blocks)} individual {resname} "
        "molecule/residue occurrence(s)."
    )

    if yes_no(
        f"Analyze ALL {len(blocks)} individual {resname} occurrences?",
        True,
    ):
        chosen_numbers = list(
            range(1, len(blocks) + 1)
        )
    else:
        _show_instance_preview(
            blocks,
            resname,
        )

        while True:
            raw = ask(
                "Select individual occurrence number(s), e.g. 1 2 5"
            )
            chosen_numbers, invalid = _parse_numbered_selection(
                raw,
                len(blocks),
            )

            if invalid or not chosen_numbers:
                if invalid:
                    print(
                        "[WARNING] Invalid occurrence selection(s): "
                        + ", ".join(invalid)
                    )
                print(
                    f"[ACTION] Enter one or more numbers from 1-{len(blocks)}."
                )
                continue

            if yes_no(
                f"Analyze {len(chosen_numbers)} selected individual "
                f"{resname} occurrence(s)?",
                True,
            ):
                break

    targets = []

    for number in chosen_numbers:
        block = blocks[number - 1]
        group = f"DYN_DRUG_{number:03d}"
        label = (
            f"{resname}_{number:03d}_resnr{block['resnr']}"
        )
        targets.append(
            {
                "number": number,
                "label": label,
                "group": group,
                "resname": resname,
                "resnr": block["resnr"],
                "indices": [
                    idx
                    for idx, _ in block["atoms"]
                ],
            }
        )

    print("\n[DRUG / GUEST TARGETS]")
    print(f"  resname             : {resname}")
    print(f"  individual targets  : {len(targets)}")

    if len(targets) <= 20:
        for target in targets:
            print(
                f"  [{target['number']:>3d}] "
                f"{target['label']}"
            )
    else:
        for target in targets[:10]:
            print(
                f"  [{target['number']:>3d}] "
                f"{target['label']}"
            )
        print("       ...")
        for target in targets[-5:]:
            print(
                f"  [{target['number']:>3d}] "
                f"{target['label']}"
            )

    return resname, targets


def prompt_dynamic_region_settings() -> dict:
    """
    Recommended defaults prioritize robust membrane-region classification
    and lateral diffusion statistics while keeping interaction compact.
    """
    print("\n" + "=" * 78)
    print(" DYNAMIC REGION ANALYSIS SETTINGS")
    print("=" * 78)

    begin_ps, end_ps, dt_ps = prompt_time_window(
        require_dt=True,
    )

    print()
    print("Recommended defaults:")
    print("  region-density slices             = 200")
    print("  HEADGROUP outer decay threshold   = 10% of normalized head density")
    print("  boundary hysteresis               = 0.10 nm")
    print("  minimum continuous residence      = 200 ps")
    print("  maximum MSD lag                   = 50% of each residence segment")
    print("  automatic D fit window            = 10%-50% of available lag range")
    print("  diffusion dimension               = XY lateral (recommended)")
    print()

    defaults = yes_no(
        "Use these recommended dynamic-region settings?",
        True,
    )

    if defaults:
        density_slices = 200
        head_decay = 0.10
        hysteresis_nm = 0.10
        min_residence_ps = 200.0
        max_lag_ps = None
        auto_fit = True
        fit_start_ps = None
        fit_end_ps = None
    else:
        while True:
            raw = ask(
                "Density slices used for automatic membrane-region boundaries",
                "200",
            )
            try:
                density_slices = int(raw)
                if density_slices < 20:
                    raise ValueError
                break
            except ValueError:
                print(
                    "[WARNING] Enter an integer >= 20."
                )

        while True:
            head_decay = ask_float(
                "Normalized HEADGROUP outer decay threshold (0-1)",
                0.10,
            )
            if 0 < head_decay < 1:
                break
            print(
                "[WARNING] Enter a value between 0 and 1."
            )

        while True:
            hysteresis_nm = ask_float(
                "Boundary hysteresis tolerance, nm",
                0.10,
            )
            if hysteresis_nm >= 0:
                break
            print(
                "[WARNING] Hysteresis must be >= 0."
            )

        while True:
            min_residence_ps = ask_float(
                "Minimum continuous residence time, ps",
                200.0,
            )
            if min_residence_ps >= 0:
                break
            print(
                "[WARNING] Minimum residence time must be >= 0."
            )

        raw = ask(
            "Maximum MSD lag time, ps; press Enter for automatic "
            "(50% of each valid segment)",
            "",
        )
        if raw.strip():
            try:
                max_lag_ps = float(raw)
                if max_lag_ps <= 0:
                    raise ValueError
            except ValueError:
                print(
                    "[WARNING] Invalid max lag; automatic mode will be used."
                )
                max_lag_ps = None
        else:
            max_lag_ps = None

        auto_fit = yes_no(
            "Use automatic diffusion fit range (10%-50% of available lag)?",
            True,
        )

        if auto_fit:
            fit_start_ps = None
            fit_end_ps = None
        else:
            while True:
                fit_start_ps = ask_float(
                    "Diffusion fit start lag, ps",
                    0.0,
                )
                fit_end_ps = ask_float(
                    "Diffusion fit end lag, ps",
                    1000.0,
                )
                if (
                    fit_start_ps >= 0
                    and fit_end_ps > fit_start_ps
                ):
                    break
                print(
                    "[WARNING] Fit end must be greater than fit start."
                )

    print("\n---------------- DIFFUSION DIMENSION ----------------")
    print(" [1] XY lateral diffusion — recommended for membrane-region comparison")
    print(" [2] 3D diffusion")
    print(" [3] Z-only diffusion")
    print()
    print(
        "XY is recommended because membrane regions are defined by Z depth; "
        "including Z can mix lateral mobility with confinement by region boundaries."
    )

    while True:
        dim_choice = ask(
            "Choose diffusion dimension",
            "1",
        )

        if dim_choice == "1":
            dimension = "XY"
            divisor = 4.0
            break
        if dim_choice == "2":
            dimension = "3D"
            divisor = 6.0
            break
        if dim_choice == "3":
            dimension = "Z"
            divisor = 2.0
            break

        print(
            "[WARNING] Please choose 1, 2, or 3."
        )

    return {
        "begin_ps": begin_ps,
        "end_ps": end_ps,
        "dt_ps": dt_ps,
        "density_slices": density_slices,
        "head_decay": head_decay,
        "hysteresis_nm": hysteresis_nm,
        "min_residence_ps": min_residence_ps,
        "max_lag_ps": max_lag_ps,
        "auto_fit": auto_fit,
        "fit_start_ps": fit_start_ps,
        "fit_end_ps": fit_end_ps,
        "dimension": dimension,
        "diffusion_divisor": divisor,
    }


def dynamic_region_indices(
    gro: Path,
    lipid_names: list[str],
    definitions: dict,
) -> tuple[list[int], dict[str, list[int]]]:
    """
    Build:
      - all membrane atoms/beads for the dynamic membrane center reference
      - combined HEADGROUP/BARRIER/TAIL atom/bead indices
    """
    membrane_indices = []
    seen = set()

    for resname in lipid_names:
        for idx in atoms_for_resname(
            gro,
            resname,
        ):
            if idx not in seen:
                seen.add(idx)
                membrane_indices.append(idx)

    if not membrane_indices:
        raise RecoverableAnalysisError(
            "The selected membrane lipid set produced an empty membrane group."
        )

    region_indices = {
        "HEADGROUP": [],
        "BARRIER": [],
        "TAIL": [],
    }

    for region in (
        "HEADGROUP",
        "BARRIER",
        "TAIL",
    ):
        for resname in lipid_names:
            region_indices[region].extend(
                validate_region_atom_names(
                    gro,
                    resname,
                    definitions[resname][region],
                    region,
                )
            )

        if not region_indices[region]:
            raise RecoverableAnalysisError(
                f"The combined {region} membrane group is empty."
            )

    return membrane_indices, region_indices


def write_dynamic_region_index(
    ndx_file: Path,
    membrane_indices: list[int],
    region_indices: dict[str, list[int]],
    drug_targets: list[dict],
):
    """
    Create one clean NDX used by density, COM and distance calculations.
    """
    ndx_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with ndx_file.open(
        "w",
        encoding="utf-8",
    ) as fh:
        write_ndx_group(
            fh,
            "DYN_MEMBRANE_CENTER",
            membrane_indices,
        )
        write_ndx_group(
            fh,
            "DYN_HEADGROUP",
            region_indices["HEADGROUP"],
        )
        write_ndx_group(
            fh,
            "DYN_BARRIER",
            region_indices["BARRIER"],
        )
        write_ndx_group(
            fh,
            "DYN_TAIL",
            region_indices["TAIL"],
        )

        for target in drug_targets:
            write_ndx_group(
                fh,
                target["group"],
                target["indices"],
            )

    print(
        f"[OK] Dynamic-region index generated: {ndx_file}"
    )


def run_dynamic_boundary_density_profiles(
    xtc: Path,
    tpr: Path,
    ndx_file: Path,
    subdir: Path,
    settings: dict,
    region_indices: dict[str, list[int]],
):
    """
    Generate normalized-region boundary source profiles using:
      gmx density
      -dens number
      -center
      -symm
      -sl <density_slices>

    Each regional curve is normalized before automatic intersection analysis,
    so unequal atom counts between HEADGROUP/BARRIER/TAIL do not dominate the
    boundary merely through amplitude.
    """
    density_settings = {
        "begin_ps": settings["begin_ps"],
        "end_ps": settings["end_ps"],
        "dt_ps": settings["dt_ps"],
        "axis_args": [],
        "slice_args": [
            "-sl",
            str(settings["density_slices"]),
        ],
        "density_args": [
            "-dens",
            "number",
        ],
        "use_center": True,
        "use_symm": True,
        "center_effective": True,
    }

    paths = {}

    for region in (
        "HEADGROUP",
        "BARRIER",
        "TAIL",
    ):
        target = {
            "label": f"DYNAMIC_{region}",
            "group": f"DYN_{region}",
            "indices": region_indices[region],
        }

        out = (
            subdir
            / f"Boundary_density_{region.lower()}.xvg"
        )

        run_density_command(
            xtc,
            tpr,
            ndx_file,
            target,
            out,
            density_settings,
            "DYN_MEMBRANE_CENTER",
        )

        paths[region] = out

    return paths


def _normalized_positive_density_profile(
    xvg_file: Path,
) -> list[tuple[float, float]]:
    """
    Fold a centered/symmetrized Z density profile to one positive membrane
    half and normalize the density to its own maximum.
    """
    rows = parse_xvg(
        xvg_file,
        min_cols=2,
    )

    if len(rows) < 5:
        raise RecoverableAnalysisError(
            f"Density profile contains too few data points: {xvg_file}"
        )

    xs = [row[0] for row in rows]
    ys = [max(0.0, row[1]) for row in rows]

    center = (
        min(xs) + max(xs)
    ) / 2.0

    positive = [
        (
            abs(x - center),
            y,
        )
        for x, y in zip(xs, ys)
        if x >= center
    ]

    positive.sort(
        key=lambda item: item[0]
    )

    ymax = max(
        y for _, y in positive
    )

    if ymax <= 0:
        raise RecoverableAnalysisError(
            f"Density profile has no positive signal: {xvg_file}"
        )

    return [
        (d, y / ymax)
        for d, y in positive
    ]


def _linear_profile_value(
    profile: list[tuple[float, float]],
    x: float,
) -> float:
    """
    Simple linear interpolation without external dependencies.
    """
    if x <= profile[0][0]:
        return profile[0][1]
    if x >= profile[-1][0]:
        return profile[-1][1]

    for i in range(1, len(profile)):
        x1, y1 = profile[i - 1]
        x2, y2 = profile[i]

        if x1 <= x <= x2:
            if x2 == x1:
                return y1
            f = (
                x - x1
            ) / (
                x2 - x1
            )
            return y1 + f * (
                y2 - y1
            )

    return profile[-1][1]


def _profile_peak_position(
    profile: list[tuple[float, float]],
) -> float:
    return max(
        profile,
        key=lambda item: item[1],
    )[0]


def _find_profile_intersection(
    profile_a: list[tuple[float, float]],
    profile_b: list[tuple[float, float]],
    peak_a: float,
    peak_b: float,
) -> tuple[float, str]:
    """
    Find normalized-profile crossing between two regional peaks.
    If no strict sign change exists, use the minimum-difference point.
    """
    lo = min(
        peak_a,
        peak_b,
    )
    hi = max(
        peak_a,
        peak_b,
    )

    grid = sorted(
        {
            x
            for x, _ in profile_a + profile_b
            if lo <= x <= hi
        }
    )

    if len(grid) < 2:
        return (
            (lo + hi) / 2.0,
            "peak_midpoint_fallback",
        )

    diffs = [
        _linear_profile_value(
            profile_a,
            x,
        )
        - _linear_profile_value(
            profile_b,
            x,
        )
        for x in grid
    ]

    roots = []

    for i in range(1, len(grid)):
        x1 = grid[i - 1]
        x2 = grid[i]
        y1 = diffs[i - 1]
        y2 = diffs[i]

        if y1 == 0:
            roots.append(x1)
            continue

        if y1 * y2 < 0:
            root = (
                x1
                - y1 * (
                    x2 - x1
                ) / (
                    y2 - y1
                )
            )
            roots.append(root)

    midpoint = (
        lo + hi
    ) / 2.0

    if roots:
        chosen = min(
            roots,
            key=lambda x: abs(
                x - midpoint
            ),
        )
        return chosen, "normalized_density_intersection"

    chosen = min(
        grid,
        key=lambda x: abs(
            _linear_profile_value(
                profile_a,
                x,
            )
            - _linear_profile_value(
                profile_b,
                x,
            )
        ),
    )

    return chosen, "minimum_normalized_difference"


def _find_head_outer_boundary(
    head_profile: list[tuple[float, float]],
    head_peak: float,
    threshold: float,
) -> tuple[float, str]:
    """
    Outer HEADGROUP/WATER boundary:
    first post-peak position where normalized head density stays below
    threshold for at least 3 consecutive bins.
    """
    post = [
        item
        for item in head_profile
        if item[0] >= head_peak
    ]

    if not post:
        return (
            head_profile[-1][0],
            "profile_end_fallback",
        )

    for i in range(
        0,
        max(0, len(post) - 2),
    ):
        window = post[i:i + 3]

        if (
            len(window) == 3
            and all(
                y <= threshold
                for _, y in window
            )
        ):
            return (
                window[0][0],
                f"head_decay_{threshold:.3f}",
            )

    return (
        post[-1][0],
        "profile_end_fallback",
    )


def write_dynamic_boundary_profile_csv(
    out_file: Path,
    profiles: dict,
):
    """
    Export one normalized positive-side density table for auditing.
    """
    grid = sorted(
        {
            x
            for profile in profiles.values()
            for x, _ in profile
        }
    )

    with out_file.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "distance_from_membrane_center_nm",
                "HEADGROUP_normalized",
                "BARRIER_normalized",
                "TAIL_normalized",
            ]
        )

        for x in grid:
            writer.writerow(
                [
                    f"{x:.10f}",
                    f"{_linear_profile_value(profiles['HEADGROUP'], x):.10f}",
                    f"{_linear_profile_value(profiles['BARRIER'], x):.10f}",
                    f"{_linear_profile_value(profiles['TAIL'], x):.10f}",
                ]
            )


def make_dynamic_boundary_plot(
    profiles: dict,
    boundaries: dict,
    out_file: Path,
):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "[INFO] matplotlib is not installed; boundary PNG skipped."
        )
        return

    fig, ax = plt.subplots(
        figsize=(7.2, 5.0)
    )

    for region in (
        "HEADGROUP",
        "BARRIER",
        "TAIL",
    ):
        x = [
            item[0]
            for item in profiles[region]
        ]
        y = [
            item[1]
            for item in profiles[region]
        ]
        ax.plot(
            x,
            y,
            linewidth=1.2,
            label=region,
        )

    for key in (
        "tail_barrier_nm",
        "barrier_head_nm",
        "head_water_nm",
    ):
        ax.axvline(
            boundaries[key],
            linestyle="--",
            linewidth=1.0,
        )

    ax.set_xlabel(
        "Distance from membrane center (nm)"
    )
    ax.set_ylabel(
        "Normalized regional density"
    )
    ax.set_title(
        "Automatic membrane-region boundaries"
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(
        out_file,
        dpi=300,
    )
    plt.close(fig)

    print(
        f"[OK] Boundary plot generated: {out_file}"
    )


def prompt_manual_dynamic_boundaries(
    defaults: dict | None = None,
) -> dict:
    """
    Manual boundary input with strict increasing-order validation.
    """
    print("\n---------------- MANUAL REGION BOUNDARIES ----------------")
    print(
        "Enter positive distances from the membrane center:"
    )
    print(
        "  TAIL/BARRIER < BARRIER/HEADGROUP < HEADGROUP/WATER"
    )

    if defaults is None:
        d1 = 1.5
        d2 = 2.3
        d3 = 3.2
    else:
        d1 = defaults["tail_barrier_nm"]
        d2 = defaults["barrier_head_nm"]
        d3 = defaults["head_water_nm"]

    while True:
        b1 = ask_float(
            "TAIL / BARRIER boundary, nm",
            d1,
        )
        b2 = ask_float(
            "BARRIER / HEADGROUP boundary, nm",
            d2,
        )
        b3 = ask_float(
            "HEADGROUP / WATER boundary, nm",
            d3,
        )

        if (
            0 < b1 < b2 < b3
        ):
            return {
                "tail_barrier_nm": b1,
                "barrier_head_nm": b2,
                "head_water_nm": b3,
                "method": "manual",
                "tb_method": "manual",
                "bh_method": "manual",
                "hw_method": "manual",
            }

        print(
            "[WARNING] Boundaries must satisfy 0 < TB < BH < HW."
        )


def determine_dynamic_boundaries(
    density_paths: dict,
    settings: dict,
    subdir: Path,
) -> tuple[dict, dict]:
    """
    Recommended automatic method:
      - normalize each regional number-density profile independently
      - TAIL/BARRIER and BARRIER/HEADGROUP from normalized intersections
      - HEADGROUP/WATER from outer head-density decay threshold
      - fallback option: midpoint between regional density peaks
      - always allow manual override
    """
    profiles = {
        region: _normalized_positive_density_profile(
            density_paths[region]
        )
        for region in (
            "HEADGROUP",
            "BARRIER",
            "TAIL",
        )
    }

    peaks = {
        region: _profile_peak_position(
            profiles[region]
        )
        for region in (
            "HEADGROUP",
            "BARRIER",
            "TAIL",
        )
    }

    print("\n" + "=" * 78)
    print(" AUTOMATIC MEMBRANE-REGION BOUNDARY DETECTION")
    print("=" * 78)
    print(
        f" TAIL peak      : {peaks['TAIL']:.4f} nm from membrane center"
    )
    print(
        f" BARRIER peak   : {peaks['BARRIER']:.4f} nm from membrane center"
    )
    print(
        f" HEADGROUP peak : {peaks['HEADGROUP']:.4f} nm from membrane center"
    )
    print()

    print("Boundary method:")
    print(" [1] Normalized density intersections — recommended")
    print(" [2] Midpoints between regional density peaks")
    print(" [3] Manual distances from membrane center")

    while True:
        choice = ask(
            "Choose boundary method",
            "1",
        )

        if choice in {
            "1",
            "2",
            "3",
        }:
            break

        print(
            "[WARNING] Please choose 1, 2, or 3."
        )

    if choice == "3":
        boundaries = prompt_manual_dynamic_boundaries()
    else:
        peak_order_ok = (
            peaks["TAIL"]
            < peaks["BARRIER"]
            < peaks["HEADGROUP"]
        )

        if not peak_order_ok:
            print()
            print(
                "[WARNING] Regional density peaks are not ordered as "
                "TAIL < BARRIER < HEADGROUP."
            )
            print(
                "This can happen when region definitions overlap strongly "
                "or when the membrane is highly distorted."
            )
            print(
                "[ACTION] Automatic names are retained, but manual override "
                "is strongly recommended."
            )

        if choice == "1":
            tb, tb_method = _find_profile_intersection(
                profiles["TAIL"],
                profiles["BARRIER"],
                peaks["TAIL"],
                peaks["BARRIER"],
            )
            bh, bh_method = _find_profile_intersection(
                profiles["BARRIER"],
                profiles["HEADGROUP"],
                peaks["BARRIER"],
                peaks["HEADGROUP"],
            )
        else:
            tb = (
                peaks["TAIL"]
                + peaks["BARRIER"]
            ) / 2.0
            bh = (
                peaks["BARRIER"]
                + peaks["HEADGROUP"]
            ) / 2.0
            tb_method = "peak_midpoint"
            bh_method = "peak_midpoint"

        hw, hw_method = _find_head_outer_boundary(
            profiles["HEADGROUP"],
            peaks["HEADGROUP"],
            settings["head_decay"],
        )

        auto = {
            "tail_barrier_nm": tb,
            "barrier_head_nm": bh,
            "head_water_nm": hw,
            "method": (
                "normalized_density_intersection"
                if choice == "1"
                else "peak_midpoint"
            ),
            "tb_method": tb_method,
            "bh_method": bh_method,
            "hw_method": hw_method,
        }

        if not (
            0
            < auto["tail_barrier_nm"]
            < auto["barrier_head_nm"]
            < auto["head_water_nm"]
        ):
            print()
            print(
                "[WARNING] Automatic boundaries are not strictly ordered."
            )
            boundaries = prompt_manual_dynamic_boundaries(
                defaults=None,
            )
        else:
            print("\n[PROPOSED MEMBRANE REGIONS]")
            print(
                f" TAIL       : 0.0000 <= |DeltaZ| < "
                f"{auto['tail_barrier_nm']:.4f} nm"
            )
            print(
                f" BARRIER    : {auto['tail_barrier_nm']:.4f} <= |DeltaZ| < "
                f"{auto['barrier_head_nm']:.4f} nm"
            )
            print(
                f" HEADGROUP  : {auto['barrier_head_nm']:.4f} <= |DeltaZ| < "
                f"{auto['head_water_nm']:.4f} nm"
            )
            print(
                f" WATER      : |DeltaZ| >= "
                f"{auto['head_water_nm']:.4f} nm"
            )

            if yes_no(
                "Accept these automatically determined boundaries?",
                True,
            ):
                boundaries = auto
            else:
                boundaries = prompt_manual_dynamic_boundaries(
                    defaults=auto,
                )

    profile_csv = (
        subdir
        / "Membrane_region_boundary_profiles.csv"
    )

    write_dynamic_boundary_profile_csv(
        profile_csv,
        profiles,
    )

    plot_file = (
        subdir
        / "Membrane_region_boundaries.png"
    )

    make_dynamic_boundary_plot(
        profiles,
        boundaries,
        plot_file,
    )

    boundary_xlsx = (
        subdir
        / "Membrane_region_boundaries.xlsx"
    )

    headers = [
        "item",
        "value_nm",
        "method",
    ]

    rows = [
        [
            "TAIL density peak",
            peaks["TAIL"],
            "normalized density maximum",
        ],
        [
            "BARRIER density peak",
            peaks["BARRIER"],
            "normalized density maximum",
        ],
        [
            "HEADGROUP density peak",
            peaks["HEADGROUP"],
            "normalized density maximum",
        ],
        [
            "TAIL/BARRIER boundary",
            boundaries["tail_barrier_nm"],
            boundaries["tb_method"],
        ],
        [
            "BARRIER/HEADGROUP boundary",
            boundaries["barrier_head_nm"],
            boundaries["bh_method"],
        ],
        [
            "HEADGROUP/WATER boundary",
            boundaries["head_water_nm"],
            boundaries["hw_method"],
        ],
        [
            "HEADGROUP outer decay threshold",
            settings["head_decay"],
            "normalized fraction",
        ],
    ]

    write_excel_table(
        boundary_xlsx,
        "Region Boundaries",
        headers,
        rows,
    )

    print(
        f"[OK] Region-boundary workbook: {boundary_xlsx}"
    )

    return boundaries, peaks


def run_gmx_traj_com_xy(
    xtc: Path,
    tpr: Path,
    ndx_file: Path,
    group_name: str,
    out_file: Path,
    settings: dict,
):
    """
    Extract unwrapped COM X/Y coordinates using legacy gmx traj:
      -com
      -nojump
      -noz

    This is used for lateral displacement after subtracting membrane COM.
    """
    cmd = [
        GMX,
        "traj",
        "-f",
        str(xtc),
        "-s",
        str(tpr),
        "-n",
        str(ndx_file),
        "-ox",
        str(out_file),
        "-com",
        "-nojump",
        "-x",
        "-y",
        "-noz",
        "-ng",
        "1",
        *gmx_time_args(
            settings["begin_ps"],
            settings["end_ps"],
            settings["dt_ps"],
        ),
    ]

    result = _animated_process(
        cmd,
        input_text=group_name + "\n",
        label=f"COM XY: {group_name}",
    )

    log_file = timestamped_log_path(
        out_file.with_suffix(
            ".log"
        )
    )

    log_file.write_text(
        log_timestamp_header() + "COMMAND\n"
        + " ".join(map(str, cmd))
        + "\n\nGROUP\n"
        + group_name
        + "\n\nSTDOUT\n"
        + (result.stdout or "")
        + "\nSTDERR\n"
        + (result.stderr or ""),
        encoding="utf-8",
    )

    if result.returncode != 0:
        raise RecoverableAnalysisError(
            f"gmx traj failed for '{group_name}'. "
            f"See {log_file}"
        )

    rows = parse_xvg(
        out_file,
        min_cols=3,
    )

    if not rows:
        raise RecoverableAnalysisError(
            f"No COM XY coordinates were generated for '{group_name}'."
        )

    return rows


def run_gmx_distance_delta_xyz(
    xtc: Path,
    tpr: Path,
    ndx_file: Path,
    membrane_group: str,
    drug_group: str,
    out_file: Path,
    settings: dict,
):
    """
    Use gmx distance -oxyz with PBC enabled to obtain minimum-image
    membrane-COM <-> drug-COM displacement components. The signed Z component
    is used for region assignment; |DeltaZ| is the membrane depth.
    """
    selection = (
        f'com of group "{membrane_group}" '
        f'plus com of group "{drug_group}"'
    )

    cmd = [
        GMX,
        "distance",
        "-f",
        str(xtc),
        "-s",
        str(tpr),
        "-n",
        str(ndx_file),
        "-oxyz",
        str(out_file),
        "-pbc",
        "-select",
        selection,
        *gmx_time_args(
            settings["begin_ps"],
            settings["end_ps"],
            settings["dt_ps"],
        ),
    ]

    result = _animated_process(
        cmd,
        input_text=None,
        label=f"Relative XYZ: {drug_group}",
    )

    log_file = timestamped_log_path(
        out_file.with_suffix(
            ".log"
        )
    )

    log_file.write_text(
        log_timestamp_header() + "COMMAND\n"
        + " ".join(map(str, cmd))
        + "\n\nSELECTION\n"
        + selection
        + "\n\nSTDOUT\n"
        + (result.stdout or "")
        + "\nSTDERR\n"
        + (result.stderr or ""),
        encoding="utf-8",
    )

    if result.returncode != 0:
        raise RecoverableAnalysisError(
            f"gmx distance failed for '{drug_group}'. "
            f"See {log_file}"
        )

    rows = parse_xvg(
        out_file,
        min_cols=4,
    )

    if not rows:
        raise RecoverableAnalysisError(
            f"No relative XYZ data were generated for '{drug_group}'."
        )

    return rows


def _series_by_rounded_time(
    rows: list[list[float]],
    needed_cols: int,
) -> dict[float, list[float]]:
    """
    GROMACS tools run with identical -b/-e/-dt should have matching times.
    Rounding protects against harmless printed floating-point differences.
    """
    result = {}

    for row in rows:
        if len(row) < needed_cols:
            continue

        key = round(
            row[0],
            6,
        )
        result[key] = row

    return result


def align_dynamic_drug_series(
    membrane_xy: list[list[float]],
    drug_xy: list[list[float]],
    delta_xyz: list[list[float]],
) -> list[dict]:
    """
    Align:
      membrane unwrapped COM XY
      drug unwrapped COM XY
      minimum-image membrane->drug XYZ

    X/Y are expressed relative to membrane COM to remove whole-membrane drift.
    Z is taken from gmx distance minimum-image displacement.
    """
    mem = _series_by_rounded_time(
        membrane_xy,
        3,
    )
    drug = _series_by_rounded_time(
        drug_xy,
        3,
    )
    delta = _series_by_rounded_time(
        delta_xyz,
        4,
    )

    common = sorted(
        set(mem)
        & set(drug)
        & set(delta)
    )

    if len(common) < 3:
        raise RecoverableAnalysisError(
            "Too few common frames among membrane COM, drug COM, "
            "and membrane-drug distance trajectories."
        )

    series = []

    for key in common:
        m = mem[key]
        d = drug[key]
        xyz = delta[key]

        series.append(
            {
                "time_ps": key,
                "x_nm": d[1] - m[1],
                "y_nm": d[2] - m[2],
                "z_nm": xyz[3],
                "abs_z_nm": abs(xyz[3]),
            }
        )

    return series


def classify_dynamic_region(
    distance_nm: float,
    boundaries: dict,
) -> str:
    if distance_nm < boundaries["tail_barrier_nm"]:
        return "TAIL"
    if distance_nm < boundaries["barrier_head_nm"]:
        return "BARRIER"
    if distance_nm < boundaries["head_water_nm"]:
        return "HEADGROUP"
    return "WATER"


def apply_dynamic_hysteresis(
    distances: list[float],
    boundaries: dict,
    hysteresis_nm: float,
) -> list[str]:
    """
    State-dependent boundary hysteresis prevents rapid HEAD/BARRIER/TAIL
    switching when a drug fluctuates around a regional boundary.
    """
    if not distances:
        return []

    states = [
        classify_dynamic_region(
            distances[0],
            boundaries,
        )
    ]

    b1 = boundaries["tail_barrier_nm"]
    b2 = boundaries["barrier_head_nm"]
    b3 = boundaries["head_water_nm"]
    h = hysteresis_nm

    for d in distances[1:]:
        previous = states[-1]
        base = classify_dynamic_region(
            d,
            boundaries,
        )

        if previous == "TAIL":
            if d >= b1 + h:
                state = base
            else:
                state = "TAIL"

        elif previous == "BARRIER":
            if d < b1 - h:
                state = base
            elif d >= b2 + h:
                state = base
            else:
                state = "BARRIER"

        elif previous == "HEADGROUP":
            if d < b2 - h:
                state = base
            elif d >= b3 + h:
                state = base
            else:
                state = "HEADGROUP"

        else:  # WATER
            if d < b3 - h:
                state = base
            else:
                state = "WATER"

        states.append(state)

    return states


def dynamic_segments(
    times: list[float],
    states: list[str],
) -> list[dict]:
    if not times or not states:
        return []

    segments = []
    start = 0

    for i in range(
        1,
        len(states) + 1,
    ):
        if (
            i == len(states)
            or states[i] != states[start]
        ):
            end = i - 1
            duration = (
                times[end]
                - times[start]
            )
            segments.append(
                {
                    "region": states[start],
                    "start_idx": start,
                    "end_idx": end,
                    "start_ps": times[start],
                    "end_ps": times[end],
                    "duration_ps": max(
                        0.0,
                        duration,
                    ),
                    "frames": end - start + 1,
                }
            )
            start = i

    return segments


def merge_short_dynamic_excursions(
    times: list[float],
    states: list[str],
    min_residence_ps: float,
) -> list[str]:
    """
    If a short intermediate segment lies between two segments of the SAME
    region, treat it as a brief boundary excursion and merge it back.
    Repeat until no such short excursion remains.
    """
    if (
        min_residence_ps <= 0
        or len(states) < 3
    ):
        return list(states)

    states = list(states)

    changed = True

    while changed:
        changed = False
        segments = dynamic_segments(
            times,
            states,
        )

        if len(segments) < 3:
            break

        for i in range(
            1,
            len(segments) - 1,
        ):
            current = segments[i]
            previous = segments[i - 1]
            following = segments[i + 1]

            if (
                current["duration_ps"]
                < min_residence_ps
                and previous["region"]
                == following["region"]
            ):
                for idx in range(
                    current["start_idx"],
                    current["end_idx"] + 1,
                ):
                    states[idx] = previous["region"]

                changed = True
                break

    return states


def write_dynamic_assignment_csv(
    path: Path,
    drug_series: list[dict],
):
    """
    Per-frame audit trail for region classification.
    """
    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "time_ps",
                "time_ns",
                "drug_label",
                "drug_resname",
                "x_relative_membrane_nm",
                "y_relative_membrane_nm",
                "z_relative_membrane_nm",
                "abs_z_nm",
                "raw_region",
                "final_region",
            ]
        )

        for entry in drug_series:
            label = entry["target"]["label"]
            resname = entry["target"]["resname"]

            for row in entry["series"]:
                writer.writerow(
                    [
                        f"{row['time_ps']:.6f}",
                        f"{row['time_ps'] / 1000.0:.9f}",
                        label,
                        resname,
                        f"{row['x_nm']:.9f}",
                        f"{row['y_nm']:.9f}",
                        f"{row['z_nm']:.9f}",
                        f"{row['abs_z_nm']:.9f}",
                        row["raw_region"],
                        row["final_region"],
                    ]
                )


def collect_dynamic_residence_rows(
    drug_series: list[dict],
    min_residence_ps: float,
) -> tuple[list[list], list[dict]]:
    """
    Return Excel rows plus segment objects used for MSD aggregation.
    """
    excel_rows = []
    usable = []

    for entry in drug_series:
        times = [
            row["time_ps"]
            for row in entry["series"]
        ]
        states = [
            row["final_region"]
            for row in entry["series"]
        ]

        segments = dynamic_segments(
            times,
            states,
        )

        for number, seg in enumerate(
            segments,
            start=1,
        ):
            use_for_msd = (
                seg["region"]
                in {
                    "HEADGROUP",
                    "BARRIER",
                    "TAIL",
                }
                and seg["duration_ps"]
                >= min_residence_ps
                and seg["frames"] >= 3
            )

            excel_rows.append(
                [
                    entry["target"]["label"],
                    entry["target"]["resname"],
                    number,
                    seg["region"],
                    seg["start_ps"],
                    seg["end_ps"],
                    seg["duration_ps"],
                    seg["duration_ps"] / 1000.0,
                    seg["frames"],
                    str(use_for_msd),
                ]
            )

            if use_for_msd:
                usable.append(
                    {
                        "target": entry["target"],
                        "series": entry["series"],
                        **seg,
                    }
                )

    return excel_rows, usable


def _median_frame_dt(
    series_collection: list[dict],
) -> float:
    intervals = []

    for entry in series_collection:
        times = [
            row["time_ps"]
            for row in entry["series"]
        ]

        for i in range(
            1,
            len(times),
        ):
            dt = (
                times[i]
                - times[i - 1]
            )
            if dt > 0:
                intervals.append(dt)

    if not intervals:
        raise RecoverableAnalysisError(
            "Could not determine trajectory frame interval."
        )

    return median(intervals)


def _dynamic_sq_displacement(
    a: dict,
    b: dict,
    dimension: str,
) -> float:
    dx = b["x_nm"] - a["x_nm"]
    dy = b["y_nm"] - a["y_nm"]
    dz = b["z_nm"] - a["z_nm"]

    if dimension == "XY":
        return dx * dx + dy * dy

    if dimension == "Z":
        return dz * dz

    return (
        dx * dx
        + dy * dy
        + dz * dz
    )


def aggregate_dynamic_region_msd(
    usable_segments: list[dict],
    region: str,
    frame_dt_ps: float,
    settings: dict,
) -> tuple[list[dict], dict]:
    """
    Multiple-time-origin MSD across all continuous residence segments belonging
    to one membrane region.

    No discontinuous segments are concatenated.
    """
    sums = {}
    counts = {}
    contributing_segments = 0
    total_residence_ps = 0.0

    for seg in usable_segments:
        if seg["region"] != region:
            continue

        start = seg["start_idx"]
        end = seg["end_idx"]
        series = seg["series"]
        n = end - start + 1

        if n < 3:
            continue

        contributing_segments += 1
        total_residence_ps += seg["duration_ps"]

        auto_max_frames = max(
            1,
            (n - 1) // 2,
        )

        if settings["max_lag_ps"] is None:
            max_frames = auto_max_frames
        else:
            requested = max(
                1,
                int(
                    settings["max_lag_ps"]
                    / frame_dt_ps
                ),
            )
            max_frames = min(
                n - 1,
                requested,
            )

        for lag in range(
            1,
            max_frames + 1,
        ):
            lag_ps = round(
                lag * frame_dt_ps,
                6,
            )

            for i in range(
                start,
                end - lag + 1,
            ):
                j = i + lag

                value = _dynamic_sq_displacement(
                    series[i],
                    series[j],
                    settings["dimension"],
                )

                sums[lag_ps] = (
                    sums.get(
                        lag_ps,
                        0.0,
                    )
                    + value
                )
                counts[lag_ps] = (
                    counts.get(
                        lag_ps,
                        0,
                    )
                    + 1
                )

    curve = [
        {
            "lag_ps": lag,
            "lag_ns": lag / 1000.0,
            "msd_nm2": sums[lag] / counts[lag],
            "origins": counts[lag],
        }
        for lag in sorted(sums)
        if counts[lag] > 0
    ]

    meta = {
        "segments": contributing_segments,
        "total_residence_ps": total_residence_ps,
        "total_residence_ns": total_residence_ps / 1000.0,
    }

    return curve, meta


def fit_dynamic_diffusion(
    curve: list[dict],
    settings: dict,
) -> dict:
    """
    Ordinary least-squares fit:
      MSD = slope * lag + intercept

    Diffusion conversion:
      XY: D = slope / 4
      3D: D = slope / 6
      Z : D = slope / 2

    slope is nm^2/ps. Multiplying D[nm^2/ps] by 1000 gives nm^2/ns,
    numerically equal to units of 1e-5 cm^2/s.

    Slope standard error is descriptive only; independent-trajectory or
    block-bootstrap uncertainty is still preferable for publication.
    """
    if len(curve) < 3:
        return {
            "success": False,
            "reason": "fewer than 3 MSD lag points",
        }

    max_lag = curve[-1]["lag_ps"]

    if settings["auto_fit"]:
        fit_start = max_lag * 0.10
        fit_end = max_lag * 0.50
    else:
        fit_start = settings["fit_start_ps"]
        fit_end = settings["fit_end_ps"]

    points = [
        row
        for row in curve
        if (
            fit_start
            <= row["lag_ps"]
            <= fit_end
        )
    ]

    if len(points) < 3:
        points = curve[
            0:max(
                3,
                min(
                    len(curve),
                    len(curve) // 2,
                ),
            )
        ]

    if len(points) < 3:
        return {
            "success": False,
            "reason": "insufficient fit points",
        }

    x = [
        row["lag_ps"]
        for row in points
    ]
    y = [
        row["msd_nm2"]
        for row in points
    ]

    xbar = mean(x)
    ybar = mean(y)

    sxx = sum(
        (xi - xbar) ** 2
        for xi in x
    )

    if sxx <= 0:
        return {
            "success": False,
            "reason": "zero variance in lag times",
        }

    slope = sum(
        (
            xi - xbar
        )
        * (
            yi - ybar
        )
        for xi, yi in zip(x, y)
    ) / sxx

    intercept = (
        ybar
        - slope * xbar
    )

    residuals = [
        yi
        - (
            slope * xi
            + intercept
        )
        for xi, yi in zip(x, y)
    ]

    if len(points) > 2:
        residual_var = sum(
            r * r
            for r in residuals
        ) / (
            len(points) - 2
        )
        slope_se = math.sqrt(
            residual_var / sxx
        )
    else:
        slope_se = float("nan")

    divisor = settings["diffusion_divisor"]

    d_nm2_ps = slope / divisor
    d_nm2_ns = d_nm2_ps * 1000.0
    d_err_nm2_ns = (
        slope_se / divisor * 1000.0
        if math.isfinite(slope_se)
        else float("nan")
    )

    y_pred = [
        slope * xi + intercept
        for xi in x
    ]

    ss_res = sum(
        (
            yi - pi
        ) ** 2
        for yi, pi in zip(y, y_pred)
    )
    ss_tot = sum(
        (
            yi - ybar
        ) ** 2
        for yi in y
    )

    if ss_tot > 0:
        r2 = 1.0 - ss_res / ss_tot
    else:
        r2 = float("nan")

    return {
        "success": True,
        "slope_nm2_ps": slope,
        "slope_se_nm2_ps": slope_se,
        "intercept_nm2": intercept,
        "D_nm2_ps": d_nm2_ps,
        "D_nm2_ns": d_nm2_ns,
        "D_1e-5_cm2_s": d_nm2_ns,
        "D_error_1e-5_cm2_s": d_err_nm2_ns,
        "fit_start_ps": min(x),
        "fit_end_ps": max(x),
        "fit_points": len(points),
        "r2": r2,
    }


def write_dynamic_msd_xvg(
    path: Path,
    region: str,
    curve: list[dict],
    fit: dict,
    settings: dict,
):
    with path.open(
        "w",
        encoding="utf-8",
    ) as fh:
        fh.write(
            "# Region-specific MSD from continuous residence segments\n"
        )
        fh.write(
            "# Discontinuous residence periods are NOT concatenated\n"
        )
        fh.write(
            f"# Diffusion dimension: {settings['dimension']}\n"
        )
        if fit.get("success"):
            fh.write(
                f"# D = {fit['D_1e-5_cm2_s']:.10g} "
                "(1e-5 cm^2/s)\n"
            )
            fh.write(
                f"# Approximate fit SE = "
                f"{fit['D_error_1e-5_cm2_s']:.10g} "
                "(1e-5 cm^2/s)\n"
            )

        fh.write(
            f'@ title "Drug MSD in {region}"\n'
        )
        fh.write(
            '@ xaxis label "Lag time (ns)"\n'
        )
        fh.write(
            '@ yaxis label "MSD (nm\\S2\\N)"\n'
        )
        fh.write(
            '@ s0 legend "MSD"\n'
        )

        for row in curve:
            fh.write(
                f"{row['lag_ns']:.10f} "
                f"{row['msd_nm2']:.10f}\n"
            )


def make_dynamic_msd_plot(
    path: Path,
    region: str,
    curve: list[dict],
    fit: dict,
):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return

    if not curve:
        return

    x = [
        row["lag_ns"]
        for row in curve
    ]
    y = [
        row["msd_nm2"]
        for row in curve
    ]

    fig, ax = plt.subplots(
        figsize=(7.2, 5.0)
    )

    ax.plot(
        x,
        y,
        linewidth=1.2,
    )

    if fit.get("success"):
        x_fit_ps = [
            fit["fit_start_ps"],
            fit["fit_end_ps"],
        ]
        y_fit = [
            (
                fit["slope_nm2_ps"] * xp
                + fit["intercept_nm2"]
            )
            for xp in x_fit_ps
        ]
        ax.plot(
            [
                xp / 1000.0
                for xp in x_fit_ps
            ],
            y_fit,
            linestyle="--",
            linewidth=1.0,
        )

    ax.set_title(
        f"Drug MSD — {region}"
    )
    ax.set_xlabel(
        "Lag time (ns)"
    )
    ax.set_ylabel(
        "MSD (nm²)"
    )
    fig.tight_layout()
    fig.savefig(
        path,
        dpi=300,
    )
    plt.close(fig)


def write_dynamic_diffusion_summary(
    subdir: Path,
    results: list[dict],
    settings: dict,
):
    """
    Final region-wise Excel summary.
    """
    path = (
        subdir
        / "Drug_region_diffusion_summary.xlsx"
    )

    headers = [
        "region",
        "dimension",
        "D_1e-5_cm2_s",
        "approx_fit_SE_1e-5_cm2_s",
        "slope_nm2_ps",
        "fit_start_ps",
        "fit_end_ps",
        "fit_points",
        "R2",
        "continuous_segments",
        "total_residence_ns",
        "MSD_curve",
    ]

    rows = []

    for result in results:
        fit = result["fit"]
        meta = result["meta"]

        rows.append(
            [
                result["region"],
                settings["dimension"],
                (
                    fit.get(
                        "D_1e-5_cm2_s",
                        "NA",
                    )
                    if fit.get("success")
                    else "NA"
                ),
                (
                    fit.get(
                        "D_error_1e-5_cm2_s",
                        "NA",
                    )
                    if fit.get("success")
                    else "NA"
                ),
                (
                    fit.get(
                        "slope_nm2_ps",
                        "NA",
                    )
                    if fit.get("success")
                    else "NA"
                ),
                (
                    fit.get(
                        "fit_start_ps",
                        "NA",
                    )
                    if fit.get("success")
                    else "NA"
                ),
                (
                    fit.get(
                        "fit_end_ps",
                        "NA",
                    )
                    if fit.get("success")
                    else "NA"
                ),
                (
                    fit.get(
                        "fit_points",
                        "NA",
                    )
                    if fit.get("success")
                    else "NA"
                ),
                (
                    fit.get(
                        "r2",
                        "NA",
                    )
                    if fit.get("success")
                    else "NA"
                ),
                meta["segments"],
                meta["total_residence_ns"],
                str(result["xvg"]),
            ]
        )

    write_excel_table(
        path,
        "Region Diffusion",
        headers,
        rows,
    )

    return path


def analyze_dynamic_drug_region_diffusion(
    gro: Path,
    tpr: Path,
    xtc: Path,
    settings: dict,
):
    """
    Full dynamic-region workflow:
      membrane density boundaries
      -> individual drug tracking
      -> dynamic region assignment
      -> hysteresis
      -> continuous residence segments
      -> multiple-time-origin MSD
      -> region-specific diffusion coefficients
    """
    subdir = DRUG_REGION_DIR
    subdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # 1. Drug/guest molecules.
    drug_resname, drug_targets = prompt_dynamic_drug_type(
        gro,
        tpr,
    )

    # 2. Membrane lipids and region definitions.
    print("\n" + "=" * 78)
    print(" MEMBRANE DEFINITION FOR DYNAMIC DRUG-REGION ANALYSIS")
    print("=" * 78)

    lipid_names = prompt_membrane_lipid_names(
        gro,
        tpr,
    )

    definitions = prompt_all_membrane_regions_batch(
        gro,
        lipid_names,
        tpr,
    )

    membrane_indices, region_indices = dynamic_region_indices(
        gro,
        lipid_names,
        definitions,
    )

    ndx_file = (
        subdir
        / "Dynamic_drug_region_analysis.ndx"
    )

    write_dynamic_region_index(
        ndx_file,
        membrane_indices,
        region_indices,
        drug_targets,
    )

    # 3. Automatic membrane-region boundaries.
    print("\n" + "=" * 78)
    print(" STEP 1 — AUTOMATIC MEMBRANE-REGION BOUNDARIES")
    print("=" * 78)

    density_paths = run_dynamic_boundary_density_profiles(
        xtc,
        tpr,
        ndx_file,
        subdir,
        settings,
        region_indices,
    )

    boundaries, peaks = determine_dynamic_boundaries(
        density_paths,
        settings,
        subdir,
    )

    # 4. Extract membrane COM XY only once.
    print("\n" + "=" * 78)
    print(" STEP 2 — INDIVIDUAL DRUG TRAJECTORY TRACKING")
    print("=" * 78)

    membrane_xy_file = (
        subdir
        / "Dynamic_membrane_center_COM_XY.xvg"
    )

    membrane_xy = run_gmx_traj_com_xy(
        xtc,
        tpr,
        ndx_file,
        "DYN_MEMBRANE_CENTER",
        membrane_xy_file,
        settings,
    )

    drug_series = []

    for number, target in enumerate(
        drug_targets,
        start=1,
    ):
        print(
            f"\n---------------- TRACKING DRUG "
            f"{number}/{len(drug_targets)}: "
            f"{target['label']} ----------------"
        )

        xy_file = (
            subdir
            / f"{target['group']}_COM_XY.xvg"
        )

        xyz_file = (
            subdir
            / f"{target['group']}_relative_XYZ.xvg"
        )

        drug_xy = run_gmx_traj_com_xy(
            xtc,
            tpr,
            ndx_file,
            target["group"],
            xy_file,
            settings,
        )

        delta_xyz = run_gmx_distance_delta_xyz(
            xtc,
            tpr,
            ndx_file,
            "DYN_MEMBRANE_CENTER",
            target["group"],
            xyz_file,
            settings,
        )

        series = align_dynamic_drug_series(
            membrane_xy,
            drug_xy,
            delta_xyz,
        )

        distances = [
            row["abs_z_nm"]
            for row in series
        ]

        raw_states = [
            classify_dynamic_region(
                d,
                boundaries,
            )
            for d in distances
        ]

        final_states = apply_dynamic_hysteresis(
            distances,
            boundaries,
            settings["hysteresis_nm"],
        )

        times = [
            row["time_ps"]
            for row in series
        ]

        final_states = merge_short_dynamic_excursions(
            times,
            final_states,
            settings["min_residence_ps"],
        )

        for row, raw_region, final_region in zip(
            series,
            raw_states,
            final_states,
        ):
            row["raw_region"] = raw_region
            row["final_region"] = final_region

        drug_series.append(
            {
                "target": target,
                "series": series,
            }
        )

    # 5. Assignment audit trail.
    assignment_csv = (
        subdir
        / "Drug_region_assignment.csv"
    )

    write_dynamic_assignment_csv(
        assignment_csv,
        drug_series,
    )

    # 6. Continuous residence segments.
    residence_rows, usable_segments = collect_dynamic_residence_rows(
        drug_series,
        settings["min_residence_ps"],
    )

    residence_xlsx = (
        subdir
        / "Drug_region_residence_segments.xlsx"
    )

    residence_headers = [
        "drug_label",
        "drug_resname",
        "segment_number",
        "region",
        "start_ps",
        "end_ps",
        "duration_ps",
        "duration_ns",
        "frames",
        "used_for_MSD",
    ]

    write_excel_table(
        residence_xlsx,
        "Residence Segments",
        residence_headers,
        residence_rows,
    )

    frame_dt_ps = _median_frame_dt(
        drug_series,
    )

    print("\n" + "=" * 78)
    print(" STEP 3 — CONTINUOUS-RESIDENCE MSD / DIFFUSION")
    print("=" * 78)
    print(
        f"Median analyzed frame interval : {frame_dt_ps:.6g} ps"
    )
    print(
        f"Minimum residence for MSD      : "
        f"{settings['min_residence_ps']:.6g} ps"
    )
    print(
        f"Diffusion dimension            : {settings['dimension']}"
    )
    print(
        f"Valid continuous segments      : {len(usable_segments)}"
    )

    results = []

    for region in (
        "HEADGROUP",
        "BARRIER",
        "TAIL",
    ):
        curve, meta = aggregate_dynamic_region_msd(
            usable_segments,
            region,
            frame_dt_ps,
            settings,
        )

        fit = fit_dynamic_diffusion(
            curve,
            settings,
        )

        xvg = (
            subdir
            / f"MSD_{drug_resname}_{region}.xvg"
        )
        png = (
            subdir
            / f"MSD_{drug_resname}_{region}.png"
        )

        write_dynamic_msd_xvg(
            xvg,
            region,
            curve,
            fit,
            settings,
        )

        make_dynamic_msd_plot(
            png,
            region,
            curve,
            fit,
        )

        result = {
            "region": region,
            "curve": curve,
            "fit": fit,
            "meta": meta,
            "xvg": xvg,
            "png": png,
        }

        results.append(result)

        print()
        print(f"[{region}]")
        print(
            f"  continuous segments : {meta['segments']}"
        )
        print(
            f"  total residence     : "
            f"{meta['total_residence_ns']:.6g} ns"
        )

        if fit.get("success"):
            print(
                f"  D[{settings['dimension']}] "
                f"= {fit['D_1e-5_cm2_s']:.6g} "
                f"(+/- {fit['D_error_1e-5_cm2_s']:.6g}) "
                "(1e-5 cm^2/s)"
            )
            print(
                f"  linear-fit R^2      : "
                f"{fit['r2']:.6g}"
            )
        else:
            print(
                f"  D                   : NOT FITTED "
                f"({fit.get('reason', 'unknown reason')})"
            )

    summary_xlsx = write_dynamic_diffusion_summary(
        subdir,
        results,
        settings,
    )

    # 7. Overall settings / method workbook.
    settings_xlsx = (
        subdir
        / "Dynamic_region_analysis_settings.xlsx"
    )

    settings_rows = [
        [
            "drug_resname",
            drug_resname,
        ],
        [
            "drug_molecules_analyzed",
            len(drug_targets),
        ],
        [
            "membrane_lipids",
            ", ".join(lipid_names),
        ],
        [
            "TAIL/BARRIER_boundary_nm",
            boundaries["tail_barrier_nm"],
        ],
        [
            "BARRIER/HEADGROUP_boundary_nm",
            boundaries["barrier_head_nm"],
        ],
        [
            "HEADGROUP/WATER_boundary_nm",
            boundaries["head_water_nm"],
        ],
        [
            "boundary_method",
            boundaries["method"],
        ],
        [
            "boundary_hysteresis_nm",
            settings["hysteresis_nm"],
        ],
        [
            "minimum_residence_ps",
            settings["min_residence_ps"],
        ],
        [
            "diffusion_dimension",
            settings["dimension"],
        ],
        [
            "frame_interval_ps",
            frame_dt_ps,
        ],
        [
            "max_MSD_lag_ps",
            (
                settings["max_lag_ps"]
                if settings["max_lag_ps"] is not None
                else "automatic: 50% of each residence segment"
            ),
        ],
        [
            "fit_mode",
            (
                "automatic: 10%-50% of available lag"
                if settings["auto_fit"]
                else (
                    f"{settings['fit_start_ps']} to "
                    f"{settings['fit_end_ps']} ps"
                )
            ),
        ],
        [
            "important_method_note",
            (
                "MSD uses continuous residence segments only; "
                "discontinuous visits are never concatenated."
            ),
        ],
    ]

    write_excel_table(
        settings_xlsx,
        "Analysis Settings",
        [
            "parameter",
            "value",
        ],
        settings_rows,
    )

    print("\n" + "=" * 78)
    print(" DYNAMIC DRUG-REGION DIFFUSION COMPLETED")
    print("=" * 78)
    print(
        f" Per-frame region assignment : {assignment_csv}"
    )
    print(
        f" Residence segments          : {residence_xlsx}"
    )
    print(
        f" Membrane boundaries         : "
        f"{subdir / 'Membrane_region_boundaries.xlsx'}"
    )
    print(
        f" Region diffusion summary    : {summary_xlsx}"
    )
    print(
        f" Analysis settings           : {settings_xlsx}"
    )
    print("=" * 78)


def run_dynamic_drug_region_menu(
    state: dict,
) -> bool:
    """
    Menu wrapper for module 7.
    """
    print("\n" + "=" * 78)
    print(" DRUG MSD / DIFFUSION BY DYNAMIC MEMBRANE REGION")
    print("=" * 78)
    print(
        "This analysis tracks EACH drug molecule independently and assigns "
        "it frame-by-frame to HEADGROUP, BARRIER, TAIL, or WATER according "
        "to its minimum-image Z distance from the membrane center."
    )
    print()
    print(
        "Membrane-region boundaries are derived from normalized regional "
        "density profiles, then continuous residence segments are detected."
    )
    print()
    print(
        "MSD is calculated only WITHIN continuous residence segments. "
        "Separate visits to the same region are never concatenated."
    )
    print()
    print(
        "Recommended output: XY lateral diffusion, because the regions "
        "themselves are defined by Z depth."
    )

    files = resolve_analysis_files(
        [
            (
                "gro",
                "gro",
                "GRO file",
            ),
            (
                "tpr",
                "tpr",
                "TPR file",
            ),
            (
                "xtc",
                "xtc",
                "XTC trajectory",
            ),
        ],
        state,
    )

    gro = files["gro"]
    tpr = files["tpr"]
    xtc = files["xtc"]

    while True:
        settings = prompt_dynamic_region_settings()

        analyze_dynamic_drug_region_diffusion(
            gro,
            tpr,
            xtc,
            settings,
        )

        action = post_calculation_menu(
            "DRUG DYNAMIC-REGION MSD / DIFFUSION",
            "Calculate another dynamic drug-region diffusion analysis",
        )

        if action == "menu":
            return False

        if action == "exit":
            return True

        reuse = yes_no(
            "Reuse the current GRO/TPR/XTC files?",
            True,
        )

        if reuse:
            continue

        files = resolve_analysis_files(
            [
                (
                    "gro",
                    "gro",
                    "GRO file",
                ),
                (
                    "tpr",
                    "tpr",
                    "TPR file",
                ),
                (
                    "xtc",
                    "xtc",
                    "XTC trajectory",
                ),
            ],
            state,
        )

        gro = files["gro"]
        tpr = files["tpr"]
        xtc = files["xtc"]





# ======================================================================
# 8) RMSD — STRUCTURAL ROOT-MEAN-SQUARE DEVIATION
# ======================================================================
#
# This module is ADDITIVE. All existing APL / thickness / SASA / order /
# MSD / density / dynamic-region functions remain unchanged.
#
# Supported AA/CG target modes:
#   1) Existing/default GROMACS index group(s)
#      e.g. Protein, Backbone, C-alpha, System
#   2) Complete residue/molecule TYPE(S) selected by resname
#   3) Specific atom(s) (AA) / bead(s) (CG) inside selected resname(s)
#   4) Specific individual molecule/residue INSTANCE(S)
#   5) Custom merged group containing multiple molecular selections,
#      with a user-defined renamed NDX group
#
# For every RMSD target:
#   - fitting group can equal the RMSD group (simple default)
#   - or one COMMON independent fitting group can be selected
#
# GROMACS 2020.x method:
#   gmx rms -s REF -f TRAJ -n INDEX -o RMSD.xvg
#
# The fitting group and RMSD calculation group do not have to be identical.
# ======================================================================


def read_ndx_groups(
    ndx_file: Path,
) -> list[dict]:
    """
    Parse a GROMACS .ndx file while preserving group order.

    Returns:
        [
          {"name": "System", "indices": [...]},
          {"name": "Protein", "indices": [...]},
          ...
        ]
    """
    groups = []
    current = None

    with ndx_file.open(
        encoding="utf-8",
        errors="replace",
    ) as fh:
        for line in fh:
            s = line.strip()

            if not s:
                continue

            if s.startswith("[") and s.endswith("]"):
                name = s[1:-1].strip()
                current = {
                    "name": name,
                    "indices": [],
                }
                groups.append(current)
                continue

            if current is None:
                continue

            for token in s.split():
                try:
                    current["indices"].append(
                        int(token)
                    )
                except ValueError:
                    pass

    return groups


def build_rmsd_default_index(
    structure_file: Path,
) -> Path:
    """
    Generate a normal GROMACS index for convenient selection of standard
    groups such as Protein / Backbone / C-alpha / System.

    This index is RMSD-local and does not overwrite the script's global
    index.ndx.
    """
    RMSD_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    ndx_file = (
        RMSD_DIR
        / "RMSD_default_groups.ndx"
    )

    run(
        [
            GMX,
            "make_ndx",
            "-f",
            str(structure_file),
            "-o",
            str(ndx_file),
        ],
        input_text="q\n",
        quiet=True,
        progress_label="RMSD: generating default GROMACS groups",
    )

    if (
        not ndx_file.exists()
        or ndx_file.stat().st_size == 0
    ):
        raise RecoverableAnalysisError(
            "Could not generate the default RMSD index."
        )

    return ndx_file


def show_numbered_ndx_groups_for_rmsd(
    ndx_file: Path,
    title: str,
) -> list[dict]:
    """
    Display all GROMACS index groups as numbered options.
    """
    groups = read_ndx_groups(
        ndx_file,
    )

    if not groups:
        raise RecoverableAnalysisError(
            f"No index groups could be read from {ndx_file}."
        )

    print("\n" + "-" * 78)
    print(f" {title}")
    print("-" * 78)
    print(
        f"{'No.':>6s}  {'Group name':<32s} {'atoms/beads':>12s}"
    )
    print(
        "  " + "-" * 60
    )

    for number, group in enumerate(
        groups,
        start=1,
    ):
        print(
            f"  [{number:>3d}]  "
            f"{group['name']:<32.32s} "
            f"{len(group['indices']):>12d}"
        )

    return groups


def prompt_rmsd_default_index_targets(
    default_ndx: Path,
) -> list[dict]:
    """
    Select one or multiple existing/default GROMACS groups by number.
    Each selected group becomes one independent RMSD target.
    """
    while True:
        groups = show_numbered_ndx_groups_for_rmsd(
            default_ndx,
            "RMSD — GROMACS INDEX GROUP SELECTION",
        )

        print()
        print(
            "Select one or more groups. Each group is calculated separately."
        )
        print(
            "Examples: Protein, Backbone, C-alpha, System, or any custom "
            "group generated by GROMACS."
        )
        print()
        print("Input examples: 2 4 5   or   2,4,5")

        raw = ask(
            "Select RMSD group number(s)"
        )

        numbers, invalid = _parse_numbered_selection(
            raw,
            len(groups),
        )

        if invalid or not numbers:
            if invalid:
                print(
                    "[WARNING] Invalid group number(s): "
                    + ", ".join(invalid)
                )
            print(
                f"[ACTION] Enter one or more numbers from 1-{len(groups)}."
            )
            continue

        selected = [
            groups[number - 1]
            for number in numbers
        ]

        print("\n[SELECTED RMSD GROUPS]")
        for number in numbers:
            group = groups[number - 1]
            print(
                f"  [{number:>3d}] {group['name']} "
                f"({len(group['indices'])} atoms/beads)"
            )

        if not yes_no(
            "Confirm these RMSD groups?",
            True,
        ):
            continue

        targets = []

        for group in selected:
            targets.append(
                {
                    "label": group["name"],
                    "group": (
                        "RMSD_CALC_"
                        + sanitize_filename(
                            group["name"]
                        )
                    ),
                    "indices": list(
                        group["indices"]
                    ),
                    "mode": "gromacs_index_group",
                    "source": group["name"],
                    "selection": (
                        f"default/existing GROMACS index group: "
                        f"{group['name']}"
                    ),
                }
            )

        return targets


def prompt_rmsd_whole_resname_targets(
    gro: Path,
    tpr: Path,
) -> list[dict]:
    """
    Select one or multiple complete residue/molecule TYPES.
    Each resname becomes one independent RMSD target.
    """
    resnames = prompt_numbered_residue_selection(
        gro,
        title="RMSD — COMPLETE RESIDUE / MOLECULE TYPE",
        instruction=(
            "Select one or more complete residue/molecule TYPES. "
            "Each selected resname receives its own RMSD curve. "
            "If a resname has many independent molecules, the combined-group "
            "RMSD describes the whole selected group; for one-molecule RMSD "
            "use mode 4 instead."
        ),
        allow_multiple=True,
        selection_name="RMSD target type",
        tpr=tpr,
    )

    targets = []

    for resname in resnames:
        indices = atoms_for_resname(
            gro,
            resname,
        )

        if not indices:
            raise RecoverableAnalysisError(
                f"No atoms/beads found for resname '{resname}'."
            )

        targets.append(
            {
                "label": resname,
                "group": (
                    "RMSD_CALC_"
                    + sanitize_filename(
                        resname
                    )
                ),
                "indices": indices,
                "mode": "whole_resname_type",
                "source": resname,
                "selection": (
                    f"all atoms/beads belonging to resname {resname}"
                ),
            }
        )

    return targets


def prompt_rmsd_atom_subset_targets(
    gro: Path,
    tpr: Path,
) -> list[dict]:
    """
    Select one or multiple resnames and then atom/bead subsets by NUMBER.

    This works for:
      AA -> atom names
      CG -> bead names
    """
    resnames = prompt_numbered_residue_selection(
        gro,
        title="RMSD — RESNAME SELECTION FOR ATOM/BEAD SUBSETS",
        instruction=(
            "Select one or more residue/molecule TYPES. "
            "For each selected type, the next step will show its atom/bead "
            "names with numbers."
        ),
        allow_multiple=True,
        selection_name="RMSD subset resname",
        tpr=tpr,
    )

    targets = []

    for resname in resnames:
        print(
            "\n" + "=" * 78
        )
        print(
            f" RMSD ATOM/BEAD SUBSET — {resname}"
        )
        print(
            "=" * 78
        )

        names = prompt_specific_atom_bead_names_for_msd(
            gro,
            tpr,
            resname,
        )

        indices = []

        for atomname in names:
            indices.extend(
                validate_region_atom_names(
                    gro,
                    resname,
                    [atomname],
                    "RMSD_SUBSET",
                )
            )

        label = (
            f"{resname}:"
            + ",".join(
                names
            )
        )

        targets.append(
            {
                "label": label,
                "group": (
                    "RMSD_CALC_"
                    + sanitize_filename(
                        resname
                        + "_"
                        + "_".join(
                            names
                        )
                    )
                ),
                "indices": indices,
                "mode": "atom_bead_subset",
                "source": resname,
                "selection": (
                    f"resname {resname}; atom/bead names: "
                    + ", ".join(
                        names
                    )
                ),
            }
        )

    return targets


def prompt_rmsd_individual_instance_targets(
    gro: Path,
    tpr: Path,
) -> list[dict]:
    """
    Select one resname and then one or multiple INDIVIDUAL molecule/residue
    instances. Each instance gets an independent RMSD curve.
    """
    selected = prompt_numbered_residue_selection(
        gro,
        title="RMSD — INDIVIDUAL MOLECULE / RESIDUE TYPE",
        instruction=(
            "Select ONE residue/molecule TYPE. "
            "You will then choose one or more individual occurrences."
        ),
        allow_multiple=False,
        selection_name="RMSD individual type",
        tpr=tpr,
    )

    resname = selected[0]

    blocks = lipid_residue_blocks(
        gro,
        resname,
    )

    if not blocks:
        raise RecoverableAnalysisError(
            f"No individual blocks named '{resname}' were found."
        )

    _show_instance_preview(
        blocks,
        resname,
    )

    while True:
        raw = ask(
            "Select individual occurrence number(s), e.g. 1 2 5"
        )

        numbers, invalid = _parse_numbered_selection(
            raw,
            len(blocks),
        )

        if invalid or not numbers:
            if invalid:
                print(
                    "[WARNING] Invalid occurrence number(s): "
                    + ", ".join(invalid)
                )
            print(
                f"[ACTION] Enter one or more numbers from 1-{len(blocks)}."
            )
            continue

        print("\n[SELECTED INDIVIDUAL RMSD TARGETS]")

        for number in numbers:
            block = blocks[number - 1]
            print(
                f"  [{number:>5d}] {resname}; "
                f"GRO resnr={block['resnr']}; "
                f"atoms/beads={len(block['atoms'])}"
            )

        if yes_no(
            "Confirm these individual RMSD targets?",
            True,
        ):
            break

    targets = []

    for number in numbers:
        block = blocks[number - 1]

        indices = [
            idx
            for idx, _ in block["atoms"]
        ]

        label = (
            f"{resname}_instance{number}_"
            f"resnr{block['resnr']}"
        )

        targets.append(
            {
                "label": label,
                "group": (
                    "RMSD_CALC_"
                    + sanitize_filename(
                        label
                    )
                ),
                "indices": indices,
                "mode": "individual_instance",
                "source": resname,
                "selection": (
                    f"{resname} occurrence {number}; "
                    f"GRO resnr {block['resnr']}; "
                    "whole individual residue/molecule"
                ),
            }
        )

    return targets




def _unique_indices_preserve_order(
    indices: list[int],
) -> list[int]:
    """
    Remove duplicate atom/bead indices while preserving first appearance.

    This is important when users merge overlapping GROMACS groups such as
    Protein + Backbone, or any other partially overlapping selections.
    """
    seen = set()
    unique = []

    for idx in indices:
        if idx not in seen:
            seen.add(idx)
            unique.append(idx)

    return unique


def prompt_custom_rmsd_group_name(
    suggested: str = "MERGED_RMSD_GROUP",
) -> str:
    """
    Ask for a new custom group name and normalize it to a GROMACS-safe,
    filesystem-safe identifier.

    Spaces and punctuation are converted to underscores so the same name can
    safely be used in:
      - the NDX group
      - output directory
      - XVG/CSV/PNG/result/log filenames
    """
    print("\n" + "-" * 78)
    print(" CUSTOM MERGED RMSD GROUP NAME")
    print("-" * 78)
    print(
        "Enter the NEW name for the merged calculation group. "
        "Example: Protein_LIG_complex"
    )

    while True:
        raw = ask(
            "New merged RMSD group name",
            suggested,
        ).strip()

        if not raw:
            print("[WARNING] Group name cannot be empty.")
            continue

        normalized = sanitize_filename(
            raw.replace(
                " ",
                "_",
            )
        )

        if not normalized:
            print(
                "[WARNING] The entered name does not contain usable "
                "letters/numbers. Please enter another name."
            )
            continue

        if normalized[0].isdigit():
            normalized = (
                "G_"
                + normalized
            )

        if normalized != raw:
            print(
                f"[INFO] Group name normalized for GROMACS/files: "
                f"{normalized}"
            )

        if yes_no(
            f"Use merged group name '{normalized}'?",
            True,
        ):
            return normalized


def prompt_rmsd_merged_default_groups(
    default_ndx: Path,
) -> tuple[list[int], str]:
    """
    Merge multiple existing/default GROMACS index groups into one selection.
    Overlapping atom indices are automatically de-duplicated.
    """
    while True:
        groups = show_numbered_ndx_groups_for_rmsd(
            default_ndx,
            "RMSD MERGED GROUP — EXISTING GROMACS GROUPS",
        )

        print()
        print(
            "Select TWO OR MORE groups to merge. "
            "Overlapping particles are counted only once."
        )
        print(
            "Example: Protein + LIG, or several custom GROMACS index groups."
        )

        raw = ask(
            "Select group numbers to merge, e.g. 2 6"
        )

        numbers, invalid = _parse_numbered_selection(
            raw,
            len(groups),
        )

        if invalid or len(numbers) < 2:
            if invalid:
                print(
                    "[WARNING] Invalid group number(s): "
                    + ", ".join(invalid)
                )
            print(
                "[ACTION] Select at least TWO valid GROMACS groups."
            )
            continue

        selected = [
            groups[number - 1]
            for number in numbers
        ]

        merged = []

        for group in selected:
            merged.extend(
                group["indices"]
            )

        merged = _unique_indices_preserve_order(
            merged
        )

        print("\n[GROUPS TO MERGE]")
        for number in numbers:
            group = groups[number - 1]
            print(
                f"  [{number:>3d}] {group['name']:<28s} "
                f"{len(group['indices'])} atoms/beads"
            )

        print(
            f"\n  Unique atoms/beads after merge: {len(merged)}"
        )

        if yes_no(
            "Confirm this merged selection?",
            True,
        ):
            description = (
                "merged existing GROMACS groups: "
                + " + ".join(
                    group["name"]
                    for group in selected
                )
            )
            return merged, description


def prompt_rmsd_merged_resname_types(
    gro: Path,
    tpr: Path,
) -> tuple[list[int], str]:
    """
    Merge all atoms/beads from TWO OR MORE complete residue/molecule types.

    Example:
        DEF + CHOL
        LIG + PROT (when PROT is a single coarse-grained resname)
    """
    while True:
        resnames = prompt_numbered_residue_selection(
            gro,
            title="RMSD MERGED GROUP — COMPLETE RESIDUE/MOLECULE TYPES",
            instruction=(
                "Select TWO OR MORE complete residue/molecule TYPES. "
                "ALL atoms/beads of ALL occurrences of the selected resnames "
                "will be merged into ONE RMSD group."
            ),
            allow_multiple=True,
            selection_name="merged RMSD type",
            tpr=tpr,
        )

        if len(resnames) < 2:
            print(
                "[WARNING] A merged multi-type RMSD group requires at least "
                "two selected residue/molecule types."
            )
            if yes_no(
                "Select the merged residue/molecule types again?",
                True,
            ):
                continue
            raise RecoverableAnalysisError(
                "Custom merged RMSD group cancelled."
            )

        indices = []

        for resname in resnames:
            indices.extend(
                atoms_for_resname(
                    gro,
                    resname,
                )
            )

        indices = _unique_indices_preserve_order(
            indices
        )

        description = (
            "merged complete resname types: "
            + " + ".join(
                resnames
            )
        )

        print(
            f"[OK] Merged {len(resnames)} residue/molecule type(s); "
            f"{len(indices)} unique atoms/beads."
        )

        return indices, description


def prompt_rmsd_instance_numbers_for_one_resname(
    gro: Path,
    resname: str,
) -> tuple[list[int], str]:
    """
    Select all or specific individual occurrences of one resname for a
    custom merged group.
    """
    blocks = lipid_residue_blocks(
        gro,
        resname,
    )

    if not blocks:
        raise RecoverableAnalysisError(
            f"No molecule/residue blocks named '{resname}' were found."
        )

    print("\n" + "-" * 78)
    print(
        f" MERGED GROUP — INDIVIDUAL INSTANCES OF {resname}"
    )
    print("-" * 78)
    print(
        f"Detected {len(blocks)} occurrence(s)."
    )
    print()
    print(" [1] Include ALL occurrences of this resname")
    print(" [2] Include selected individual occurrences only")

    while True:
        choice = ask(
            f"Choose how {resname} enters the merged group",
            "2",
        )

        if choice in {
            "1",
            "2",
        }:
            break

        print(
            "[WARNING] Please choose 1 or 2."
        )

    if choice == "1":
        numbers = list(
            range(
                1,
                len(blocks) + 1,
            )
        )
    else:
        _show_instance_preview(
            blocks,
            resname,
        )

        while True:
            raw = ask(
                f"Select {resname} occurrence number(s), e.g. 1 2 5"
            )

            numbers, invalid = _parse_numbered_selection(
                raw,
                len(blocks),
            )

            if invalid or not numbers:
                if invalid:
                    print(
                        "[WARNING] Invalid occurrence number(s): "
                        + ", ".join(invalid)
                    )
                print(
                    f"[ACTION] Enter one or more numbers from "
                    f"1-{len(blocks)}."
                )
                continue

            break

    indices = []
    labels = []

    for number in numbers:
        block = blocks[
            number - 1
        ]

        indices.extend(
            idx
            for idx, _ in block["atoms"]
        )

        labels.append(
            f"{resname}#{number}(resnr={block['resnr']})"
        )

    return (
        indices,
        ", ".join(
            labels
        ),
    )


def prompt_rmsd_merged_individual_instances(
    gro: Path,
    tpr: Path,
) -> tuple[list[int], str]:
    """
    Merge multiple selected molecule/residue INSTANCES, including instances
    from different resnames.

    Example:
        DEF occurrence #1 + DEF occurrence #2
        Protein-resname occurrence #1 + LIG occurrence #1
        selected lipids from several residue types
    """
    resnames = prompt_numbered_residue_selection(
        gro,
        title="RMSD MERGED GROUP — RESNAME TYPES CONTAINING INSTANCES",
        instruction=(
            "Select one or more residue/molecule TYPES. "
            "For each selected type, you will then choose which individual "
            "occurrences enter the ONE merged RMSD group."
        ),
        allow_multiple=True,
        selection_name="merged-instance type",
        tpr=tpr,
    )

    all_indices = []
    descriptions = []
    total_instances = 0

    for resname in resnames:
        indices, description = (
            prompt_rmsd_instance_numbers_for_one_resname(
                gro,
                resname,
            )
        )

        all_indices.extend(
            indices
        )
        descriptions.append(
            description
        )

        # description is one comma-separated item per selected instance.
        total_instances += (
            description.count(
                ", "
            )
            + 1
            if description
            else 0
        )

    all_indices = _unique_indices_preserve_order(
        all_indices
    )

    if total_instances < 2:
        print(
            "[WARNING] Only one individual molecule/residue instance was "
            "selected. The merged-group mode is intended for >=2 instances."
        )

        if not yes_no(
            "Continue anyway?",
            False,
        ):
            raise RecoverableAnalysisError(
                "Custom merged RMSD instance group cancelled."
            )

    description = (
        "merged individual molecule/residue instances: "
        + " | ".join(
            descriptions
        )
    )

    return all_indices, description


def prompt_rmsd_custom_merged_target(
    gro: Path,
    tpr: Path,
    default_ndx: Path,
) -> list[dict]:
    """
    NEW RMSD TARGET MODE:
    create ONE user-named group by merging multiple molecular selections.

    Merge sources:
      [1] multiple existing GROMACS index groups
      [2] multiple complete residue/molecule types
      [3] multiple selected individual molecule/residue instances

    The returned target is one RMSD group and therefore produces one RMSD
    time series for the entire merged selection.
    """
    print("\n" + "=" * 78)
    print(" RMSD — CUSTOM MERGED GROUP")
    print("=" * 78)
    print(
        "Create ONE RMSD calculation group from multiple molecular selections, "
        "then rename the merged group."
    )
    print()
    print(" [1] Merge existing/default GROMACS index groups")
    print(
        "     Example: Protein + LIG or several existing custom groups."
    )
    print()
    print(" [2] Merge complete residue/molecule TYPES")
    print(
        "     Example: DEF + CHOL; all occurrences of each selected resname "
        "are included."
    )
    print()
    print(" [3] Merge selected INDIVIDUAL molecule/residue instances")
    print(
        "     Example: DEF #1 + DEF #2 + LIG #1."
    )

    while True:
        mode = ask(
            "Choose merged-group source",
            "2",
        )

        if mode in {
            "1",
            "2",
            "3",
        }:
            break

        print(
            "[WARNING] Please choose 1, 2, or 3."
        )

    if mode == "1":
        indices, selection = prompt_rmsd_merged_default_groups(
            default_ndx,
        )
        source_mode = "merged_gromacs_groups"
        suggested = "Merged_GROMACS_Group"

    elif mode == "2":
        indices, selection = prompt_rmsd_merged_resname_types(
            gro,
            tpr,
        )
        source_mode = "merged_resname_types"
        suggested = "Merged_Molecules"

    else:
        indices, selection = prompt_rmsd_merged_individual_instances(
            gro,
            tpr,
        )
        source_mode = "merged_individual_instances"
        suggested = "Merged_Instances"

    if not indices:
        raise RecoverableAnalysisError(
            "The custom merged RMSD group is empty."
        )

    custom_name = prompt_custom_rmsd_group_name(
        suggested,
    )

    print("\n[CUSTOM MERGED RMSD GROUP]")
    print(
        f"  New group name       : {custom_name}"
    )
    print(
        f"  Unique atoms/beads   : {len(indices)}"
    )
    print(
        f"  Selection definition : {selection}"
    )
    print()
    print(
        "[NOTE] This produces ONE RMSD curve for the entire merged selection. "
        "It is not an arithmetic average of separate molecule RMSD curves."
    )

    if not yes_no(
        "Confirm this custom merged RMSD group?",
        True,
    ):
        return prompt_rmsd_custom_merged_target(
            gro,
            tpr,
            default_ndx,
        )

    return [
        {
            "label": custom_name,
            "group": custom_name,
            "indices": indices,
            "mode": source_mode,
            "source": custom_name,
            "selection": selection,
            "custom_group_name": custom_name,
        }
    ]


def prompt_rmsd_target_mode(
    gro: Path,
    tpr: Path,
    default_ndx: Path,
) -> list[dict]:
    """
    High-level, user-friendly RMSD target selector.

    Original modes 1-4 are preserved.
    Mode 5 adds one user-named merged RMSD group.
    """
    print("\n" + "=" * 78)
    print(" RMSD — CALCULATION TARGET")
    print("=" * 78)
    print()
    print(" [1] GROMACS/default index group(s)")
    print(
        "     Best for Protein, Backbone, C-alpha, System, etc."
    )
    print()
    print(" [2] Complete residue/molecule TYPE(S)")
    print(
        "     Each selected type is calculated SEPARATELY."
    )
    print()
    print(" [3] Specific atom/bead subset(s) inside selected resname(s)")
    print(
        "     AA: choose atoms by number; CG: choose beads by number."
    )
    print()
    print(" [4] Specific individual molecule/residue INSTANCE(S)")
    print(
        "     Each selected individual instance is calculated SEPARATELY."
    )
    print()
    print(" [5] CUSTOM MERGED GROUP — multiple molecules/groups -> ONE RMSD")
    print(
        "     Merge multiple GROMACS groups, complete resname types, or "
        "selected individual molecules, then rename the merged group."
    )
    print("=" * 78)

    while True:
        mode = ask(
            "Choose RMSD target mode",
            "1",
        )

        if mode == "1":
            return prompt_rmsd_default_index_targets(
                default_ndx,
            )

        if mode == "2":
            return prompt_rmsd_whole_resname_targets(
                gro,
                tpr,
            )

        if mode == "3":
            return prompt_rmsd_atom_subset_targets(
                gro,
                tpr,
            )

        if mode == "4":
            return prompt_rmsd_individual_instance_targets(
                gro,
                tpr,
            )

        if mode == "5":
            return prompt_rmsd_custom_merged_target(
                gro,
                tpr,
                default_ndx,
            )

        print(
            "[WARNING] Please choose 1, 2, 3, 4, or 5."
        )


def prompt_rmsd_common_fit_group(
    gro: Path,
    tpr: Path,
    default_ndx: Path,
) -> dict:
    """
    Build ONE common fitting group for all RMSD targets.

    This supports common workflows such as:
      - fit Backbone, calculate whole Protein RMSD
      - fit one CG protein backbone-like bead set, calculate another subset
      - fit one molecular scaffold, calculate a flexible substituent group
    """
    print("\n" + "=" * 78)
    print(" RMSD — COMMON FITTING GROUP")
    print("=" * 78)
    print()
    print(" [1] One existing/default GROMACS index group")
    print(" [2] One or more complete residue/molecule types combined")
    print(" [3] Specific atom/bead subset inside ONE selected resname")
    print(" [4] One specific individual molecule/residue instance")
    print(" [5] Custom merged group (multiple groups/molecules, renamed)")

    while True:
        mode = ask(
            "Choose fitting-group source",
            "1",
        )

        if mode in {
            "1",
            "2",
            "3",
            "4",
            "5",
        }:
            break

        print(
            "[WARNING] Please choose 1, 2, 3, 4, or 5."
        )

    if mode == "1":
        groups = show_numbered_ndx_groups_for_rmsd(
            default_ndx,
            "RMSD — SELECT ONE COMMON FITTING GROUP",
        )

        while True:
            raw = ask(
                "Select ONE fitting-group number"
            )

            numbers, invalid = _parse_numbered_selection(
                raw,
                len(groups),
            )

            if (
                invalid
                or len(numbers) != 1
            ):
                print(
                    "[WARNING] Select exactly one valid group number."
                )
                continue

            group = groups[
                numbers[0] - 1
            ]

            if yes_no(
                f"Use '{group['name']}' as the common fitting group?",
                True,
            ):
                return {
                    "label": group["name"],
                    "indices": list(
                        group["indices"]
                    ),
                    "selection": (
                        f"default GROMACS index group: "
                        f"{group['name']}"
                    ),
                }

    if mode == "2":
        resnames = prompt_numbered_residue_selection(
            gro,
            title="RMSD FIT — COMPLETE RESIDUE/MOLECULE TYPES",
            instruction=(
                "Select one or more complete residue/molecule TYPES. "
                "They will be combined into ONE fitting group."
            ),
            allow_multiple=True,
            selection_name="RMSD fit type",
            tpr=tpr,
        )

        indices = []

        for resname in resnames:
            indices.extend(
                atoms_for_resname(
                    gro,
                    resname,
                )
            )

        return {
            "label": "+".join(
                resnames
            ),
            "indices": indices,
            "selection": (
                "combined complete resnames: "
                + ", ".join(
                    resnames
                )
            ),
        }

    if mode == "3":
        selected = prompt_numbered_residue_selection(
            gro,
            title="RMSD FIT — RESNAME FOR ATOM/BEAD SUBSET",
            instruction=(
                "Select ONE residue/molecule type, then choose its "
                "fit atoms/beads by number."
            ),
            allow_multiple=False,
            selection_name="RMSD fit resname",
            tpr=tpr,
        )

        resname = selected[0]

        names = prompt_specific_atom_bead_names_for_msd(
            gro,
            tpr,
            resname,
        )

        indices = []

        for atomname in names:
            indices.extend(
                validate_region_atom_names(
                    gro,
                    resname,
                    [atomname],
                    "RMSD_FIT",
                )
            )

        return {
            "label": (
                resname
                + ":"
                + ",".join(
                    names
                )
            ),
            "indices": indices,
            "selection": (
                f"resname {resname}; fit atom/bead names: "
                + ", ".join(
                    names
                )
            ),
        }

    if mode == "5":
        merged_target = prompt_rmsd_custom_merged_target(
            gro,
            tpr,
            default_ndx,
        )[0]

        return {
            "label": merged_target["label"],
            "indices": merged_target["indices"],
            "selection": (
                f"custom merged fitting group "
                f"'{merged_target['label']}': "
                f"{merged_target['selection']}"
            ),
            "group_name": merged_target["group"],
        }

    # mode 4
    selected = prompt_numbered_residue_selection(
        gro,
        title="RMSD FIT — INDIVIDUAL MOLECULE / RESIDUE TYPE",
        instruction=(
            "Select ONE residue/molecule type, then choose ONE individual "
            "occurrence for fitting."
        ),
        allow_multiple=False,
        selection_name="RMSD fit individual type",
        tpr=tpr,
    )

    resname = selected[0]

    blocks = lipid_residue_blocks(
        gro,
        resname,
    )

    _show_instance_preview(
        blocks,
        resname,
    )

    while True:
        raw = ask(
            "Select ONE individual occurrence number"
        )

        numbers, invalid = _parse_numbered_selection(
            raw,
            len(blocks),
        )

        if (
            invalid
            or len(numbers) != 1
        ):
            print(
                "[WARNING] Select exactly one valid occurrence number."
            )
            continue

        number = numbers[0]
        block = blocks[
            number - 1
        ]

        if yes_no(
            (
                f"Use {resname} occurrence #{number} "
                f"(GRO resnr {block['resnr']}) as fitting group?"
            ),
            True,
        ):
            return {
                "label": (
                    f"{resname}_instance{number}_"
                    f"resnr{block['resnr']}"
                ),
                "indices": [
                    idx
                    for idx, _ in block["atoms"]
                ],
                "selection": (
                    f"{resname} occurrence {number}; "
                    f"GRO resnr {block['resnr']}"
                ),
            }


def prompt_rmsd_settings() -> dict:
    """
    Ask RMSD settings with simple defaults that work for AA and CG.

    GROMACS default structural RMSD normally uses:
      - fit = rot+trans
      - mass weighting = yes
      - PBC check = yes
    """
    begin_ps, end_ps, dt_ps = prompt_time_window(
        require_dt=True,
    )

    print("\n" + "-" * 78)
    print(" RMSD FITTING METHOD")
    print("-" * 78)
    print(" [1] Rotation + translation  (recommended structural RMSD)")
    print(" [2] Translation only")
    print(" [3] No fitting")

    while True:
        choice = ask(
            "Choose RMSD fitting method",
            "1",
        )

        if choice == "1":
            fit = "rot+trans"
            break

        if choice == "2":
            fit = "translation"
            break

        if choice == "3":
            fit = "none"
            break

        print(
            "[WARNING] Please choose 1, 2, or 3."
        )

    print("\n" + "-" * 78)
    print(" RMSD WEIGHTING / PBC")
    print("-" * 78)
    print(
        "Mass-weighted RMSD is the GROMACS default and is valid when the TPR "
        "contains correct masses for AA atoms or CG beads."
    )

    mass_weighted = yes_no(
        "Use mass-weighted fitting/RMSD (-mw yes)?",
        True,
    )

    pbc_check = yes_no(
        "Use GROMACS PBC check (-pbc yes)?",
        True,
    )

    return {
        "begin_ps": begin_ps,
        "end_ps": end_ps,
        "dt_ps": dt_ps,
        "fit": fit,
        "mass_weighted": mass_weighted,
        "pbc_check": pbc_check,
    }


def prompt_rmsd_reference_structure(
    tpr: Path,
) -> Path:
    """
    Default to the selected TPR as the reference structure.

    An alternative reference is supported for advanced users, but atom
    ordering/selection compatibility must match the trajectory/index.
    """
    print("\n" + "-" * 78)
    print(" RMSD REFERENCE STRUCTURE")
    print("-" * 78)
    print(
        f"Default reference: {tpr}"
    )
    print(
        "This normally corresponds to the starting structure stored in the TPR."
    )

    if yes_no(
        "Use the selected TPR as the RMSD reference structure?",
        True,
    ):
        return tpr

    print(
        "[NOTE] An alternative reference must have compatible atom ordering "
        "for the selected RMSD/fit groups."
    )

    return ask_existing_path(
        "Alternative reference structure (.tpr/.gro/.pdb/etc.)"
    )


def write_rmsd_index(
    ndx_file: Path,
    fit_indices: list[int],
    calc_indices: list[int],
    fit_group_name: str = "RMSD_FIT",
    calc_group_name: str = "RMSD_CALC",
):
    """
    Write a minimal two-group RMSD index.

    The RMSD calculation group keeps the target's actual/custom name.
    Therefore a user-created merged group such as:
        Protein_LIG_complex
    appears literally in the generated NDX as:
        [ Protein_LIG_complex ]

    The fitting group remains separately named so the two interactive
    selections are deterministic even when fitting and calculation use the
    same particle indices.
    """
    if not fit_indices:
        raise RecoverableAnalysisError(
            "RMSD fitting group is empty."
        )

    if not calc_indices:
        raise RecoverableAnalysisError(
            "RMSD calculation group is empty."
        )

    fit_group_name = sanitize_filename(
        fit_group_name.replace(" ", "_")
    )
    calc_group_name = sanitize_filename(
        calc_group_name.replace(" ", "_")
    )

    if not fit_group_name:
        fit_group_name = "RMSD_FIT"

    if not calc_group_name:
        calc_group_name = "RMSD_CALC"

    if fit_group_name == calc_group_name:
        fit_group_name = (
            fit_group_name
            + "_FIT"
        )

    ndx_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with ndx_file.open(
        "w",
        encoding="utf-8",
    ) as fh:
        write_ndx_group(
            fh,
            fit_group_name,
            _unique_indices_preserve_order(
                fit_indices
            ),
        )

        write_ndx_group(
            fh,
            calc_group_name,
            _unique_indices_preserve_order(
                calc_indices
            ),
        )

    print(
        f"[OK] RMSD index generated: {ndx_file}"
    )
    print(
        f"     Fit group  : [ {fit_group_name} ] "
        f"({len(_unique_indices_preserve_order(fit_indices))} atoms/beads)"
    )
    print(
        f"     RMSD group : [ {calc_group_name} ] "
        f"({len(_unique_indices_preserve_order(calc_indices))} atoms/beads)"
    )

    return fit_group_name, calc_group_name


def run_rmsd_command(
    xtc: Path,
    reference: Path,
    ndx_file: Path,
    raw_xvg: Path,
    log_file: Path,
    settings: dict,
    fit_group_name: str = "RMSD_FIT",
    calc_group_name: str = "RMSD_CALC",
):
    """
    Run one gmx rms calculation.

    Compatible with the classic GROMACS 2020.x interactive interface:
      first selection  = fitting group
      second selection = RMSD calculation group
    """
    cmd = [
        GMX,
        "rms",
        "-s",
        str(reference),
        "-f",
        str(xtc),
        "-n",
        str(ndx_file),
        "-o",
        str(raw_xvg),
        "-what",
        "rmsd",
        "-fit",
        settings["fit"],
        "-mw",
        (
            "yes"
            if settings["mass_weighted"]
            else "no"
        ),
        "-pbc",
        (
            "yes"
            if settings["pbc_check"]
            else "no"
        ),
        *gmx_time_args(
            settings["begin_ps"],
            settings["end_ps"],
            settings["dt_ps"],
        ),
    ]

    result = _animated_process(
        cmd,
        input_text=(
            fit_group_name
            + "\n"
            + calc_group_name
            + "\n"
        ),
        label=f"RMSD: {raw_xvg.stem}",
    )

    log_file.write_text(
        log_timestamp_header() + "COMMAND\n"
        + " ".join(
            map(
                str,
                cmd,
            )
        )
        + "\n\nFIT GROUP\n"
        + fit_group_name
        + "\n"
        + "\nRMSD CALCULATION GROUP\n"
        + calc_group_name
        + "\n"
        + "\nSTDOUT\n"
        + (
            result.stdout
            or ""
        )
        + "\nSTDERR\n"
        + (
            result.stderr
            or ""
        ),
        encoding="utf-8",
    )

    if result.returncode != 0:
        if result.stdout:
            print(
                result.stdout
            )
        if result.stderr:
            print(
                result.stderr,
                file=sys.stderr,
            )

        raise RecoverableAnalysisError(
            f"gmx rms failed. Full log: {log_file}"
        )

    if (
        not raw_xvg.exists()
        or raw_xvg.stat().st_size == 0
    ):
        raise RecoverableAnalysisError(
            f"gmx rms did not generate {raw_xvg}"
        )

    return log_file


def convert_rmsd_raw_to_outputs(
    raw_xvg: Path,
    out_xvg: Path,
    csv_file: Path,
    target: dict,
    fit_description: str,
    reference: Path,
    settings: dict,
) -> list[tuple[float, float]]:
    """
    Convert GROMACS RMSD output time from ps to ns and create a clean CSV/XVG.
    """
    rows = parse_xvg(
        raw_xvg,
        min_cols=2,
    )

    if not rows:
        raise RecoverableAnalysisError(
            f"No RMSD data could be read from {raw_xvg}."
        )

    pairs = [
        (
            row[0] / 1000.0,
            row[1],
        )
        for row in rows
    ]

    with out_xvg.open(
        "w",
        encoding="utf-8",
    ) as fh:
        fh.write(
            "# RMSD calculated with GROMACS gmx rms\n"
        )
        fh.write(
            f"# Target: {target['label']}\n"
        )
        fh.write(
            f"# Target selection: {target['selection']}\n"
        )
        fh.write(
            f"# Fitting group: {fit_description}\n"
        )
        fh.write(
            f"# Reference: {reference}\n"
        )
        fh.write(
            f"# Fit method: {settings['fit']}\n"
        )
        fh.write(
            f"# Mass weighted: {settings['mass_weighted']}\n"
        )
        fh.write(
            f"# PBC check: {settings['pbc_check']}\n"
        )
        fh.write(
            '@ title "RMSD vs Time"\n'
        )
        fh.write(
            '@ xaxis label "Time (ns)"\n'
        )
        fh.write(
            '@ yaxis label "RMSD (nm)"\n'
        )
        fh.write(
            '@ s0 legend "RMSD"\n'
        )

        for time_ns, rmsd_nm in pairs:
            fh.write(
                f"{time_ns:.10f} "
                f"{rmsd_nm:.10f}\n"
            )

    with csv_file.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        writer = csv.writer(
            fh
        )

        writer.writerow(
            [
                "Time_ns",
                "RMSD_nm",
            ]
        )

        for time_ns, rmsd_nm in pairs:
            writer.writerow(
                [
                    f"{time_ns:.10f}",
                    f"{rmsd_nm:.10f}",
                ]
            )

    return pairs


def write_rmsd_result_file(
    result_file: Path,
    target: dict,
    fit_description: str,
    fit_count: int,
    reference: Path,
    settings: dict,
    pairs: list[tuple[float, float]],
    ndx_file: Path,
    xvg_file: Path,
    csv_file: Path,
    log_file: Path,
):
    """
    Write a readable summary for one RMSD target.
    """
    values = [
        value
        for _, value in pairs
    ]

    avg = mean(
        values
    )

    sd = (
        stdev(
            values
        )
        if len(
            values
        ) > 1
        else 0.0
    )

    lines = [
        "=" * 78,
        "GROMACS RMSD RESULT",
        "=" * 78,
        f"Target mode            : {target['mode']}",
        f"RMSD target            : {target['label']}",
        f"RMSD target selection  : {target['selection']}",
        f"RMSD particle count    : {len(target['indices'])}",
        "",
        f"Fitting group          : {fit_description}",
        f"Fitting particle count : {fit_count}",
        f"Reference structure    : {reference}",
        f"Fit method             : {settings['fit']}",
        f"Mass weighted          : {settings['mass_weighted']}",
        f"PBC check              : {settings['pbc_check']}",
        "",
        f"Trajectory begin       : {settings['begin_ps']} ps",
        f"Trajectory end         : "
        + (
            f"{settings['end_ps']} ps"
            if settings["end_ps"] is not None
            else "trajectory end"
        ),
        f"Sampling interval      : {settings['dt_ps']} ps",
        "",
        f"RMSD data points       : {len(values)}",
        f"Mean RMSD              : {avg:.10f} nm",
        f"SD                     : {sd:.10f} nm",
        f"Mean +/- SD            : {avg:.10f} +/- {sd:.10f} nm",
        f"Minimum RMSD           : {min(values):.10f} nm",
        f"Maximum RMSD           : {max(values):.10f} nm",
        "",
        f"Index file             : {ndx_file}",
        f"RMSD curve             : {xvg_file}",
        f"RMSD CSV               : {csv_file}",
        f"Command log            : {log_file}",
        "=" * 78,
        "",
    ]

    result_file.write_text(
        "\n".join(
            lines
        ),
        encoding="utf-8",
    )

    return {
        "mean": avg,
        "sd": sd,
        "min": min(
            values
        ),
        "max": max(
            values
        ),
        "points": len(
            values
        ),
    }


def append_rmsd_summary_excel(
    target: dict,
    fit_description: str,
    fit_count: int,
    reference: Path,
    settings: dict,
    stats: dict,
    ndx_file: Path,
    xvg_file: Path,
    csv_file: Path,
    result_file: Path,
    log_file: Path,
):
    """
    Persistent Excel summary across RMSD calculations.
    """
    summary = (
        RMSD_DIR
        / "RMSD_summary.xlsx"
    )

    cache = (
        RMSD_DIR
        / ".RMSD_summary_data.json"
    )

    headers = [
        "target_mode",
        "target",
        "target_selection",
        "target_particle_count",
        "fit_group",
        "fit_particle_count",
        "reference_structure",
        "fit_method",
        "mass_weighted",
        "pbc_check",
        "begin_ps",
        "end_ps",
        "sampling_interval_ps",
        "mean_RMSD_nm",
        "SD_nm",
        "minimum_RMSD_nm",
        "maximum_RMSD_nm",
        "data_points",
        "index_file",
        "RMSD_xvg",
        "RMSD_csv",
        "result_file",
        "log_file",
    ]

    row = [
        target["mode"],
        target["label"],
        target["selection"],
        len(
            target["indices"]
        ),
        fit_description,
        fit_count,
        str(
            reference
        ),
        settings["fit"],
        str(
            settings["mass_weighted"]
        ),
        str(
            settings["pbc_check"]
        ),
        settings["begin_ps"],
        (
            settings["end_ps"]
            if settings["end_ps"] is not None
            else "trajectory_end"
        ),
        settings["dt_ps"],
        stats["mean"],
        stats["sd"],
        stats["min"],
        stats["max"],
        stats["points"],
        str(
            ndx_file
        ),
        str(
            xvg_file
        ),
        str(
            csv_file
        ),
        str(
            result_file
        ),
        str(
            log_file
        ),
    ]

    append_excel_summary(
        summary,
        cache,
        "RMSD Summary",
        headers,
        row,
    )


def make_rmsd_plot(
    xvg_file: Path,
    png_file: Path,
    title: str,
):
    """
    Generate a publication-friendly RMSD time-series PNG when matplotlib
    is available.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "[INFO] matplotlib is not installed; RMSD PNG plotting skipped. "
            "The XVG/CSV files remain available."
        )
        return

    rows = parse_xvg(
        xvg_file,
        min_cols=2,
    )

    if not rows:
        return

    x = [
        row[0]
        for row in rows
    ]
    y = [
        row[1]
        for row in rows
    ]

    fig, ax = plt.subplots(
        figsize=(
            7.2,
            5.0,
        )
    )

    ax.plot(
        x,
        y,
        linewidth=1.2,
    )

    ax.set_title(
        title
    )
    ax.set_xlabel(
        "Time (ns)"
    )
    ax.set_ylabel(
        "RMSD (nm)"
    )

    fig.tight_layout()

    fig.savefig(
        png_file,
        dpi=300,
    )

    plt.close(
        fig
    )

    print(
        f"[OK] RMSD plot generated: {png_file}"
    )


def run_rmsd_batch(
    gro: Path,
    tpr: Path,
    xtc: Path,
    reference: Path,
    default_ndx: Path,
    targets: list[dict],
    settings: dict,
):
    """
    Run all selected RMSD targets sequentially.

    The user chooses once whether:
      A) each target is fitted on itself, or
      B) one common fit group is used for every target.
    """
    RMSD_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("\n" + "=" * 78)
    print(" RMSD FIT GROUP")
    print("=" * 78)
    print()
    print(
        "GROMACS can use the SAME particles for fitting and RMSD, "
        "or use a DIFFERENT fitting group."
    )
    print()
    print(
        "Example: fit Protein Backbone, then calculate RMSD of the whole Protein."
    )

    same_fit = yes_no(
        "Use each RMSD target itself as its fitting group?",
        True,
    )

    common_fit = None

    if not same_fit:
        common_fit = prompt_rmsd_common_fit_group(
            gro,
            tpr,
            default_ndx,
        )

    definition_file = (
        RMSD_DIR
        / "RMSD_target_definitions.xlsx"
    )

    definition_headers = [
        "target",
        "target_mode",
        "target_selection",
        "target_particle_count",
        "RMSD_index_group_name",
        "fit_mode",
        "fit_group",
        "fit_particle_count",
    ]

    definition_rows = []

    for target in targets:
        if same_fit:
            fit_description = (
                "same as RMSD target"
            )
            fit_count = len(
                target["indices"]
            )
        else:
            fit_description = (
                common_fit["selection"]
            )
            fit_count = len(
                common_fit["indices"]
            )

        definition_rows.append(
            [
                target["label"],
                target["mode"],
                target["selection"],
                len(
                    target["indices"]
                ),
                target.get(
                    "group",
                    target["label"],
                ),
                (
                    "same_as_target"
                    if same_fit
                    else "common_separate_fit_group"
                ),
                fit_description,
                fit_count,
            ]
        )

    write_excel_table(
        definition_file,
        "RMSD Targets",
        definition_headers,
        definition_rows,
    )

    print("\n" + "=" * 78)
    print(
        f" STARTING RMSD BATCH: {len(targets)} TARGET(S)"
    )
    print("=" * 78)

    results = []

    for number, target in enumerate(
        targets,
        start=1,
    ):
        print(
            f"\n---------------- RMSD {number}/{len(targets)}: "
            f"{target['label']} ----------------"
        )

        if same_fit:
            fit_indices = target[
                "indices"
            ]
            fit_description = (
                "same as RMSD target: "
                + target["selection"]
            )
        else:
            fit_indices = common_fit[
                "indices"
            ]
            fit_description = common_fit[
                "selection"
            ]

        safe = sanitize_filename(
            target["label"]
        )

        target_dir = (
            RMSD_DIR
            / safe
        )

        target_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        ndx_file = (
            target_dir
            / f"RMSD_{safe}_index.ndx"
        )

        raw_xvg = (
            target_dir
            / f"RMSD_{safe}_raw_ps.xvg"
        )

        xvg_file = (
            target_dir
            / f"RMSD_{safe}_vs_time.xvg"
        )

        csv_file = (
            target_dir
            / f"RMSD_{safe}_vs_time.csv"
        )

        png_file = (
            target_dir
            / f"RMSD_{safe}_vs_time.png"
        )

        result_file = (
            target_dir
            / f"RMSD_{safe}_result.txt"
        )

        log_file = timestamped_log_path(
            target_dir
            / f"RMSD_{safe}.log"
        )

        fit_group_name, calc_group_name = write_rmsd_index(
            ndx_file,
            fit_indices,
            target["indices"],
            fit_group_name="RMSD_FIT",
            calc_group_name=target.get(
                "group",
                target["label"],
            ),
        )

        run_rmsd_command(
            xtc,
            reference,
            ndx_file,
            raw_xvg,
            log_file,
            settings,
            fit_group_name=fit_group_name,
            calc_group_name=calc_group_name,
        )

        pairs = convert_rmsd_raw_to_outputs(
            raw_xvg,
            xvg_file,
            csv_file,
            target,
            fit_description,
            reference,
            settings,
        )

        stats = write_rmsd_result_file(
            result_file,
            target,
            fit_description,
            len(
                fit_indices
            ),
            reference,
            settings,
            pairs,
            ndx_file,
            xvg_file,
            csv_file,
            log_file,
        )

        append_rmsd_summary_excel(
            target,
            fit_description,
            len(
                fit_indices
            ),
            reference,
            settings,
            stats,
            ndx_file,
            xvg_file,
            csv_file,
            result_file,
            log_file,
        )

        make_rmsd_plot(
            xvg_file,
            png_file,
            (
                "RMSD — "
                + target["label"]
            ),
        )

        results.append(
            {
                "target": target[
                    "label"
                ],
                "mean": stats[
                    "mean"
                ],
                "sd": stats[
                    "sd"
                ],
                "xvg": xvg_file,
                "csv": csv_file,
                "png": png_file,
                "result": result_file,
                "log": log_file,
            }
        )

        print(
            f"[OK] RMSD target completed: {target['label']}"
        )
        print(
            f"     Mean +/- SD = "
            f"{stats['mean']:.6f} +/- "
            f"{stats['sd']:.6f} nm"
        )

    print("\n" + "=" * 78)
    print(" RMSD BATCH SUMMARY")
    print("=" * 78)

    for row in results:
        print(
            f" {row['target']:<28.28s}: "
            f"{row['mean']:.6f} +/- "
            f"{row['sd']:.6f} nm"
        )

    print()
    print(
        f"  Target definitions : {definition_file}"
    )
    print(
        f"  Global Excel summary: "
        f"{RMSD_DIR / 'RMSD_summary.xlsx'}"
    )

    return results


def run_rmsd_menu(
    state: dict,
) -> bool:
    """
    Interactive RMSD controller.

    Workflow matches the style of the existing script:
      - automatic file detection
      - numbered selections
      - input validation
      - optional default settings
      - independent output folders
      - completion menu
      - retry/recovery handled by run_parameter_with_recovery()
    """
    print("\n" + "=" * 78)
    print(" RMSD — ROOT-MEAN-SQUARE DEVIATION")
    print("=" * 78)
    print(
        "Method: gmx rms compares trajectory structures with a reference "
        "structure after optional least-squares fitting."
    )
    print()
    print(
        "The fitting group and RMSD calculation group can be the same or "
        "different."
    )
    print(
        "Multiple molecules/groups can also be merged into ONE custom-named "
        "RMSD group."
    )
    print()
    print(
        "Compatible with all-atom atoms and coarse-grained beads because "
        "the script reads actual GRO/TPR names instead of hard-coding a "
        "specific force field."
    )

    files = resolve_analysis_files(
        [
            (
                "gro",
                "gro",
                "GRO file",
            ),
            (
                "tpr",
                "tpr",
                "TPR file",
            ),
            (
                "xtc",
                "xtc",
                "XTC trajectory",
            ),
        ],
        state,
    )

    gro = files[
        "gro"
    ]
    tpr = files[
        "tpr"
    ]
    xtc = files[
        "xtc"
    ]

    while True:
        default_ndx = build_rmsd_default_index(
            tpr,
        )

        targets = prompt_rmsd_target_mode(
            gro,
            tpr,
            default_ndx,
        )

        reference = prompt_rmsd_reference_structure(
            tpr,
        )

        settings = prompt_rmsd_settings()

        run_rmsd_batch(
            gro,
            tpr,
            xtc,
            reference,
            default_ndx,
            targets,
            settings,
        )

        action = post_calculation_menu(
            "RMSD",
            "Calculate another RMSD target / fitting definition",
        )

        if action == "menu":
            return False

        if action == "exit":
            return True

        reuse = yes_no(
            "Reuse the current GRO/TPR/XTC files for another RMSD calculation?",
            True,
        )

        if reuse:
            continue

        files = resolve_analysis_files(
            [
                (
                    "gro",
                    "gro",
                    "GRO file",
                ),
                (
                    "tpr",
                    "tpr",
                    "TPR file",
                ),
                (
                    "xtc",
                    "xtc",
                    "XTC trajectory",
                ),
            ],
            state,
        )

        gro = files[
            "gro"
        ]
        tpr = files[
            "tpr"
        ]
        xtc = files[
            "xtc"
        ]


# ======================================================================
# 9) MEMBRANE PERMEATION MOLECULE COUNT
# ======================================================================
#
# Purpose
# -------
# Quantify, as a function of simulation time:
#   - number of permeant molecules inside the membrane slab
#   - number of permeant molecules below the lower leaflet
#   - optionally both
#
# Membrane boundaries are NOT fixed user-entered Z values.
#
# Workflow
# --------
# 1) User selects all membrane lipid residue types.
# 2) For each lipid, user selects ONE boundary/headgroup marker atom/bead
#    by NUMBER (e.g. CHOL-ROH, HSPC-PO4).
# 3) The membrane is centered for every frame using gmx trjconv.
# 4) Marker Z positions are extracted for every frame.
# 5) A robust 1D two-cluster split identifies lower and upper leaflet marker
#    planes independently in every frame:
#        Z_lower(t), Z_upper(t)
# 6) Permeant molecule COM Z positions are extracted with gmx trajectory.
# 7) Classification:
#        membrane:
#            Z_lower(t) <= Z_COM(t) <= Z_upper(t)
#        below lower leaflet:
#            Z_COM(t) < Z_lower(t)
#
# A user-selected outward buffer may expand the membrane slab:
#        effective lower = Z_lower - buffer
#        effective upper = Z_upper + buffer
#
# Outputs include time series, mean +/- SD, membrane Z-range time series,
# instantaneous below-membrane occupancy, and cumulative first-passage count.
# ======================================================================


def _permeation_all_atom_indices(
    gro: Path,
) -> list[int]:
    indices = []

    for block in read_gro_residue_blocks(
        gro
    ):
        indices.extend(
            idx
            for idx, _ in block["atoms"]
        )

    return indices


def prompt_permeation_boundary_markers(
    gro: Path,
    tpr: Path,
    lipid_names: list[str],
) -> dict[str, str]:
    """
    Select ONE representative outer/headgroup marker per membrane lipid type.

    Numbered selection makes this easy:
        [1] CHOL-ROH [2] CHOL-C1 ...
    """
    markers = {}

    print("\n" + "=" * 78)
    print(" PERMEATION — MEMBRANE BOUNDARY MARKERS")
    print("=" * 78)
    print(
        "Choose ONE representative outer/headgroup atom or bead for EACH "
        "membrane lipid species."
    )
    print(
        "These markers are used to locate the lower and upper leaflet planes "
        "independently in every trajectory frame."
    )

    for number, resname in enumerate(
        lipid_names,
        start=1,
    ):
        print(
            f"\nLipid {number}/{len(lipid_names)}: {resname}"
        )

        marker = prompt_numbered_atom_bead_selection(
            gro,
            resname,
            title=f"PERMEATION BOUNDARY MARKER — {resname}",
            instruction=(
                "Select ONE headgroup/interface marker near the aqueous "
                "boundary of this lipid. Examples: CHOL-ROH or a phospholipid "
                "phosphate/phosphorus marker when appropriate."
            ),
            allow_multiple=False,
            min_count=1,
            tpr=tpr,
            require_once_per_residue=True,
        )[0]

        markers[resname] = marker

    print("\n[MEMBRANE BOUNDARY MARKERS]")
    for resname in lipid_names:
        print(
            f"  {resname:<14s}: {resname}-{markers[resname]}"
        )

    return markers


def prompt_permeation_target(
    gro: Path,
    tpr: Path,
) -> tuple[str, str]:
    """
    Select the permeating molecule/residue type and COM grouping method.

    whole_res_com:
        best when one permeant molecule corresponds to one residue block.

    whole_mol_com:
        uses TPR molecular grouping and is useful when a molecule contains
        multiple residues but the selected resname identifies those molecules.
    """
    selected = prompt_numbered_residue_selection(
        gro,
        title="PERMEATION — TARGET MOLECULE SELECTION",
        instruction=(
            "Select ONE permeating molecule/residue TYPE whose molecule count "
            "will be tracked relative to the membrane."
        ),
        allow_multiple=False,
        selection_name="permeating molecule",
        tpr=tpr,
    )

    resname = selected[0]

    print("\nTarget position definition")
    print(" [1] Residue COM — recommended when one molecule = one residue")
    print(" [2] Whole-molecule COM from TPR — for multi-residue molecules")

    while True:
        choice = ask(
            "Choose target COM definition",
            "1",
        )

        if choice == "1":
            return resname, "whole_res_com"

        if choice == "2":
            return resname, "whole_mol_com"

        print(
            "[WARNING] Please choose 1 or 2."
        )


def prompt_permeation_count_mode() -> str:
    print("\n" + "-" * 78)
    print(" PERMEATION — COUNT OUTPUT")
    print("-" * 78)
    print(" [1] Molecules INSIDE the membrane only")
    print(" [2] Molecules BELOW the lower leaflet only")
    print(" [3] BOTH membrane-inside and below-lower counts")

    while True:
        choice = ask(
            "Choose molecule-count output",
            "3",
        )

        if choice == "1":
            return "membrane"

        if choice == "2":
            return "below"

        if choice == "3":
            return "both"

        print(
            "[WARNING] Please choose 1, 2, or 3."
        )


def write_permeation_index(
    gro: Path,
    lipid_names: list[str],
    markers: dict[str, str],
    target_resname: str,
    out_file: Path,
):
    """
    Build one clean index containing:
      PERM_SYSTEM
      PERM_MEMBRANE
      PERM_MARKERS
      PERM_TARGET
    """
    all_indices = _permeation_all_atom_indices(
        gro
    )

    membrane_indices = []

    for resname in lipid_names:
        membrane_indices.extend(
            atoms_for_resname(
                gro,
                resname,
            )
        )

    membrane_indices = list(
        dict.fromkeys(
            membrane_indices
        )
    )

    marker_indices = []

    for resname in lipid_names:
        marker_indices.extend(
            validate_region_atom_names(
                gro,
                resname,
                [markers[resname]],
                "PERMEATION_BOUNDARY_MARKER",
            )
        )

    marker_indices = list(
        dict.fromkeys(
            marker_indices
        )
    )

    target_indices = atoms_for_resname(
        gro,
        target_resname,
    )

    if not all_indices:
        raise RecoverableAnalysisError(
            "Could not build PERM_SYSTEM group."
        )

    if len(marker_indices) < 2:
        raise RecoverableAnalysisError(
            "Too few membrane boundary markers were selected."
        )

    if not target_indices:
        raise RecoverableAnalysisError(
            f"No atoms/beads found for permeating target '{target_resname}'."
        )

    out_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with out_file.open(
        "w",
        encoding="utf-8",
    ) as fh:
        write_ndx_group(
            fh,
            "PERM_SYSTEM",
            all_indices,
        )
        write_ndx_group(
            fh,
            "PERM_MEMBRANE",
            membrane_indices,
        )
        write_ndx_group(
            fh,
            "PERM_MARKERS",
            marker_indices,
        )
        write_ndx_group(
            fh,
            "PERM_TARGET",
            target_indices,
        )

    print(
        f"[OK] Permeation index generated: {out_file}"
    )
    print(
        f"     PERM_MEMBRANE : {len(membrane_indices)} atoms/beads"
    )
    print(
        f"     PERM_MARKERS  : {len(marker_indices)} marker atoms/beads"
    )
    print(
        f"     PERM_TARGET   : {len(target_indices)} atoms/beads"
    )

    return {
        "system": all_indices,
        "membrane": membrane_indices,
        "markers": marker_indices,
        "target": target_indices,
    }


def prepare_permeation_centered_trajectory(
    xtc: Path,
    tpr: Path,
    ndx_file: Path,
    out_file: Path,
    begin_ps: float,
    end_ps: float | None,
    dt_ps: float,
):
    """
    Center the selected membrane in every frame and keep molecules whole.

    This removes whole-membrane Z drift before membrane-range detection.
    """
    cmd = [
        GMX,
        "trjconv",
        "-f",
        str(xtc),
        "-s",
        str(tpr),
        "-n",
        str(ndx_file),
        "-o",
        str(out_file),
        "-pbc",
        "mol",
        "-center",
        "-ur",
        "compact",
        *gmx_time_args(
            begin_ps,
            end_ps,
            dt_ps,
        ),
    ]

    result = _animated_process(
        cmd,
        input_text=(
            "PERM_MEMBRANE\n"
            "PERM_SYSTEM\n"
        ),
        label="Permeation: centering membrane trajectory",
    )

    log_file = timestamped_log_path(
        out_file.with_suffix(
            ".log"
        )
    )

    log_file.write_text(
        log_timestamp_header() + "COMMAND\n"
        + " ".join(
            map(
                str,
                cmd,
            )
        )
        + "\n\nCENTER GROUP\nPERM_MEMBRANE\n"
        + "\nOUTPUT GROUP\nPERM_SYSTEM\n"
        + "\nSTDOUT\n"
        + (
            result.stdout
            or ""
        )
        + "\nSTDERR\n"
        + (
            result.stderr
            or ""
        ),
        encoding="utf-8",
    )

    if result.returncode != 0:
        raise RecoverableAnalysisError(
            f"Membrane-centering trajectory step failed. See {log_file}"
        )

    if (
        not out_file.exists()
        or out_file.stat().st_size == 0
    ):
        raise RecoverableAnalysisError(
            "Centered permeation trajectory was not generated."
        )

    return log_file


def run_permeation_z_trajectory(
    centered_xtc: Path,
    tpr: Path,
    ndx_file: Path,
    selection: str,
    out_file: Path,
    label: str,
):
    """
    Extract only Z coordinates using gmx trajectory.

    For marker atoms:
        group "PERM_MARKERS"

    For target molecule COMs:
        whole_res_com of group "PERM_TARGET"
        or
        whole_mol_com of group "PERM_TARGET"
    """
    cmd = [
        GMX,
        "trajectory",
        "-f",
        str(centered_xtc),
        "-s",
        str(tpr),
        "-n",
        str(ndx_file),
        "-ox",
        str(out_file),
        "-xvg",
        "none",
        "-nox",
        "-noy",
        "-z",
        "-select",
        selection,
    ]

    result = _animated_process(
        cmd,
        label=label,
    )

    log_file = timestamped_log_path(
        out_file.with_suffix(
            ".log"
        )
    )

    log_file.write_text(
        log_timestamp_header() + "COMMAND\n"
        + " ".join(
            map(
                str,
                cmd,
            )
        )
        + "\n\nSELECTION\n"
        + selection
        + "\n\nSTDOUT\n"
        + (
            result.stdout
            or ""
        )
        + "\nSTDERR\n"
        + (
            result.stderr
            or ""
        ),
        encoding="utf-8",
    )

    if result.returncode != 0:
        raise RecoverableAnalysisError(
            f"{label} failed. See {log_file}"
        )

    rows = parse_xvg(
        out_file,
        min_cols=2,
    )

    if not rows:
        raise RecoverableAnalysisError(
            f"No Z-coordinate data were generated for: {selection}"
        )

    return rows, log_file


def split_permeation_leaflet_planes(
    z_values: list[float],
) -> tuple[float, float]:
    """
    Robust 1D k-means (k=2) followed by median leaflet positions.

    The two clusters correspond to lower- and upper-leaflet boundary markers.
    Median positions reduce sensitivity to local undulations or a few outliers.
    """
    values = [
        float(z)
        for z in z_values
        if math.isfinite(
            float(z)
        )
    ]

    if len(values) < 2:
        raise RecoverableAnalysisError(
            "Too few valid membrane marker Z coordinates in one frame."
        )

    c1 = min(
        values
    )
    c2 = max(
        values
    )

    if abs(
        c2 - c1
    ) < 1e-8:
        raise RecoverableAnalysisError(
            "Upper/lower leaflet markers cannot be separated in one frame."
        )

    g1 = []
    g2 = []

    for _ in range(
        100
    ):
        g1 = []
        g2 = []

        for z in values:
            if abs(
                z - c1
            ) <= abs(
                z - c2
            ):
                g1.append(
                    z
                )
            else:
                g2.append(
                    z
                )

        if not g1 or not g2:
            raise RecoverableAnalysisError(
                "Membrane marker clustering produced an empty leaflet."
            )

        n1 = mean(
            g1
        )
        n2 = mean(
            g2
        )

        if (
            abs(
                n1 - c1
            ) < 1e-8
            and abs(
                n2 - c2
            ) < 1e-8
        ):
            break

        c1 = n1
        c2 = n2

    lower_cluster = (
        g1
        if mean(
            g1
        ) < mean(
            g2
        )
        else g2
    )

    upper_cluster = (
        g2
        if lower_cluster is g1
        else g1
    )

    lower = median(
        lower_cluster
    )
    upper = median(
        upper_cluster
    )

    if upper <= lower:
        raise RecoverableAnalysisError(
            "Invalid leaflet marker ordering."
        )

    return lower, upper


def align_permeation_marker_target_rows(
    marker_rows: list[list[float]],
    target_rows: list[list[float]],
) -> list[tuple[float, list[float], list[float]]]:
    """
    Align marker and target-COM rows by rounded trajectory time.
    """
    markers = {
        round(
            row[0],
            6,
        ): row[1:]
        for row in marker_rows
        if len(
            row
        ) >= 2
    }

    targets = {
        round(
            row[0],
            6,
        ): row[1:]
        for row in target_rows
        if len(
            row
        ) >= 2
    }

    common = sorted(
        set(
            markers
        )
        & set(
            targets
        )
    )

    if len(
        common
    ) < 2:
        raise RecoverableAnalysisError(
            "Too few common frames between membrane markers and target COMs."
        )

    return [
        (
            time_ps,
            markers[time_ps],
            targets[time_ps],
        )
        for time_ps in common
    ]


def calculate_permeation_count_series(
    aligned_rows,
    boundary_buffer_nm: float,
    initial_loading_molecules: int,
):
    """
    Classify every target-molecule COM relative to the dynamic membrane slab
    and calculate time-dependent permeation percentages.

    Operational permeation rate (%) requested by the workflow:
        cumulative unique molecules that newly crossed below the lower leaflet
        ---------------------------------------------------------------------- x 100
                       user-entered initial loading molecules

    The script ALSO retains instantaneous molecule counts because they are the
    raw data used to derive the percentages.

    Percentage fields:
      membrane_fraction_percent
          instantaneous molecules currently inside the membrane / N0 * 100

      below_fraction_percent
          instantaneous molecules currently below the lower leaflet / N0 * 100

      cumulative_permeation_percent
          cumulative unique newly crossed molecules / N0 * 100

    The final value of cumulative_permeation_percent is reported as the
    FINAL PERMEATION RATE (%) for this operational definition.
    """
    if initial_loading_molecules <= 0:
        raise RecoverableAnalysisError(
            "Initial loading molecule count must be greater than zero."
        )

    records = []

    initially_below = None
    ever_crossed = None

    expected_targets = None

    for frame_no, (
        time_ps,
        marker_z,
        target_z,
    ) in enumerate(
        aligned_rows
    ):
        lower_raw, upper_raw = split_permeation_leaflet_planes(
            marker_z
        )

        lower = (
            lower_raw
            - boundary_buffer_nm
        )
        upper = (
            upper_raw
            + boundary_buffer_nm
        )

        if expected_targets is None:
            expected_targets = len(
                target_z
            )

            initially_below = [
                z < lower
                for z in target_z
            ]

            ever_crossed = [
                False
                for _ in target_z
            ]

        elif len(
            target_z
        ) != expected_targets:
            raise RecoverableAnalysisError(
                "The number of target COM positions changed between frames. "
                "Use a static target selection."
            )

        inside = []
        below = []

        for z in target_z:
            inside.append(
                lower <= z <= upper
            )
            below.append(
                z < lower
            )

        for i, is_below in enumerate(
            below
        ):
            if (
                not initially_below[i]
                and is_below
            ):
                ever_crossed[i] = True

        membrane_count = sum(
            inside
        )
        below_count = sum(
            below
        )
        cumulative_crossed = sum(
            ever_crossed
        )

        records.append(
            {
                "time_ps": time_ps,
                "time_ns": time_ps / 1000.0,
                "lower_raw_nm": lower_raw,
                "upper_raw_nm": upper_raw,
                "lower_nm": lower,
                "upper_nm": upper,
                "thickness_nm": upper - lower,
                "membrane_count": membrane_count,
                "below_count": below_count,
                "above_count": (
                    expected_targets
                    - membrane_count
                    - below_count
                ),
                "cumulative_crossed_to_lower": cumulative_crossed,
                "initial_loading_molecules": initial_loading_molecules,
                "membrane_fraction_percent": (
                    membrane_count
                    / initial_loading_molecules
                    * 100.0
                ),
                "below_fraction_percent": (
                    below_count
                    / initial_loading_molecules
                    * 100.0
                ),
                "cumulative_permeation_percent": (
                    cumulative_crossed
                    / initial_loading_molecules
                    * 100.0
                ),
                "total_target_molecules": expected_targets,
            }
        )

    return records



def _permeation_mean_sd(
    values: list[float],
) -> tuple[float, float]:
    if not values:
        return float(
            "nan"
        ), float(
            "nan"
        )

    avg = mean(
        values
    )

    sd = (
        stdev(
            values
        )
        if len(
            values
        ) > 1
        else 0.0
    )

    return avg, sd


def write_permeation_outputs(
    records: list[dict],
    count_mode: str,
    target_resname: str,
    lipid_names: list[str],
    markers: dict[str, str],
    target_com_mode: str,
    boundary_buffer_nm: float,
    initial_loading_molecules: int,
    outdir: Path,
):
    """
    Write molecule-count raw data AND derived permeation-rate outputs.

    Final permeation rate (%) =
        final cumulative unique newly crossed molecules / initial loading * 100

    Outputs:
      - count vs time CSV/XVG/PNG
      - permeation rate vs time CSV/XVG/PNG
      - dynamic membrane Z boundaries
      - TXT result including FINAL PERMEATION RATE (%)
      - timestamped summary LOG including FINAL PERMEATION RATE (%)
      - Excel summary
    """
    if not records:
        raise RecoverableAnalysisError(
            "No permeation records are available for output."
        )

    csv_file = (
        outdir
        / "permeation_molecule_count_vs_time.csv"
    )

    xvg_file = (
        outdir
        / "permeation_molecule_count_vs_time.xvg"
    )

    rate_csv = (
        outdir
        / "permeation_rate_vs_time.csv"
    )

    rate_xvg = (
        outdir
        / "permeation_rate_vs_time.xvg"
    )

    boundary_csv = (
        outdir
        / "membrane_z_range_vs_time.csv"
    )

    stats_file = (
        outdir
        / "permeation_rate_result.txt"
    )

    summary_xlsx = (
        outdir
        / "permeation_rate_summary.xlsx"
    )

    summary_log = timestamped_log_path(
        outdir
        / "permeation_rate_summary.log"
    )

    # --------------------------------------------------------------
    # Count + rate combined raw table
    # --------------------------------------------------------------
    with csv_file.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        writer = csv.writer(
            fh
        )

        writer.writerow(
            [
                "Time_ns",
                "Lower_membrane_Z_nm",
                "Upper_membrane_Z_nm",
                "Dynamic_membrane_thickness_nm",
                "Molecules_inside_membrane",
                "Molecules_below_lower_leaflet",
                "Molecules_above_upper_leaflet",
                "Cumulative_unique_crossed_to_lower",
                "Initial_loading_molecules",
                "Inside_membrane_percent_of_initial",
                "Below_lower_percent_of_initial",
                "Cumulative_permeation_rate_percent",
                "Detected_target_molecules",
            ]
        )

        for row in records:
            writer.writerow(
                [
                    f"{row['time_ns']:.10f}",
                    f"{row['lower_nm']:.10f}",
                    f"{row['upper_nm']:.10f}",
                    f"{row['thickness_nm']:.10f}",
                    row["membrane_count"],
                    row["below_count"],
                    row["above_count"],
                    row["cumulative_crossed_to_lower"],
                    row["initial_loading_molecules"],
                    f"{row['membrane_fraction_percent']:.10f}",
                    f"{row['below_fraction_percent']:.10f}",
                    f"{row['cumulative_permeation_percent']:.10f}",
                    row["total_target_molecules"],
                ]
            )

    # --------------------------------------------------------------
    # Dedicated permeability percentage table
    # --------------------------------------------------------------
    with rate_csv.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        writer = csv.writer(
            fh
        )

        writer.writerow(
            [
                "Time_ns",
                "Cumulative_unique_permeated_molecules",
                "Initial_loading_molecules",
                "Permeation_rate_percent",
                "Instantaneous_below_lower_percent",
                "Instantaneous_inside_membrane_percent",
            ]
        )

        for row in records:
            writer.writerow(
                [
                    f"{row['time_ns']:.10f}",
                    row["cumulative_crossed_to_lower"],
                    initial_loading_molecules,
                    f"{row['cumulative_permeation_percent']:.10f}",
                    f"{row['below_fraction_percent']:.10f}",
                    f"{row['membrane_fraction_percent']:.10f}",
                ]
            )

    # --------------------------------------------------------------
    # Membrane dynamic Z range
    # --------------------------------------------------------------
    with boundary_csv.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        writer = csv.writer(
            fh
        )

        writer.writerow(
            [
                "Time_ns",
                "Lower_marker_plane_raw_nm",
                "Upper_marker_plane_raw_nm",
                "Effective_lower_boundary_nm",
                "Effective_upper_boundary_nm",
                "Effective_membrane_Z_range",
                "Thickness_nm",
            ]
        )

        for row in records:
            writer.writerow(
                [
                    f"{row['time_ns']:.10f}",
                    f"{row['lower_raw_nm']:.10f}",
                    f"{row['upper_raw_nm']:.10f}",
                    f"{row['lower_nm']:.10f}",
                    f"{row['upper_nm']:.10f}",
                    (
                        f"[{row['lower_nm']:.6f}, "
                        f"{row['upper_nm']:.6f}]"
                    ),
                    f"{row['thickness_nm']:.10f}",
                ]
            )

    # --------------------------------------------------------------
    # Count XVG
    # --------------------------------------------------------------
    with xvg_file.open(
        "w",
        encoding="utf-8",
    ) as fh:
        fh.write(
            "# Membrane permeation molecule count vs time\n"
        )
        fh.write(
            f"# Target resname: {target_resname}\n"
        )
        fh.write(
            f"# Initial loading molecules: {initial_loading_molecules}\n"
        )
        fh.write(
            f"# Target position: {target_com_mode}\n"
        )
        fh.write(
            f"# Boundary buffer: {boundary_buffer_nm} nm\n"
        )
        fh.write(
            '@ title "Membrane permeation molecule count vs Time"\n'
        )
        fh.write(
            '@ xaxis label "Time (ns)"\n'
        )
        fh.write(
            '@ yaxis label "Molecule count"\n'
        )

        legends = []

        if count_mode in {
            "membrane",
            "both",
        }:
            legends.append(
                "Inside membrane"
            )

        if count_mode in {
            "below",
            "both",
        }:
            legends.append(
                "Below lower leaflet"
            )

        legends.append(
            "Cumulative unique crossed to lower"
        )

        for i, legend in enumerate(
            legends
        ):
            fh.write(
                f'@ s{i} legend "{legend}"\n'
            )

        for row in records:
            values = [
                f"{row['time_ns']:.10f}"
            ]

            if count_mode in {
                "membrane",
                "both",
            }:
                values.append(
                    str(
                        row["membrane_count"]
                    )
                )

            if count_mode in {
                "below",
                "both",
            }:
                values.append(
                    str(
                        row["below_count"]
                    )
                )

            values.append(
                str(
                    row["cumulative_crossed_to_lower"]
                )
            )

            fh.write(
                " ".join(
                    values
                )
                + "\n"
            )

    # --------------------------------------------------------------
    # Permeation-rate XVG
    # --------------------------------------------------------------
    with rate_xvg.open(
        "w",
        encoding="utf-8",
    ) as fh:
        fh.write(
            "# Operational permeation rate (%) vs time\n"
        )
        fh.write(
            "# Definition: cumulative unique newly crossed molecules / "
            "initial loading molecules * 100\n"
        )
        fh.write(
            f"# Initial loading molecules N0: {initial_loading_molecules}\n"
        )
        fh.write(
            '@ title "Membrane permeation rate vs Time"\n'
        )
        fh.write(
            '@ xaxis label "Time (ns)"\n'
        )
        fh.write(
            '@ yaxis label "Permeation rate (%)"\n'
        )
        fh.write(
            '@ s0 legend "Cumulative permeation rate"\n'
        )
        fh.write(
            '@ s1 legend "Below lower leaflet / initial"\n'
        )
        fh.write(
            '@ s2 legend "Inside membrane / initial"\n'
        )

        for row in records:
            fh.write(
                f"{row['time_ns']:.10f} "
                f"{row['cumulative_permeation_percent']:.10f} "
                f"{row['below_fraction_percent']:.10f} "
                f"{row['membrane_fraction_percent']:.10f}\n"
            )

    # --------------------------------------------------------------
    # Statistics
    # --------------------------------------------------------------
    membrane_values = [
        row["membrane_count"]
        for row in records
    ]

    below_values = [
        row["below_count"]
        for row in records
    ]

    membrane_percent_values = [
        row["membrane_fraction_percent"]
        for row in records
    ]

    below_percent_values = [
        row["below_fraction_percent"]
        for row in records
    ]

    cumulative_rate_values = [
        row["cumulative_permeation_percent"]
        for row in records
    ]

    lower_values = [
        row["lower_nm"]
        for row in records
    ]

    upper_values = [
        row["upper_nm"]
        for row in records
    ]

    thickness_values = [
        row["thickness_nm"]
        for row in records
    ]

    mem_avg, mem_sd = _permeation_mean_sd(
        membrane_values
    )

    below_avg, below_sd = _permeation_mean_sd(
        below_values
    )

    mem_pct_avg, mem_pct_sd = _permeation_mean_sd(
        membrane_percent_values
    )

    below_pct_avg, below_pct_sd = _permeation_mean_sd(
        below_percent_values
    )

    cumulative_pct_avg, cumulative_pct_sd = _permeation_mean_sd(
        cumulative_rate_values
    )

    lower_avg, lower_sd = _permeation_mean_sd(
        lower_values
    )

    upper_avg, upper_sd = _permeation_mean_sd(
        upper_values
    )

    thick_avg, thick_sd = _permeation_mean_sd(
        thickness_values
    )

    final_cumulative = records[-1][
        "cumulative_crossed_to_lower"
    ]

    final_permeation_rate = records[-1][
        "cumulative_permeation_percent"
    ]

    final_below_percent = records[-1][
        "below_fraction_percent"
    ]

    detected_target_molecules = records[0][
        "total_target_molecules"
    ]

    marker_text = ", ".join(
        f"{resname}-{markers[resname]}"
        for resname in lipid_names
    )

    mismatch_note = ""

    if detected_target_molecules != initial_loading_molecules:
        mismatch_note = (
            "\nWARNING: user-entered initial loading does not equal the "
            "number of target molecule COMs detected in the trajectory.\n"
            f"  Initial loading input : {initial_loading_molecules}\n"
            f"  Detected target COMs : {detected_target_molecules}\n"
            "The requested user-entered initial loading is retained as the "
            "denominator for the reported permeation rate.\n"
        )

    stats_text = (
        "=" * 78
        + "\nMEMBRANE PERMEATION RATE RESULT\n"
        + "=" * 78
        + "\n"
        + f"Result generated at           : {log_timestamp_text()}\n"
        + f"Target molecule/resname       : {target_resname}\n"
        + f"Target COM definition         : {target_com_mode}\n"
        + f"Initial loading molecules N0  : {initial_loading_molecules}\n"
        + f"Detected target molecules     : {detected_target_molecules}\n"
        + f"Membrane lipid types          : {', '.join(lipid_names)}\n"
        + f"Boundary markers              : {marker_text}\n"
        + f"Boundary buffer               : {boundary_buffer_nm:.6f} nm\n"
        + "\n"
        + "PERMEATION DEFINITION\n"
        + "  Permeation rate (%) = cumulative unique molecules that newly\n"
        + "  crossed below the lower leaflet / initial loading N0 * 100\n"
        + "\n"
        + f"Final permeated molecules     : {final_cumulative} molecules\n"
        + f"FINAL PERMEATION RATE         : {final_permeation_rate:.6f} %\n"
        + f"Final below-lower fraction    : {final_below_percent:.6f} %\n"
        + "\n"
        + f"Time-averaged cumulative rate : "
        + f"{cumulative_pct_avg:.6f} +/- {cumulative_pct_sd:.6f} %\n"
        + f"Inside-membrane count         : "
        + f"{mem_avg:.6f} +/- {mem_sd:.6f} molecules\n"
        + f"Inside-membrane / initial     : "
        + f"{mem_pct_avg:.6f} +/- {mem_pct_sd:.6f} %\n"
        + f"Below-lower-leaflet count     : "
        + f"{below_avg:.6f} +/- {below_sd:.6f} molecules\n"
        + f"Below-lower / initial         : "
        + f"{below_pct_avg:.6f} +/- {below_pct_sd:.6f} %\n"
        + "\n"
        + f"Lower membrane boundary       : "
        + f"{lower_avg:.6f} +/- {lower_sd:.6f} nm\n"
        + f"Upper membrane boundary       : "
        + f"{upper_avg:.6f} +/- {upper_sd:.6f} nm\n"
        + f"Dynamic membrane thickness    : "
        + f"{thick_avg:.6f} +/- {thick_sd:.6f} nm\n"
        + mismatch_note
        + "=" * 78
        + "\n"
    )

    stats_file.write_text(
        stats_text,
        encoding="utf-8",
    )

    # Timestamped LOG required by the user.
    summary_log.write_text(
        log_timestamp_header()
        + stats_text,
        encoding="utf-8",
    )

    # --------------------------------------------------------------
    # Excel summary
    # --------------------------------------------------------------
    headers = [
        "result_generated_at",
        "target_resname",
        "target_COM_definition",
        "initial_loading_molecules",
        "detected_target_molecules",
        "final_permeated_unique_molecules",
        "FINAL_PERMEATION_RATE_PERCENT",
        "final_below_lower_percent",
        "time_average_cumulative_permeation_percent",
        "time_SD_cumulative_permeation_percent",
        "membrane_lipids",
        "boundary_markers",
        "boundary_buffer_nm",
        "inside_membrane_mean_count",
        "inside_membrane_SD",
        "inside_membrane_mean_percent_initial",
        "inside_membrane_SD_percent_initial",
        "below_lower_mean_count",
        "below_lower_SD",
        "below_lower_mean_percent_initial",
        "below_lower_SD_percent_initial",
        "lower_boundary_mean_nm",
        "lower_boundary_SD_nm",
        "upper_boundary_mean_nm",
        "upper_boundary_SD_nm",
        "membrane_thickness_mean_nm",
        "membrane_thickness_SD_nm",
        "time_points",
    ]

    rows = [
        [
            log_timestamp_text(),
            target_resname,
            target_com_mode,
            initial_loading_molecules,
            detected_target_molecules,
            final_cumulative,
            final_permeation_rate,
            final_below_percent,
            cumulative_pct_avg,
            cumulative_pct_sd,
            ", ".join(
                lipid_names
            ),
            marker_text,
            boundary_buffer_nm,
            mem_avg,
            mem_sd,
            mem_pct_avg,
            mem_pct_sd,
            below_avg,
            below_sd,
            below_pct_avg,
            below_pct_sd,
            lower_avg,
            lower_sd,
            upper_avg,
            upper_sd,
            thick_avg,
            thick_sd,
            len(
                records
            ),
        ]
    ]

    write_excel_table(
        summary_xlsx,
        "Permeation Rate",
        headers,
        rows,
    )

    # --------------------------------------------------------------
    # Plots
    # --------------------------------------------------------------
    try:
        import matplotlib.pyplot as plt

        times = [
            row["time_ns"]
            for row in records
        ]

        # Count plot
        fig, ax = plt.subplots(
            figsize=(
                7.2,
                5.0,
            )
        )

        if count_mode in {
            "membrane",
            "both",
        }:
            ax.plot(
                times,
                membrane_values,
                linewidth=1.2,
                label="Inside membrane",
            )

        if count_mode in {
            "below",
            "both",
        }:
            ax.plot(
                times,
                below_values,
                linewidth=1.2,
                label="Below lower leaflet",
            )

        ax.plot(
            times,
            [
                row["cumulative_crossed_to_lower"]
                for row in records
            ],
            linewidth=1.2,
            label="Cumulative unique crossed",
        )

        ax.set_xlabel(
            "Time (ns)"
        )
        ax.set_ylabel(
            "Molecule count"
        )
        ax.set_title(
            "Membrane permeation molecule count"
        )
        ax.legend()
        fig.tight_layout()
        fig.savefig(
            outdir
            / "permeation_molecule_count_vs_time.png",
            dpi=300,
        )
        plt.close(
            fig
        )

        # Permeation rate plot
        fig, ax = plt.subplots(
            figsize=(
                7.2,
                5.0,
            )
        )

        ax.plot(
            times,
            cumulative_rate_values,
            linewidth=1.2,
            label="Cumulative permeation rate",
        )

        ax.set_xlabel(
            "Time (ns)"
        )
        ax.set_ylabel(
            "Permeation rate (%)"
        )
        ax.set_title(
            "Membrane permeation rate"
        )
        ax.legend()
        fig.tight_layout()
        fig.savefig(
            outdir
            / "permeation_rate_vs_time.png",
            dpi=300,
        )
        plt.close(
            fig
        )

        # Dynamic membrane Z-range plot
        fig, ax = plt.subplots(
            figsize=(
                7.2,
                5.0,
            )
        )

        ax.plot(
            times,
            lower_values,
            linewidth=1.2,
            label="Lower membrane boundary",
        )

        ax.plot(
            times,
            upper_values,
            linewidth=1.2,
            label="Upper membrane boundary",
        )

        ax.set_xlabel(
            "Time (ns)"
        )
        ax.set_ylabel(
            "Z (nm)"
        )
        ax.set_title(
            "Dynamic membrane Z range"
        )
        ax.legend()
        fig.tight_layout()
        fig.savefig(
            outdir
            / "membrane_z_range_vs_time.png",
            dpi=300,
        )
        plt.close(
            fig
        )

    except ImportError:
        print(
            "[INFO] matplotlib is not installed; permeation PNG plots skipped."
        )

    # Explicit final screen output requested by the user.
    print(
        stats_text
    )

    print("\n" + "*" * 78)
    print(
        f" FINAL PERMEATION RATE = {final_permeation_rate:.6f} %"
    )
    print(
        f" Final permeated molecules = {final_cumulative} / "
        f"{initial_loading_molecules}"
    )
    print("*" * 78)

    return {
        "count_csv": csv_file,
        "count_xvg": xvg_file,
        "rate_csv": rate_csv,
        "rate_xvg": rate_xvg,
        "boundary_csv": boundary_csv,
        "stats": stats_file,
        "summary_log": summary_log,
        "summary_xlsx": summary_xlsx,
        "final_permeation_rate_percent": final_permeation_rate,
        "final_permeated_molecules": final_cumulative,
    }



def run_permeation_count_menu(
    state: dict,
) -> bool:
    """
    Interactive membrane-permeation molecule-count + permeation-rate workflow.
    """
    print("\n" + "=" * 78)
    print(" MEMBRANE PERMEATION RATE (%)")
    print("=" * 78)
    print(
        "This module automatically tracks the membrane Z range in every "
        "trajectory frame, counts target molecules, and converts cumulative "
        "transmembrane events into a user-defined initial-dose percentage."
    )
    print()
    print(
        "Operational final permeation rate (%) = "
        "cumulative unique newly crossed molecules / initial loading * 100."
    )
    print()
    print(
        "The membrane itself is centered before analysis, so whole-membrane "
        "translation does not create false permeation counts."
    )

    files = resolve_analysis_files(
        [
            (
                "gro",
                "gro",
                "GRO file",
            ),
            (
                "tpr",
                "tpr",
                "TPR file",
            ),
            (
                "xtc",
                "xtc",
                "XTC trajectory",
            ),
        ],
        state,
    )

    gro = files[
        "gro"
    ]
    tpr = files[
        "tpr"
    ]
    xtc = files[
        "xtc"
    ]

    while True:
        begin_ps, end_ps, dt_ps = prompt_time_window(
            require_dt=True,
        )

        lipid_names = prompt_membrane_lipid_names(
            gro,
            tpr,
        )

        markers = prompt_permeation_boundary_markers(
            gro,
            tpr,
            lipid_names,
        )

        target_resname, target_com_mode = prompt_permeation_target(
            gro,
            tpr,
        )

        # User explicitly enters the initial loading amount.
        gro_counts = residue_counts_from_gro(
            gro
        )

        detected_residue_blocks = int(
            gro_counts.get(
                target_resname,
                0,
            )
        )

        print("\n" + "-" * 78)
        print(" INITIAL MOLECULE LOADING")
        print("-" * 78)
        print(
            "Enter the initial target-molecule loading N0 used as the "
            "denominator of the permeation rate."
        )

        if detected_residue_blocks > 0:
            print(
                f"[INFO] GRO currently contains {detected_residue_blocks} "
                f"residue block(s) named {target_resname}. "
                "This is shown only as a convenient reference."
            )

        initial_loading_molecules = ask_positive_integer(
            "Initial loading molecule count N0",
            (
                detected_residue_blocks
                if detected_residue_blocks > 0
                else None
            ),
        )

        print(
            f"[CONFIRMED] Initial loading N0 = "
            f"{initial_loading_molecules} molecule(s)"
        )

        count_mode = prompt_permeation_count_mode()

        boundary_buffer_nm = ask_float(
            "Optional outward membrane-boundary buffer (nm; 0 = marker planes)",
            0.0,
        )

        PERMEATION_DIR.mkdir(
            parents=True,
            exist_ok=True,
        )

        # Timestamped subfolder avoids overwriting repeated analyses.
        run_stamp = log_timestamp()

        run_tag = sanitize_filename(
            target_resname
        )

        subdir = (
            PERMEATION_DIR
            / f"{run_tag}_{run_stamp}"
        )

        subdir.mkdir(
            parents=True,
            exist_ok=True,
        )

        ndx_file = (
            subdir
            / "permeation_analysis_index.ndx"
        )

        write_permeation_index(
            gro,
            lipid_names,
            markers,
            target_resname,
            ndx_file,
        )

        centered_xtc = (
            subdir
            / "permeation_membrane_centered.xtc"
        )

        center_log = prepare_permeation_centered_trajectory(
            xtc,
            tpr,
            ndx_file,
            centered_xtc,
            begin_ps,
            end_ps,
            dt_ps,
        )

        marker_xvg = (
            subdir
            / "membrane_boundary_marker_Z_raw.xvg"
        )

        marker_rows, marker_log = run_permeation_z_trajectory(
            centered_xtc,
            tpr,
            ndx_file,
            'group "PERM_MARKERS"',
            marker_xvg,
            "Permeation: membrane-marker Z positions",
        )

        target_xvg = (
            subdir
            / "permeant_COM_Z_raw.xvg"
        )

        target_selection = (
            f'{target_com_mode} of group "PERM_TARGET"'
        )

        target_rows, target_log = run_permeation_z_trajectory(
            centered_xtc,
            tpr,
            ndx_file,
            target_selection,
            target_xvg,
            "Permeation: target molecule COM Z positions",
        )

        aligned = align_permeation_marker_target_rows(
            marker_rows,
            target_rows,
        )

        records = calculate_permeation_count_series(
            aligned,
            boundary_buffer_nm,
            initial_loading_molecules,
        )

        actual_detected = records[0][
            "total_target_molecules"
        ]

        if actual_detected != initial_loading_molecules:
            print("\n[WARNING]")
            print(
                f"User-entered initial loading N0 = "
                f"{initial_loading_molecules}"
            )
            print(
                f"Trajectory target COMs detected = "
                f"{actual_detected}"
            )
            print(
                "The user-entered N0 will be retained as the denominator, "
                "as requested."
            )

        outputs = write_permeation_outputs(
            records,
            count_mode,
            target_resname,
            lipid_names,
            markers,
            target_com_mode,
            boundary_buffer_nm,
            initial_loading_molecules,
            subdir,
        )

        print("\n[PERMEATION OUTPUTS]")
        print(
            f"  Count/raw table    : {outputs['count_csv']}"
        )
        print(
            f"  Count XVG          : {outputs['count_xvg']}"
        )
        print(
            f"  Permeation-rate CSV: {outputs['rate_csv']}"
        )
        print(
            f"  Permeation-rate XVG: {outputs['rate_xvg']}"
        )
        print(
            f"  Membrane Z range   : {outputs['boundary_csv']}"
        )
        print(
            f"  Result TXT         : {outputs['stats']}"
        )
        print(
            f"  Timestamped LOG    : {outputs['summary_log']}"
        )
        print(
            f"  Excel summary      : {outputs['summary_xlsx']}"
        )
        print(
            f"  Centering log      : {center_log}"
        )
        print(
            f"  Marker log         : {marker_log}"
        )
        print(
            f"  Target COM log     : {target_log}"
        )

        print("\n" + "#" * 78)
        print(
            f" FINAL PERMEATION RATE (%) : "
            f"{outputs['final_permeation_rate_percent']:.6f} %"
        )
        print(
            f" FINAL PERMEATED MOLECULES : "
            f"{outputs['final_permeated_molecules']} / "
            f"{initial_loading_molecules}"
        )
        print("#" * 78)

        action = post_calculation_menu(
            "MEMBRANE PERMEATION RATE",
            "Calculate another permeating molecule / membrane system",
        )

        if action == "menu":
            return False

        if action == "exit":
            return True

        reuse = yes_no(
            "Reuse the current GRO/TPR/XTC files?",
            True,
        )

        if reuse:
            continue

        files = resolve_analysis_files(
            [
                (
                    "gro",
                    "gro",
                    "GRO file",
                ),
                (
                    "tpr",
                    "tpr",
                    "TPR file",
                ),
                (
                    "xtc",
                    "xtc",
                    "XTC trajectory",
                ),
            ],
            state,
        )

        gro = files[
            "gro"
        ]
        tpr = files[
            "tpr"
        ]
        xtc = files[
            "xtc"
        ]



def record_recoverable_error(
    parameter_name: str,
    exc: Exception,
    unexpected: bool = False,
) -> Path | None:
    """Append error information to a log without closing the program."""
    try:
        OUTDIR.mkdir(parents=True, exist_ok=True)
        log_file = timestamped_log_path(OUTDIR / "analysis_error.log")

        with log_file.open("w", encoding="utf-8") as fh:
            fh.write(log_timestamp_header())
            fh.write("=" * 78 + "\n")
            fh.write(f"Parameter: {parameter_name}\n")
            fh.write(
                "Error type: "
                + ("UNEXPECTED" if unexpected else "RECOVERABLE")
                + "\n"
            )
            fh.write(f"Message: {exc}\n")

            if unexpected:
                fh.write("\nTraceback:\n")
                fh.write(traceback.format_exc())

            fh.write("=" * 78 + "\n")

        return log_file

    except Exception:
        return None


def error_recovery_menu(
    parameter_name: str,
    exc: Exception,
    unexpected: bool = False,
) -> str:
    """
    Keep the program alive after a failed calculation step.
    """
    log_file = record_recoverable_error(
        parameter_name,
        exc,
        unexpected=unexpected,
    )

    print()
    print("=" * 78)
    print(" CALCULATION ERROR — PROGRAM IS STILL RUNNING")
    print("=" * 78)
    print(f"Parameter : {parameter_name}")
    print(f"Problem   : {exc}")

    if unexpected:
        print(
            "Type      : unexpected Python/runtime error "
            "(details were saved for troubleshooting)"
        )

    if log_file is not None:
        print(f"Error log : {log_file}")

    print()
    print(
        "The program has NOT been closed. "
        "Correct the input/file/setting and try again."
    )
    print()
    print(" [1] Retry the CURRENT parameter")
    print(" [2] Return to the MAIN MENU")
    print(" [3] Exit the program")
    print("=" * 78)

    while True:
        choice = ask("Choose recovery action", "1")

        if choice == "1":
            return "retry"
        if choice == "2":
            return "menu"
        if choice == "3":
            return "exit"

        print("[WARNING] Please choose 1, 2, or 3.")


def run_parameter_with_recovery(
    parameter_name: str,
    runner,
    state: dict,
) -> bool:
    """
    Execute one analysis module inside a recovery loop.

    True means the user explicitly chose to exit the whole program.
    """
    while True:
        try:
            return runner(state)

        except RecoverableAnalysisError as exc:
            action = error_recovery_menu(
                parameter_name,
                exc,
                unexpected=False,
            )

        except KeyboardInterrupt:
            action = error_recovery_menu(
                parameter_name,
                RecoverableAnalysisError(
                    "Keyboard interruption detected. "
                    "The current parameter can be retried."
                ),
                unexpected=False,
            )

        except Exception as exc:
            action = error_recovery_menu(
                parameter_name,
                exc,
                unexpected=True,
            )

        if action == "retry":
            print(
                f"\n[RETRY] Restarting {parameter_name}. "
                "Previously selected file paths remain available as defaults."
            )
            continue

        if action == "menu":
            return False

        return True



def print_main_menu():
    print("\n======================================================================")
    print(" GROMACS MEMBRANE ANALYSIS — MAIN MENU")
    print("======================================================================")
    print(" [1] APL calculation")
    print("     Area per lipid vs Time")
    print()
    print(" [2] Membrane thickness calculation")
    print("     Membrane thickness vs Time")
    print()
    print(" [3] Membrane SASA calculation")
    print("     Membrane SASA vs Time")
    print()
    print(" [4] Lipid order parameter calculation")
    print("     Batch gmx order, Z-axis only")
    print()
    print(" [5] Lipid / molecule MSD and diffusion coefficient")
    print("     Residue/molecule, specific atom/bead, whole-membrane regions,")
    print("     or selected molecule/residue in dynamic membrane regions")
    print()
    print(" [6] Molecular / residue density distribution")
    print("     Whole resname type(s), atom/bead subset(s), or individual instance(s)")
    print()
    print(" [7] Drug MSD / diffusion by dynamic membrane region")
    print("     Dynamic HEADGROUP / BARRIER / TAIL assignment + residence-segment MSD")
    print()
    print(" [8] RMSD / structural deviation calculation")
    print("     Default groups, resnames, atom/bead subsets, or individual molecules")
    print("     Same or separate fitting group; AA and CG compatible")
    print()
    print(" [9] Membrane permeation rate (%) / molecule-count analysis")
    print("     Dynamic membrane Z range + molecules inside/below lower leaflet")
    print("     Initial loading N0 + permeation-rate curve + final rate (%)")
    print()
    print(" [0] Exit")
    print("======================================================================")


def main():
    if not check_environment():
        print("\n[EXIT] Program finished by user request.")
        return

    APL_DIR.mkdir(parents=True, exist_ok=True)
    THICK_DIR.mkdir(parents=True, exist_ok=True)
    SASA_DIR.mkdir(parents=True, exist_ok=True)
    ORDER_DIR.mkdir(parents=True, exist_ok=True)
    MSD_DIR.mkdir(parents=True, exist_ok=True)
    DENSITY_DIR.mkdir(parents=True, exist_ok=True)
    DRUG_REGION_DIR.mkdir(parents=True, exist_ok=True)
    RMSD_DIR.mkdir(parents=True, exist_ok=True)
    PERMEATION_DIR.mkdir(parents=True, exist_ok=True)

    # Remember the last-used paths so returning to another parameter does
    # not force the user to type the same file names repeatedly.
    state = {}

    print("======================================================================")
    print(" GROMACS MEMBRANE STRUCTURAL ANALYSIS v20")
    print(" Parameter-selection / menu-driven version")
    print("======================================================================")
    print(
        "Choose the parameter FIRST. The program then scans the CURRENT folder "
        "for the input files required by that calculation."
    )

    while True:
        print_main_menu()
        choice = ask("Select calculation type", "1")

        if choice == "0":
            print("\n[EXIT] Program finished.")
            return

        if choice == "1":
            exit_program = run_parameter_with_recovery(
                "APL",
                run_apl_menu,
                state,
            )
        elif choice == "2":
            exit_program = run_parameter_with_recovery(
                "MEMBRANE THICKNESS",
                run_thickness_menu,
                state,
            )
        elif choice == "3":
            exit_program = run_parameter_with_recovery(
                "MEMBRANE SASA",
                run_sasa_menu,
                state,
            )
        elif choice == "4":
            exit_program = run_parameter_with_recovery(
                "LIPID ORDER PARAMETER",
                run_order_menu,
                state,
            )
        elif choice == "5":
            exit_program = run_parameter_with_recovery(
                "MSD / DIFFUSION",
                run_msd_menu,
                state,
            )
        elif choice == "6":
            exit_program = run_parameter_with_recovery(
                "DENSITY DISTRIBUTION",
                run_density_menu,
                state,
            )
        elif choice == "7":
            exit_program = run_parameter_with_recovery(
                "DRUG DYNAMIC-REGION DIFFUSION",
                run_dynamic_drug_region_menu,
                state,
            )
        elif choice == "8":
            exit_program = run_parameter_with_recovery(
                "RMSD",
                run_rmsd_menu,
                state,
            )
        elif choice == "9":
            exit_program = run_parameter_with_recovery(
                "MEMBRANE PERMEATION RATE",
                run_permeation_count_menu,
                state,
            )
        else:
            print("[WARNING] Please choose 0, 1, 2, 3, 4, 5, 6, 7, 8, or 9.")
            continue

        if exit_program:
            print("\n[EXIT] Program finished.")
            return


if __name__ == "__main__":
    main()
