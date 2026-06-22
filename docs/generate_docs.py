#!/usr/bin/env python3
"""
ClassyMC Documentation Generator
==================================
Parses Fortran source files to extract ProcessIO parameters and generates
Sphinx RST documentation under docs/auto/.

Usage:
    cd docs/
    python generate_docs.py

The script reads:
  - Script_*.f90   → maps input keywords to Fortran type names
  - Each type's source file → extracts description and ProcessIO parameters

Output:
  docs/auto/<category>/index.rst   (category overview + type list)
  docs/auto/<category>/<type>.rst  (per-type parameter reference)
  docs/auto/set_commands.rst       (global 'set' command reference)
"""

import re
import sys
from pathlib import Path
from textwrap import dedent, indent

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).parent
SRC_DIR = SCRIPT_DIR.parent / "src"
AUTO_DIR = SCRIPT_DIR / "auto"

# ---------------------------------------------------------------------------
# Category definitions
# ---------------------------------------------------------------------------
# Each entry describes one input category.
#   id          directory name under auto/
#   title       RST section title
#   script      Script_*.f90 file that dispatches this category
#   intro       introductory paragraph for the category index
#   creation    how the type is selected in the input script
#   modify      how parameters are set after creation (None = in .clFF file)
#   ifdef_note  dict mapping ifdef-guard → explanatory note

CATEGORIES = [
    {
        "id": "forcefields",
        "title": "Forcefield Types",
        "script": "Script_FieldType.f90",
        "intro": (
            "Forcefield types define the energy model used to compute interactions "
            "between atoms.  They are specified in the forcefield file (``.clFF``) "
            "using the ``forcefieldtype`` block.  Multiple forcefields may be "
            "created and assigned to different simulation boxes."
        ),
        "creation": dedent("""\
            .. code-block:: none

               forcefieldtype
                 <type_keyword>
               end_forcefieldtype

               forcefield 1
                 <parameters ...>
               end_forcefield
            """),
        "modify": None,
        "ifdef_note": {
            "AENET":      "Requires the AENet library (build with ``make AENET=1``).",
            "LAMMPS":     "Requires the LAMMPS library interface (build with ``make LAMMPS=1``).",
            "EMBPYTHON":  "Requires the embedded Python build (build with ``make PYTHON=1``).",
        },
    },
    {
        "id": "moves",
        "title": "Monte Carlo Move Types",
        "script": "Script_MCMoves.f90",
        "intro": (
            "Monte Carlo move types define the trial moves attempted during "
            "simulation.  Each move is assigned a relative *probability weight*; "
            "moves are chosen proportional to their weight."
        ),
        "creation": dedent("""\
            .. code-block:: none

               create moves
                 <type_keyword> <probability_weight>
               end_create
            """),
        "modify": "``modify move <N> <parameter> <value>``",
        "ifdef_note": {
            "EMBPYTHON": "Requires the embedded Python build (build with ``make PYTHON=1``).",
        },
    },
    {
        "id": "sampling",
        "title": "Sampling / Acceptance Rule Types",
        "script": "Script_Sampling.f90",
        "intro": (
            "Sampling types control the acceptance criterion applied to trial moves.  "
            "Only one sampling object may exist per simulation."
        ),
        "creation": dedent("""\
            .. code-block:: none

               samplingtype <type_keyword>
            """),
        "modify": "``modify sampling <parameter> <value>``",
        "ifdef_note": {},
    },
    {
        "id": "boxes",
        "title": "Simulation Box Types",
        "script": "Script_SimBoxes.f90",
        "intro": (
            "Simulation box types define the geometry and boundary conditions of "
            "each simulation box.  Multiple boxes can coexist (e.g. for GEMC)."
        ),
        "creation": dedent("""\
            .. code-block:: none

               create boxes
                 <type_keyword>
               end_create

               # or, load coordinates directly:
               create boxes
                 fromfile \"config.clssy\"
               end_create
            """),
        "modify": "``modify box <N> <parameter> <value>``",
        "ifdef_note": {},
    },
    {
        "id": "neighborlists",
        "title": "Neighbor List Types",
        "script": "Script_NeighType.f90",
        "intro": (
            "Neighbor list types determine how the short-range pair interaction "
            "candidate list is built and updated.  Each box can have one or more "
            "neighbor lists."
        ),
        "creation": dedent("""\
            .. code-block:: none

               neighlist <box_number> <type_keyword> <n_lists>
            """),
        "modify": "``modify box <N> neighlist <M> <parameter> <value>``",
        "ifdef_note": {},
    },
    {
        "id": "constraints",
        "title": "Constraint Types",
        "script": "Script_Constraint.f90",
        "intro": (
            "Constraints are applied per-box and reject moves that violate "
            "specified geometric or energetic conditions.  Multiple constraints "
            "can be stacked in a single ``create constraint`` block."
        ),
        "creation": dedent("""\
            .. code-block:: none

               create constraint <box_number>
                 <type_keyword> [arguments ...]
               end_create
            """),
        "modify": None,
        "ifdef_note": {
            "EMBPYTHON": "Requires the embedded Python build (build with ``make PYTHON=1``).",
        },
    },
    {
        "id": "analysis",
        "title": "Analysis Types",
        "script": "Script_AnalysisType.f90",
        "intro": (
            "Analysis objects collect statistics or observables during the "
            "simulation.  They are created in the ``create analysis`` block "
            "and updated at every accepted move."
        ),
        "creation": dedent("""\
            .. code-block:: none

               create analysis
                 <type_keyword> [arguments ...]
               end_create
            """),
        "modify": "``modify analysis <N> <parameter> <value>``",
        "ifdef_note": {
            "EMBPYTHON": "Requires the embedded Python build (build with ``make PYTHON=1``).",
        },
    },
    {
        "id": "trajectories",
        "title": "Trajectory Writer Types",
        "script": "Script_TrajType.f90",
        "intro": (
            "Trajectory writers dump simulation coordinates to disk at a "
            "specified frequency.  Multiple writers with different formats "
            "can be active simultaneously."
        ),
        "creation": dedent("""\
            .. code-block:: none

               create trajectory
                 <type_keyword> <box_number> <frequency> \"filename\"
               end_create
            """),
        "modify": None,
        "ifdef_note": {
            "EMBPYTHON": "Requires the embedded Python build (build with ``make PYTHON=1``).",
        },
    },
    {
        "id": "bonds",
        "title": "Bond Potential Types",
        "script": "Script_BondType.f90",
        "intro": (
            "Bond potential types define the intramolecular bond energy "
            "function.  They are specified in the forcefield file (``.clFF``) "
            "inside a ``bonddef`` block."
        ),
        "creation": dedent("""\
            .. code-block:: none

               bonddef
                 <type_keyword> [parameters ...]
               end_bonddef
            """),
        "modify": None,
        "ifdef_note": {},
    },
    {
        "id": "angles",
        "title": "Angle Potential Types",
        "script": "Script_AngleType.f90",
        "intro": (
            "Angle potential types define the intramolecular angle energy "
            "function.  They are specified in the forcefield file (``.clFF``) "
            "inside an ``angledef`` block."
        ),
        "creation": dedent("""\
            .. code-block:: none

               angledef
                 <type_keyword> [parameters ...]
               end_angledef
            """),
        "modify": None,
        "ifdef_note": {},
    },
    {
        "id": "torsions",
        "title": "Torsion Potential Types",
        "script": "Script_TorsionType.f90",
        "intro": (
            "Torsion potential types define the intramolecular dihedral energy "
            "function.  They are specified in the forcefield file (``.clFF``) "
            "inside a ``torsiondef`` block."
        ),
        "creation": dedent("""\
            .. code-block:: none

               torsiondef
                 <type_keyword> [parameters ...]
               end_torsiondef
            """),
        "modify": None,
        "ifdef_note": {},
    },
]

# ---------------------------------------------------------------------------
# Fortran source parsing utilities
# ---------------------------------------------------------------------------

def _strip_comment(line):
    """Return the line with leading Fortran comment marker stripped."""
    s = line.strip()
    if s.startswith("!"):
        return s.lstrip("!").strip()
    return None


def extract_file_description(content):
    """
    Extract the leading comment block from a Fortran source file.
    Returns a list of non-empty description lines, stopping before any
    'Modifiable Parameters' or parameter-list section.
    """
    lines = content.splitlines()
    description = []
    for line in lines:
        s = line.strip()
        # Stop once we hit a non-comment, non-empty token (module/type/etc.)
        if s and not s.startswith("!") and not s.startswith("#define"):
            break
        if s.startswith("!"):
            text = s.lstrip("!").strip()
            # Stop before modifiable parameters sections
            if re.search(r"modifiable\s+param", text, re.IGNORECASE):
                break
            # Skip pure separator lines like !=====... or !-----
            if text and not re.fullmatch(r"[=\-\*\s]+", text):
                description.append(text)
    return description


def parse_modifiable_params_block(content):
    """
    Look for a '! Modifiable Parameters' comment block in the file header
    and extract the listed parameters.  Returns a list of (keyword, desc) tuples.
    """
    params = []
    lines = content.splitlines()
    in_block = False
    current_kw = None
    current_desc_lines = []

    for line in lines:
        s = line.strip()
        if not s.startswith("!"):
            if in_block:
                break
            continue

        text = s.lstrip("!").strip()

        # Detect the start of the block
        if re.search(r"modifiable\s+param", text, re.IGNORECASE):
            in_block = True
            continue

        if in_block:
            # Empty comment line → flush current param
            if not text:
                if current_kw:
                    params.append((current_kw, " ".join(current_desc_lines).strip()))
                current_kw = None
                current_desc_lines = []
                continue
            # New parameter definition: "  keyword (type) => description"
            m = re.match(r"(\w[\w,/]+)\s*(?:\([^)]*\))?\s*(?:=>|:|-)\s*(.*)", text)
            if m:
                if current_kw:
                    params.append((current_kw, " ".join(current_desc_lines).strip()))
                current_kw = m.group(1).lower()
                current_desc_lines = [m.group(2).strip()] if m.group(2).strip() else []
            else:
                # Continuation line
                if current_kw:
                    current_desc_lines.append(text)

    if current_kw:
        params.append((current_kw, " ".join(current_desc_lines).strip()))

    return params


def parse_processio_keywords(content):
    """
    Extract parameters from a ProcessIO subroutine.

    Finds the subroutine body, then extracts each case statement along with:
      - keywords   (list of strings from the case label)
      - type_hint  ('integer', 'float', 'logical', or '')
      - body_lines (non-comment code lines in the case body)
      - comment    (any Fortran comment immediately before the case line)

    Returns a list of dicts.
    """
    # Locate the ProcessIO subroutine
    sub_re = re.compile(
        r"subroutine\s+\w*[Pp]rocess[Ii][Oo]\w*\s*\([^)]*\)"
        r"(.*?)"
        r"^\s*end\s+subroutine\b",
        re.DOTALL | re.MULTILINE,
    )
    match = sub_re.search(content)
    if not match:
        return []

    body = match.group(1)
    lines = body.splitlines()

    results = []
    i = 0
    pending_comment = ""

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # Collect consecutive comment lines as a pending description
        if stripped.startswith("!"):
            pending_comment = stripped.lstrip("!").strip()
            i += 1
            continue

        # Detect case("...") entry
        case_m = re.match(r'\s*case\s*\(\s*(.+?)\s*\)', line, re.IGNORECASE)
        if case_m:
            case_content = case_m.group(1)
            # Skip 'default' and numeric/selector cases
            if "default" in case_content.lower():
                pending_comment = ""
                i += 1
                continue

            keywords = [
                kw.strip().strip('"').strip("'")
                for kw in re.findall(r'["\']([^"\']+)["\']', case_content)
            ]
            if not keywords:
                pending_comment = ""
                i += 1
                continue

            # Collect the case body until the next case/end select
            body_lines = []
            i += 1
            while i < len(lines):
                bl = lines[i].strip()
                if re.match(r"case\s*\(", bl, re.IGNORECASE) or re.match(r"end\s+select", bl, re.IGNORECASE):
                    break
                body_lines.append(bl)
                i += 1

            body_text = " ".join(body_lines)

            # Infer parameter type from read statement
            type_hint = ""
            if re.search(r"read\s*\([^)]+\)\s*intVal", body_text, re.IGNORECASE):
                type_hint = "integer"
            elif re.search(r"read\s*\([^)]+\)\s*realVal", body_text, re.IGNORECASE):
                type_hint = "float"
            elif re.search(r"read\s*\([^)]+\)\s*logicVal", body_text, re.IGNORECASE):
                type_hint = "logical"
            elif re.search(r"GetXCommand", body_text, re.IGNORECASE):
                # Multiple calls → mixed / complex
                n_reads = len(re.findall(r"GetXCommand", body_text, re.IGNORECASE))
                type_hint = f"{'integer / float' if n_reads > 2 else 'value'}"

            results.append(
                {
                    "keywords": keywords,
                    "type_hint": type_hint,
                    "comment": pending_comment,
                }
            )
            pending_comment = ""
        else:
            pending_comment = ""
            i += 1

    return results


# Files that are known deprecated duplicates and should never be used as
# the canonical source for a type or module.
_DEPRECATED_FILENAMES = {
    "Sampling_JunkCollection.f90",  # superseded by Sampling_Metropolis.f90
    "FF_AENetOld.f90",              # superseded by FF_AENet.f90
    "FF_LJ_Cut_NoNei.f90",          # specialised deprecated variant
    "FF_LJ_Cut.f90",                # superseded by EP-style (FF_EP_LJ_Cut.f90)
    "FF_LJ_Ele_Long.f90",           # superseded
    "FF_Pedone.f90",                # superseded by FF_EP_Pedone_Cut.f90
    "FF_TosiFumi_Cut.f90",          # superseded by FF_EP_TosiFumi_Cut.f90
    "FF_LJ_Shift.f90",              # superseded by EP style
    "FF_LJ_Ele_Cut.f90",            # superseded by EP style
    "Analysis_Forces.f90",          # contains a duplicate 'rdf' type definition
}


def _all_f90_files(preferred_first=True):
    """Yield all .f90 files from src/ and python/, sorted by path,
    with deprecated/junk files filtered out."""
    dirs = [SRC_DIR, SRC_DIR.parent / "python"]
    files = []
    for d in dirs:
        if d.exists():
            files.extend(d.glob("*.f90"))
    return [f for f in sorted(files) if f.name not in _DEPRECATED_FILENAMES]


def parse_use_statements(content):
    """Parse 'use ModuleName, only: TypeName' lines.

    Returns a dict mapping TypeName (case-preserved) → ModuleName.
    """
    use_map = {}
    use_re = re.compile(
        r"^\s*use\s+(\w+)\s*,\s*only\s*:\s*(.+)",
        re.IGNORECASE | re.MULTILINE,
    )
    for m in use_re.finditer(content):
        module_name = m.group(1)
        for token in m.group(2).split(","):
            # Handle 'LocalName => ExportedName' aliases
            t = token.split("=>")[0].strip()
            if t:
                use_map[t] = module_name
    return use_map


def find_module_file(module_name):
    """Search all .f90 files for 'module <module_name>' definition."""
    pattern = re.compile(
        r"^\s*module\s+" + re.escape(module_name) + r"\b",
        re.IGNORECASE | re.MULTILINE,
    )
    for f in _all_f90_files():
        try:
            content = f.read_text(errors="replace")
        except Exception:
            continue
        if pattern.search(content):
            return f, content
    return None, None


def find_type_source_file(type_name, module_name=None):
    """Search .f90 files for the definition of `type_name`.

    If `module_name` is provided, the file containing that module is
    tried first (this avoids picking up duplicate/junk definitions).

    Handles patterns like:
      type, public, extends(Parent) :: TypeName
      type :: TypeName
    """
    type_re = re.compile(
        r"^\s*type\s*(?:,\s*[^:]+)?::\s*" + re.escape(type_name) + r"\b",
        re.IGNORECASE | re.MULTILINE,
    )

    # Preferred: the file that owns the module we know the type comes from
    if module_name:
        preferred_file, preferred_content = find_module_file(module_name)
        if preferred_content and type_re.search(preferred_content):
            return preferred_file, preferred_content
        # Fall back to searching all files if the module file doesn't define it
        if preferred_content:
            return preferred_file, preferred_content

    # General search – return first match
    for f in _all_f90_files():
        try:
            content = f.read_text(errors="replace")
        except Exception:
            continue
        if type_re.search(content):
            return f, content
    return None, None


def parse_script_file(script_path):
    """
    Parse a Script_*.f90 file to build a list of dicts:
      { keyword, type_name, conditional, module_name }

    Only active (non-commented-out) case entries are returned.
    The ``module_name`` field is filled from the ``use ModuleName, only: TypeName``
    statements at the top of the script file, so we know which source file is
    authoritative for each type.
    """
    try:
        content = script_path.read_text(errors="replace")
    except FileNotFoundError:
        return []

    # Build TypeName → ModuleName from 'use' statements in this Script file
    type_to_module = parse_use_statements(content)

    mappings = []
    lines = content.splitlines()
    current_ifdef = None
    i = 0

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # Track conditional compilation blocks
        m = re.match(r"#ifdef\s+(\w+)", stripped)
        if m:
            current_ifdef = m.group(1)
            i += 1
            continue
        if re.match(r"#endif", stripped):
            current_ifdef = None
            i += 1
            continue

        # Skip commented-out lines
        if stripped.startswith("!"):
            i += 1
            continue

        # Look for  case("keyword")  on this line
        case_m = re.match(r'\s*case\s*\(\s*(.+?)\s*\)', line, re.IGNORECASE)
        if case_m:
            case_content = case_m.group(1)
            if "default" in case_content.lower():
                i += 1
                continue

            keywords = [
                kw.strip().strip('"').strip("'")
                for kw in re.findall(r'["\']([^"\']+)["\']', case_content)
            ]
            if not keywords:
                i += 1
                continue

            # Look ahead for the allocate statement (within 5 lines),
            # skipping blank and commented lines.
            type_name = None
            for j in range(i + 1, min(i + 6, len(lines))):
                look = lines[j].strip()
                if not look or look.startswith("!"):
                    continue   # skip comments and blanks
                al_m = re.search(
                    r"allocate\s*\(\s*(\w+)\s*::", look, re.IGNORECASE
                )
                if al_m:
                    type_name = al_m.group(1)
                    break
                # Stop searching once we encounter another case() or end select
                if re.match(r"case\s*\(", look, re.IGNORECASE) or \
                   re.match(r"end\s+select", look, re.IGNORECASE):
                    break

            if type_name:
                # Find which module this type comes from (case-insensitive key lookup)
                module_name = None
                for k, v in type_to_module.items():
                    if k.lower() == type_name.lower():
                        module_name = v
                        break

                for kw in keywords:
                    mappings.append(
                        {
                            "keyword": kw,
                            "type_name": type_name,
                            "conditional": current_ifdef,
                            "module_name": module_name,
                        }
                    )

        i += 1

    return mappings


# ---------------------------------------------------------------------------
# RST generation helpers
# ---------------------------------------------------------------------------

def rst_title(text, char="="):
    return f"{text}\n{char * len(text)}\n"


def rst_subtitle(text, char="-"):
    return f"{text}\n{char * len(text)}\n"


def _safe_id(name):
    """Convert a string to a safe RST/file-system identifier."""
    return re.sub(r"[^a-z0-9_]", "_", name.lower()).strip("_")


def build_param_table(params, header_cols=("Keyword", "Type", "Description")):
    """
    Build a reStructuredText list-table from a list of param dicts.
    Each dict must have: keywords, type_hint, comment
    """
    if not params:
        return "*No configurable parameters.*\n"

    lines = [
        ".. list-table::",
        "   :header-rows: 1",
        "   :widths: 25 12 63",
        "",
    ]

    # Header
    lines.append(f"   * - {header_cols[0]}")
    for col in header_cols[1:]:
        lines.append(f"     - {col}")

    for p in params:
        kw_display = " / ".join(f"``{k}``" for k in p["keywords"])
        type_str = p["type_hint"] if p["type_hint"] else "—"
        desc = p["comment"] if p["comment"] else "*(no description available)*"
        lines.append(f"   * - {kw_display}")
        lines.append(f"     - {type_str}")
        lines.append(f"     - {desc}")

    return "\n".join(lines) + "\n"


def build_type_rst(cat, entry, description_lines, modifiable_params, header_params):
    """
    Generate the RST content for a single type page.

    cat              : category dict
    entry            : {'keyword', 'type_name', 'conditional'}
    description_lines: list of strings from the file header
    modifiable_params: list of param dicts from ProcessIO
    header_params    : list of (kw, desc) from 'Modifiable Parameters' block in header
    """
    kw = entry["keyword"]
    type_name = entry["type_name"]
    cond = entry["conditional"]

    lines = []

    # Page title
    title = f"``{kw}``"
    lines.append(rst_title(title))
    lines.append("")

    # Conditional build warning
    if cond and cond in cat.get("ifdef_note", {}):
        lines.append(f".. warning::")
        lines.append("")
        lines.append(f"   {cat['ifdef_note'][cond]}")
        lines.append("")

    # Description from file header
    if description_lines:
        desc_text = " ".join(description_lines)
        # Wrap at ~80 chars (simple split)
        lines.append(desc_text)
        lines.append("")
    else:
        lines.append(f"*{type_name} — no description available.*")
        lines.append("")

    # Creation syntax
    lines.append(rst_subtitle("Usage"))
    lines.append("")
    lines.append(cat["creation"])
    lines.append("")

    # Modify command
    if cat["modify"]:
        lines.append(rst_subtitle("Modify Command"))
        lines.append("")
        lines.append(f"After creation, parameters are set with:")
        lines.append("")
        lines.append(f".. code-block:: none")
        lines.append("")
        # Indent the modify command
        lines.append(f"   {cat['modify'].strip('`')}")
        lines.append("")

    # Prefer the 'Modifiable Parameters' block from the header if present
    if header_params:
        lines.append(rst_subtitle("Parameters"))
        lines.append("")
        lines.append(".. list-table::")
        lines.append("   :header-rows: 1")
        lines.append("   :widths: 25 75")
        lines.append("")
        lines.append("   * - Keyword")
        lines.append("     - Description")
        for kw_p, desc_p in header_params:
            lines.append(f"   * - ``{kw_p}``")
            lines.append(f"     - {desc_p if desc_p else '*(see source)*'}")
        lines.append("")
    elif modifiable_params:
        lines.append(rst_subtitle("Parameters"))
        lines.append("")
        lines.append(build_param_table(modifiable_params))

    lines.append(f".. rubric:: Source")
    lines.append("")
    lines.append(f"Fortran type: ``{type_name}``")
    lines.append("")

    return "\n".join(lines)


def build_category_index(cat, entries_with_meta):
    """
    Build the category index RST page listing all types.

    entries_with_meta: list of (entry_dict, rst_filename) pairs
    """
    title = cat["title"]
    lines = []
    lines.append(rst_title(title))
    lines.append("")
    lines.append(cat["intro"])
    lines.append("")

    lines.append(rst_subtitle("Creation Syntax"))
    lines.append("")
    lines.append(cat["creation"])
    lines.append("")

    if cat["modify"]:
        lines.append(rst_subtitle("Modifying Parameters"))
        lines.append("")
        lines.append(
            f"After creation, parameters are set using {cat['modify']}."
        )
        lines.append("")

    lines.append(rst_subtitle("Available Types"))
    lines.append("")
    lines.append(".. toctree::")
    lines.append("   :maxdepth: 1")
    lines.append("")
    for entry, rst_file in entries_with_meta:
        kw = entry["keyword"]
        fname = Path(rst_file).stem
        cond_tag = f" *(requires {entry['conditional']})*" if entry["conditional"] else ""
        lines.append(f"   {fname}")
    lines.append("")

    # Quick-reference table
    lines.append(".. list-table:: Quick Reference")
    lines.append("   :header-rows: 1")
    lines.append("   :widths: 30 15 55")
    lines.append("")
    lines.append("   * - Keyword")
    lines.append("     - Optional Build")
    lines.append("     - Type")
    for entry, rst_file in entries_with_meta:
        kw = entry["keyword"]
        cond = entry["conditional"] if entry["conditional"] else "—"
        tn = entry["type_name"]
        lines.append(f"   * - ``{kw}``")
        lines.append(f"     - {cond}")
        lines.append(f"     - ``{tn}``")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 'set' command extraction (Script_Main.f90)
# ---------------------------------------------------------------------------

def generate_set_commands_rst():
    """Parse setCommand in Script_Main.f90 and generate set_commands.rst."""
    script_main = SRC_DIR / "Script_Main.f90"
    if not script_main.exists():
        return None

    content = script_main.read_text(errors="replace")

    # Find the setCommand subroutine
    sub_re = re.compile(
        r"subroutine\s+setCommand\s*\([^)]*\)(.*?)end\s+subroutine",
        re.DOTALL | re.IGNORECASE,
    )
    m = sub_re.search(content)
    if not m:
        return None

    body = m.group(1)

    # Extract case entries + comments
    params = parse_processio_keywords(content)
    # The setCommand select case uses position-2 keyword, same structure

    lines = []
    lines.append(rst_title("Global ``set`` Commands"))
    lines.append("")
    lines.append(
        "Global simulation parameters are set with the ``set`` command at the "
        "top level of the simulation input script."
    )
    lines.append("")
    lines.append(".. code-block:: none")
    lines.append("")
    lines.append("   set <keyword> <value>")
    lines.append("")

    # Manually enumerate the set commands with descriptions (parsed from source)
    SET_PARAMS = [
        ("cycles",          "integer", "Total number of Monte Carlo cycles to run."),
        ("moves",           "integer", "Number of MC trial moves attempted per cycle."),
        ("rng_seed",        "integer", "Random number generator seed.  Use a negative value to seed from the system clock."),
        ("energycheck",     "integer", "Frequency (in cycles) to perform a full energy recomputation for validation.  0 = disabled."),
        ("screenecho",      "logical", "If ``.true.``, print progress to stdout; if ``.false.``, suppress screen output."),
        ("screenfrequency", "integer", "How often (in cycles) to print status to screen."),
        ("configfrequency", "integer", "How often (in cycles) to write configuration files."),
        ("neighskin",       "float",   "Extra distance (Å) added to the neighbor-list cutoff.  The list is rebuilt when an atom moves more than ``neighskin/2``."),
        ("energyunits",     "string",  "Output energy unit string (e.g. ``kcal/mol``, ``kJ/mol``, ``K``)."),
        ("distunits",       "string",  "Output distance unit string (e.g. ``ang``, ``nm``)."),
        ("angleunits",      "string",  "Output angle unit string (e.g. ``deg``, ``rad``)."),
        ("pressureunits",   "string",  "Input/output pressure unit string (e.g. ``atm``, ``bar``, ``Pa``)."),
        ("energytol",       "float",   "Energy convergence tolerance used by the minimizer."),
        ("forcetol",        "float",   "Force convergence tolerance used by the minimizer."),
        ("learnrate",       "float",   "Learning rate (step size) for gradient-descent minimization."),
    ]

    lines.append(".. list-table::")
    lines.append("   :header-rows: 1")
    lines.append("   :widths: 22 12 66")
    lines.append("")
    lines.append("   * - Keyword")
    lines.append("     - Type")
    lines.append("     - Description")
    for kw, typ, desc in SET_PARAMS:
        lines.append(f"   * - ``{kw}``")
        lines.append(f"     - {typ}")
        lines.append(f"     - {desc}")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Box modify-command reference (Box_SimpleBox.f90)
# ---------------------------------------------------------------------------

def generate_box_params_section():
    """
    All box types inherit SimpleBox, which provides the core modify commands.
    Returns an RST snippet with the standard box parameters.
    """
    box_file = SRC_DIR / "Box_SimpleBox.f90"
    if not box_file.exists():
        return ""

    content = box_file.read_text(errors="replace")
    params = parse_processio_keywords(content)

    # Supplement with known descriptions where comments are absent
    KNOWN = {
        "buildfreq":          ("integer", "Frequency (moves) at which the neighbor list is rebuilt."),
        "chempotential":      ("integer float", "Chemical potential for a molecule type: ``chempotential <moltype> <value>``."),
        "rebuildsensitivity": ("float",   "Fractional displacement threshold that triggers a neighbor-list rebuild."),
        "energycalc":         ("integer", "Index of the forcefield object to use for energy calculations."),
        "energyrecompute":    ("logical", "Force a full energy recomputation every step (debug only; very slow)."),
        "neighcut":           ("integer float", "Override the cutoff for a specific neighbor list: ``neighcut <list_id> <rcut>``."),
        "neighlist":          ("integer …", "Delegate further settings to a neighbor-list object: ``neighlist <list_id> <param> <value>``."),
        "pressure":           ("float",   "Target pressure for NPT simulations."),
        "temperature":        ("float",   "Simulation temperature (also sets ``beta = 1/T``)."),
    }

    for p in params:
        kw = p["keywords"][0].lower()
        if kw in KNOWN:
            if not p["type_hint"]:
                p["type_hint"] = KNOWN[kw][0]
            if not p["comment"]:
                p["comment"] = KNOWN[kw][1]

    return build_param_table(params)


# ---------------------------------------------------------------------------
# Main generation loop
# ---------------------------------------------------------------------------

def process_category(cat):
    """Generate all RST files for one category. Returns list of generated files."""
    cat_dir = AUTO_DIR / cat["id"]
    cat_dir.mkdir(parents=True, exist_ok=True)

    script_path = SRC_DIR / cat["script"]
    raw_entries = parse_script_file(script_path)

    if not raw_entries:
        print(f"  [!] No entries found in {cat['script']}")
        return []

    # De-duplicate by type_name (keep first keyword occurrence as canonical)
    seen_types = {}
    for entry in raw_entries:
        tn = entry["type_name"]
        if tn not in seen_types:
            seen_types[tn] = entry
        else:
            # Append alternate keyword to the canonical entry's type info
            pass

    # Group alternate keywords under the same type
    by_type = {}
    for entry in raw_entries:
        tn = entry["type_name"]
        if tn not in by_type:
            by_type[tn] = {"primary": entry, "alternates": []}
        else:
            by_type[tn]["alternates"].append(entry["keyword"])

    generated = []
    entries_with_meta = []

    for type_name, info in by_type.items():
        primary = info["primary"]
        alternates = info["alternates"]

        kw = primary["keyword"]
        rst_stem = _safe_id(kw)
        rst_path = cat_dir / f"{rst_stem}.rst"

        # Find source file – use module_name hint to select the canonical file
        module_name = primary.get("module_name")
        src_file, src_content = find_type_source_file(type_name, module_name)
        if src_content is None:
            print(f"  [!] Source file not found for type {type_name}")
            desc_lines = []
            processio_params = []
            header_params = []
        else:
            desc_lines = extract_file_description(src_content)
            processio_params = parse_processio_keywords(src_content)
            header_params = parse_modifiable_params_block(src_content)

        # Add alternate keywords to the entry metadata
        if alternates:
            primary = dict(primary)
            all_kws = [primary["keyword"]] + alternates
            primary["_all_keywords"] = all_kws

        rst_content = build_type_rst(
            cat, primary, desc_lines, processio_params, header_params
        )

        # Append alternate keyword aliases
        if alternates:
            alias_note = (
                f"\n.. note::\n\n"
                f"   This type can also be selected using the keyword(s): "
                + ", ".join(f"``{a}``" for a in alternates) + ".\n"
            )
            rst_content += alias_note

        rst_path.write_text(rst_content)
        generated.append(rst_path)
        entries_with_meta.append((primary, rst_path))

        src_note = f" ({src_file.name})" if src_file else " (source not found)"
        print(f"  + {kw:30s} → {type_name}{src_note}")

    # Generate category index
    index_content = build_category_index(cat, entries_with_meta)
    index_path = cat_dir / "index.rst"
    index_path.write_text(index_content)
    generated.append(index_path)

    return generated


def generate_set_commands():
    """Write the set_commands.rst file."""
    content = generate_set_commands_rst()
    if content:
        out = AUTO_DIR / "set_commands.rst"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(content)
        print(f"  + set_commands.rst")
        return out
    return None


def generate_box_params():
    """
    The box parameters RST is embedded inside the boxes/simplebox.rst by the
    build_type_rst function.  This function is kept for standalone reference.
    """
    box_section = generate_box_params_section()
    return box_section


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    AUTO_DIR.mkdir(parents=True, exist_ok=True)

    all_files = []

    print("Generating set command reference...")
    f = generate_set_commands()
    if f:
        all_files.append(f)

    for cat in CATEGORIES:
        print(f"\nCategory: {cat['title']}")
        files = process_category(cat)
        all_files.extend(files)

    print(f"\nDone — {len(all_files)} RST files written to {AUTO_DIR}/")


if __name__ == "__main__":
    main()
