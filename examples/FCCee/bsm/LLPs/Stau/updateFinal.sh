#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Configuration
# ============================================================

# Folder containing the MadGraph production folders.
#
# Example:
#   MADGRAPH_BASE="/data/dust/user/creusett/SummerSchool/Signal/mc-pythia-main"
#
# Inside this folder:
#
#   MADGRAPH_BASE/
#       run_01/
#       run_02/
#       run_03/
#       ...
#
MADGRAPH_BASE="/afs/desy.de/user/m/mrandria/summer_project/madgraph/MG5_aMC_v3_5_16"

# Name of the banner file inside each MadGraph production folder.
BANNER_NAME="run_01_tag_1_banner.txt"

# Python file containing procDictAdd.
#
# Example:
#   PROC_DICT_FILE="/data/dust/user/creusett/SummerSchool/Signal/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/analysis_final.py"
#
PROC_DICT_FILE="/afs/desy.de/user/m/mrandria/summer_project/Signal/FCCAnalyses/examples/FCCee/bsm/LLPs/Stau/analysis_final1.py"

# ============================================================
# Sanity checks
# ============================================================

if [[ ! -d "$MADGRAPH_BASE" ]]; then
    echo "ERROR: MadGraph base directory does not exist:"
    echo "  $MADGRAPH_BASE"
    exit 1
fi

if [[ ! -f "$PROC_DICT_FILE" ]]; then
    echo "ERROR: Python file does not exist:"
    echo "  $PROC_DICT_FILE"
    exit 1
fi

# ============================================================
# Python updater
# ============================================================

python3 - "$MADGRAPH_BASE" "$PROC_DICT_FILE" "$BANNER_NAME" <<'PY'

import sys
import re
from pathlib import Path


MADGRAPH_BASE = Path(sys.argv[1])
PROC_DICT_FILE = Path(sys.argv[2])
BANNER_NAME = sys.argv[3]


# ============================================================
# Read values from a MadGraph banner
# ============================================================

def read_madgraph_banner(banner_path):

    text = banner_path.read_text()

    # --------------------------------------------------------
    # Number of events
    #
    # Example:
    #
    # #  Number of Events        :       100000
    # --------------------------------------------------------

    match = re.search(
        r"#\s*Number of Events\s*:\s*(\d+)",
        text
    )

    if not match:
        raise RuntimeError(
            "Could not find 'Number of Events'"
        )

    number_of_events = int(match.group(1))


    # --------------------------------------------------------
    # Cross section
    #
    # Example:
    #
    # #  Integrated weight (pb)  :       0.08362259999999999
    #
    # This value is ALREADY in pb.
    #
    # --------------------------------------------------------

    match = re.search(
        r"#\s*Integrated weight\s*\(pb\)\s*:\s*"
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)"
        r"(?:[eE][-+]?\d+)?)",
        text
    )

    if not match:
        raise RuntimeError(
            "Could not find 'Integrated weight (pb)'"
        )

    cross_section_pb = float(match.group(1))


    # --------------------------------------------------------
    # Sum of weights
    #
    # For the generated unweighted sample this is normally the
    # same as Number of Events.
    #
    # We use the number of events from MGGenerationInfo.
    # --------------------------------------------------------

    sum_of_weights = number_of_events


    return (
        number_of_events,
        sum_of_weights,
        cross_section_pb,
    )


# ============================================================
# Format floating-point values
# ============================================================

def format_cross_section(value):
    """
    Use scientific notation, suitable for the Python dictionary.
    """

    return f"{value:.8e}"


# ============================================================
# Read the Python file
# ============================================================

text = PROC_DICT_FILE.read_text()

proc_match = re.search(
    r"procDictAdd\s*=\s*\{(.*?)\n\}",
    text,
    re.DOTALL
)

if not proc_match:
    raise RuntimeError("Could not find procDictAdd block")

proc_text = proc_match.group(1)


# ============================================================
# Find procDictAdd entries
# ============================================================
#
# Expected structure:
#
# 'FCCee_150_stau_1m_ctau_ecm_365': {
#     "numberOfEvents": 100000,
#     "sumOfWeights": 100000,
#     "crossSection": 3.61349000e-11 * 1e9,
#     "kfactor": 1.0,
#     "matchingEfficiency": 1.0
# },
#
# Only uncommented entries are processed.
#
# The three fields modified are:
#
#     numberOfEvents
#     sumOfWeights
#     crossSection
#
# Everything else remains untouched.
# ============================================================


# Match a complete dictionary entry.
#
# The expression deliberately captures the contents of the
# entry and then modifies only the three desired fields.

entry_pattern = re.compile(
    r"""(?P<comment>^\s*\#\s*)?
    (?P<prefix>
        ['"]
        (?P<name>FCCee_[^'"]+)
        ['"]
        \s*:\s*\{
    )
    (?P<body>.*?)
    (?P<suffix>\}
    ,?)
    """,
    re.VERBOSE | re.MULTILINE
)


updated = 0
skipped = 0
failed = 0

def update_field(body, field_name, new_value):

    pattern = re.compile(
        rf'("{re.escape(field_name)}"\s*:\s*)'
        rf'[^,}}]+'
    )

    match = pattern.search(body)

    if not match:
        raise RuntimeError(
            f"Could not find field '{field_name}'"
        )

    return (
        body[:match.start()]
        + match.group(1)
        + new_value
        + body[match.end():]
    )


# ============================================================
# Process every dictionary entry
# ============================================================

def replace_entry(match):

    global updated, skipped, failed

    full_entry = match.group(0)
    name = match.group("name")
    body = match.group("body")

    if match.group("comment"):
        skipped += 1
        return full_entry

    # --------------------------------------------------------
    # Do not touch entries whose first non-whitespace character
    # is '#'.
    #
    # This mainly protects against commented entries.
    # --------------------------------------------------------

    lines = full_entry.splitlines()

    # --------------------------------------------------------
    # Look for corresponding MadGraph directory
    # --------------------------------------------------------

    mg_folder = MADGRAPH_BASE / name
    banner = mg_folder / "Events" / "run_01" / BANNER_NAME

    if not banner.is_file():

        print(
            f"WARNING: MadGraph banner not found for '{name}'"
        )
        print(f"         Expected: {banner}")
        print("         Entry left unchanged.")
        print()

        failed += 1

        return full_entry


    # --------------------------------------------------------
    # Extract MadGraph values
    # --------------------------------------------------------

    try:

        number_of_events, sum_of_weights, cross_section = (
            read_madgraph_banner(banner)
        )

    except Exception as exc:

        print(
            f"ERROR: Could not read '{name}': {exc}"
        )
        print("       Entry left unchanged.")
        print()

        failed += 1

        return full_entry


    # --------------------------------------------------------
    # Update ONLY the three requested fields
    # --------------------------------------------------------

    try:

        body = update_field(
            body,
            "numberOfEvents",
            str(number_of_events)
        )

        body = update_field(
            body,
            "sumOfWeights",
            str(sum_of_weights)
        )

        body = update_field(
            body,
            "crossSection",
            format_cross_section(cross_section)
        )

    except Exception as exc:

        print(
            f"ERROR updating '{name}': {exc}"
        )
        print("       Entry left unchanged.")
        print()

        failed += 1

        return full_entry


    # --------------------------------------------------------
    # Reconstruct entry
    # --------------------------------------------------------

    new_entry = (
        match.group("prefix")
        + body
        + match.group("suffix")
    )


    print(
        f"UPDATED {name}"
    )

    print(
        f"    numberOfEvents = {number_of_events}"
    )

    print(
        f"    sumOfWeights   = {sum_of_weights}"
    )

    print(
        f"    crossSection   = {cross_section:.8e} pb"
    )

    print()

    updated += 1

    return new_entry


# ============================================================
# Apply replacements
# ============================================================

new_proc_text = entry_pattern.sub(
    replace_entry,
    proc_text
)

new_text = (
    text[:proc_match.start(1)]
    + new_proc_text
    + text[proc_match.end(1):]
)


# ============================================================
# Write file
# ============================================================

PROC_DICT_FILE.write_text(new_text)


# ============================================================
# Summary
# ============================================================

print("============================================")
print("procDictAdd update complete")
print("============================================")
print(f"Updated : {updated}")
print(f"Skipped : {skipped}")
print(f"Failed : {failed}")
print(f"File    : {PROC_DICT_FILE}")
print()

PY

echo "Done."