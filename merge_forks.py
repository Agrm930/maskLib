#!/usr/bin/env python3
"""
maskLib Fork Merge Summary
==========================
Documents all fork branches, what they changed, and how to apply them.

Run with:  python merge_forks.py [--setup-remotes] [--apply]
"""

import subprocess, sys

FORKS = {
    "probvar":      "https://github.com/probvar/maskLib.git",
    "sebastienLeg": "https://github.com/sebastienLeg/maskLib.git",
    "chakramlab":   "https://github.com/chakramlab/maskLib.git",
}

MERGE_PLAN = [
    {
        "source":   "probvar/chunny-masklib-changes",
        "strategy": "merge",
        "files":    ["src/maskLib/microwaveLib.py"],
        "summary":  "Adds waffle_bumpbond() — bump-bond-aware fill that avoids existing structures.",
        "status":   "MERGED",
    },
    {
        "source":   "probvar/fluxonium (cherry-pick 08c2748)",
        "strategy": "cherry-pick",
        "files":    ["src/maskLib/fluxoniumLib.py"],
        "summary":  "Shorted (x=0) and open (x=0.6) reference structures in test chips 9–12.",
        "status":   "MERGED",
    },
    {
        "source":   "chakramlab/Tom_SQUILL",
        "strategy": "merge (conflict resolved: fluxoniumLib ulayer_edge — took Tom's version)",
        "files":    [
            "src/maskLib/junctionLib.py",
            "src/maskLib/qubitLib.py",
            "src/maskLib/fluxoniumLib.py",
            "src/maskLib/utilities.py",
            "src/maskLib/MaskLib.py",
            "src/maskLib/Entities.py",
            "src/maskLib/dcLib.py",
            "src/maskLib/markerLib.py",
            "src/maskLib/mmWaveLib.py",
            "src/maskLib/microwaveLib.py",
            "src/maskLib/point_logger.py (new)",
            "src/maskLib/check_off_grid_points.py (new)",
            "src/maskLib/check_snap_grid.py (new)",
        ],
        "summary": (
            "Major expansion by Tom DiNapoli (ChakramLab/SQUILL/UT-Austin):\n"
            "   junctionLib: JContact_slot rewritten with unified polygon approach (fixes grid-snapping);\n"
            "                new: JContact_tab, JSingleProbePadLeads, FlagPads, TPads,\n"
            "                     Transmon3DWithShunt, CrossAlignMark\n"
            "   qubitLib:    new Snailmon3D, Fluxonium3D, Transmon3D_leads,\n"
            "                    add_dose_array, add_JJ_dose_array\n"
            "   utilities:   snap_to_grid() — numpy-aware grid snapping for all point types\n"
            "   DRC tools:   point_logger.py, check_off_grid_points.py, check_snap_grid.py\n"
            "   fluxoniumLib: expanded M1_pads right-pad geometry + bandage layer support"
        ),
        "status": "MERGED",
    },
    {
        "source":   "probvar/fluxonium-45-degree-bandaging",
        "strategy": "SKIPPED",
        "files":    ["src/maskLib/fluxoniumLib.py"],
        "summary":  "Experimental 45° angled bandage junction. Author notes 'not currently working'.",
        "status":   "SKIPPED — experimental/incomplete",
    },
    {
        "source":   "chakramlab/tom-initial-mods",
        "strategy": "SKIPPED",
        "files":    [],
        "summary":  "Superseded by chakramlab/Tom_SQUILL (Tom_SQUILL is a superset).",
        "status":   "SKIPPED — superseded by Tom_SQUILL",
    },
    {
        "source":   "sebastienLeg/master",
        "strategy": "N/A",
        "files":    [],
        "summary":  "No commits ahead of origin/master. Nothing to merge.",
        "status":   "UP TO DATE",
    },
    {
        "source":   "chakramlab/master",
        "strategy": "N/A",
        "files":    [],
        "summary":  "No commits ahead of origin/master. Nothing to merge.",
        "status":   "UP TO DATE",
    },
]


def run(cmd, check=True):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if check and result.returncode != 0:
        print(f"ERROR: {result.stderr.strip()}")
        sys.exit(1)
    return result.stdout.strip()


def setup_remotes():
    print("\n=== Setting up fork remotes ===")
    existing = run("git remote").splitlines()
    for name, url in FORKS.items():
        if name not in existing:
            run(f"git remote add {name} {url}")
            print(f"  Added remote: {name} -> {url}")
        else:
            print(f"  Remote already exists: {name}")
    print("\n  Fetching all remotes...")
    run("git fetch --all")
    print("  Done.")


def print_summary():
    icons = {
        "MERGED":                          "✅",
        "SKIPPED — experimental/incomplete": "⚠️ ",
        "SKIPPED — superseded by Tom_SQUILL": "⏭️ ",
        "UP TO DATE":                      "✔️ ",
    }
    print("\n" + "=" * 62)
    print("  FORK MERGE SUMMARY  (master-merged vs origin/master)")
    print("=" * 62)
    for item in MERGE_PLAN:
        icon = icons.get(item["status"], "❓")
        print(f"\n{icon}  {item['source']}")
        print(f"   Status  : {item['status']}")
        print(f"   Strategy: {item['strategy']}")
        if item["files"]:
            print(f"   Files   : {item['files'][0]}")
            for f in item["files"][1:]:
                print(f"              {f}")
        print(f"   Change  : {item['summary']}")

    print("\n" + "=" * 62)
    print("  WHAT'S NEW IN master-merged")
    print("=" * 62)
    print("""
  microwaveLib.py  + waffle_bumpbond() — bump-bond fill avoiding existing geometry
  junctionLib.py   + unified-polygon JContact_slot (no grid snapping artifacts)
                   + JContact_tab, JSingleProbePadLeads, FlagPads, TPads
                   + Transmon3DWithShunt, CrossAlignMark
  qubitLib.py      + Snailmon3D, Fluxonium3D, Transmon3D_leads
                   + add_dose_array, add_JJ_dose_array
  utilities.py     + snap_to_grid(pt, grid) — handles scalar/tuple/array inputs
  fluxoniumLib.py  + StandardTestChip chips 9–12: shorted/open reference structures
                   + Expanded M1_pads right-pad + bandage layer geometry
  point_logger.py  NEW — log and export coordinate points for DRC
  check_off_grid_points.py  NEW — flag points outside chip boundary
  check_snap_grid.py        NEW — verify all points land on SQUILL snap grid
""")


def apply_merges():
    print("\n=== Applying merges to current branch ===")
    print("  Ensure you have run --setup-remotes first.\n")

    print("  Step 1: probvar/chunny-masklib-changes (waffle_bumpbond)...")
    run("git merge probvar/chunny-masklib-changes --no-ff "
        "-m 'Merge: add waffle_bumpbond from chunny-masklib-changes'")

    print("  Step 2: probvar/fluxonium shorted/open config (08c2748)...")
    run("git cherry-pick 08c2748")

    print("  Step 3: chakramlab/Tom_SQUILL (SNAIL, 3D qubits, DRC tools)...")
    ret = subprocess.run("git merge chakramlab/Tom_SQUILL --no-ff "
                         "-m 'Merge chakramlab/Tom_SQUILL'",
                         shell=True, capture_output=True, text=True)
    if ret.returncode != 0:
        print("  ⚠️  Merge conflict in fluxoniumLib.py — auto-resolving (take theirs)...")
        subprocess.run("git checkout --theirs src/maskLib/fluxoniumLib.py", shell=True)
        subprocess.run("git add src/maskLib/fluxoniumLib.py", shell=True)
        subprocess.run("git commit --no-edit", shell=True)
        print("  Conflict resolved.")
    else:
        print("  Clean merge.")

    print("\n  All merges applied. Run `git log --oneline -6` to verify.")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="maskLib fork merge helper")
    p.add_argument("--setup-remotes", action="store_true")
    p.add_argument("--apply", action="store_true")
    args = p.parse_args()

    print_summary()

    if args.setup_remotes:
        setup_remotes()
    if args.apply:
        apply_merges()
    elif not args.setup_remotes:
        print("Tip: run with --setup-remotes to add all fork remotes,")
        print("     or   --apply          to apply all merges to a fresh clone.")
