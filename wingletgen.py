"""
generate_winglets_openvsp.py

Creates a base wing and then generates multiple variants with different
winglet types using the OpenVSP Python API. Each variant is saved as a .vsp3.

Requires OpenVSP python bindings (importable as `openvsp`).
"""

import os
import math
import itertools
import openvsp as vsp

OUT_DIR = "output_winglets"
os.makedirs(OUT_DIR, exist_ok=True)

def create_base_wing():
    """
    Create a simple base wing geometry and return its geom_id.
    Adjust base wing parameters here as desired.
    """
    wing_id = vsp.AddGeom("WING")
    # Set some reasonable baseline planform parms
    span_pid = vsp.GetParm(wing_id, "Span", "Plan")
    root_chord_pid = vsp.GetParm(wing_id, "Root_Chord", "Plan")
    tip_chord_pid = vsp.GetParm(wing_id, "Tip_Chord", "Plan")
    sweep_pid = vsp.GetParm(wing_id, "Sweep", "Plan")
    dihedral_pid = vsp.GetParm(wing_id, "Dihedral", "XForm")
    # Set values (use SetParmValUpdate to force update)
    vsp.SetParmValUpdate(span_pid, 20.0)        # 20 m span baseline
    vsp.SetParmValUpdate(root_chord_pid, 3.5)
    vsp.SetParmValUpdate(tip_chord_pid, 1.0)
    vsp.SetParmValUpdate(sweep_pid, 20.0)       # deg
    vsp.SetParmValUpdate(dihedral_pid, 3.0)    # deg
    vsp.Update()
    return wing_id

def make_winglet(wing_parent_id, category, params):
    """
    Add a winglet geometry (another WING geom) attached to the tip of the base wing.
    category: string (one of the 8 categories)
    params: dict with keys span, tip_chord, sweep, cant (dihedral), twist, extra
    """
    # Add a new wing geom as a child of the base wing so it is grouped in the model.
    wl_id = vsp.AddGeom("WING", wing_parent_id)

    # Most wing parameters live in Plan or XForm groups:
    span_pid = vsp.GetParm(wl_id, "Span", "Plan")
    root_chord_pid = vsp.GetParm(wl_id, "Root_Chord", "Plan")
    tip_chord_pid = vsp.GetParm(wl_id, "Tip_Chord", "Plan")
    sweep_pid = vsp.GetParm(wl_id, "Sweep", "Plan")
    dihedral_pid = vsp.GetParm(wl_id, "Dihedral", "XForm")
    twist_pid = vsp.GetParm(wl_id, "Twist", "XSec_1")  # twist is often per-xsec; this is a best-effort
    # Note: specific parm names/locations may differ by OpenVSP version; check API if needed.

    # Set basic geometry values
    vsp.SetParmValUpdate(span_pid, params["span"])
    # Make the winglet rooted with small root chord to blend to the wing tip
    vsp.SetParmValUpdate(root_chord_pid, max(0.05, params["root_chord"])) 
    vsp.SetParmValUpdate(tip_chord_pid, params["tip_chord"])
    vsp.SetParmValUpdate(sweep_pid, params["sweep"])
    vsp.SetParmValUpdate(dihedral_pid, params["cant"])
    # Try set twist if available - use SetParmValUpdate by name is safer for twist group indexing
    try:
        vsp.SetParmValUpdate(twist_pid, params["twist"])
    except Exception:
        # Twist parm location varies; ignore if not found
        pass

    # Apply category-specific modifications (approximate; many winglet types are conceptual)
    if category == "traditional":
        # Traditional: simple vertical/dihedral tip with small sweep
        # No special extra changes; keep small toe-in/out via "Sweep"
        pass

    elif category == "blended":
        # Blended winglet: emulate blend by increasing root chord and a small additional bend
        # Increase root chord to smooth transition
        vsp.SetParmValUpdate(root_chord_pid, params["root_chord"] * 1.6)

    elif category == "fences":
        # Wingtip fences: represent as two small vertical plates.
        # We approximate fence by creating a narrow wing (thin span, large chord) rotated slightly.
        # Create a second tiny "fence" wing attached to same parent, offset slightly
        fence = vsp.AddGeom("WING", wing_parent_id)
        f_span = params["span"] * 0.25
        f_root = params["root_chord"] * 0.4
        f_tip = params["tip_chord"] * 0.2
        f_span_pid = vsp.GetParm(fence, "Span", "Plan")
        f_root_pid = vsp.GetParm(fence, "Root_Chord", "Plan")
        f_tip_pid = vsp.GetParm(fence, "Tip_Chord", "Plan")
        f_dihedral_pid = vsp.GetParm(fence, "Dihedral", "XForm")
        vsp.SetParmValUpdate(f_span_pid, f_span)
        vsp.SetParmValUpdate(f_root_pid, f_root)
        vsp.SetParmValUpdate(f_tip_pid, f_tip)
        vsp.SetParmValUpdate(f_dihedral_pid, params["cant"] + 85.0)  # make near-vertical

    elif category == "sharklet":
        # Sharklet (Airbus style): taller, more swept, slender
        vsp.SetParmValUpdate(span_pid, params["span"] * 1.2)
        vsp.SetParmValUpdate(sweep_pid, params["sweep"] + 8.0)

    elif category == "split_scimitar":
        # Split Scimitar: emulate by creating two stacked small winglets: one up, one down (or fore/aft)
        sc_up = vsp.AddGeom("WING", wing_parent_id)
        sc_dn = vsp.AddGeom("WING", wing_parent_id)
        for sc, sign in [(sc_up, 1.0), (sc_dn, -1.0)]:
            s_span_pid = vsp.GetParm(sc, "Span", "Plan")
            s_root_pid = vsp.GetParm(sc, "Root_Chord", "Plan")
            s_tip_pid = vsp.GetParm(sc, "Tip_Chord", "Plan")
            s_dihedral_pid = vsp.GetParm(sc, "Dihedral", "XForm")
            vsp.SetParmValUpdate(s_span_pid, params["span"] * 0.6)
            vsp.SetParmValUpdate(s_root_pid, params["root_chord"] * 0.5)
            vsp.SetParmValUpdate(s_tip_pid, params["tip_chord"] * 0.35)
            vsp.SetParmValUpdate(s_dihedral_pid, params["cant"] + sign * 35.0)

    elif category == "raked":
        # Raked: long span, shallow cant (almost continuation of wing)
        vsp.SetParmValUpdate(span_pid, params["span"] * 1.8)
        vsp.SetParmValUpdate(dihedral_pid, params["cant"] * 0.25)
        vsp.SetParmValUpdate(sweep_pid, params["sweep"] - 10.0)

    elif category == "spiroid":
        # Spiroid is a loop / ring — hard to do exactly in parametric wings.
        # Approximate by creating a highly swept circular tip: increase chord & sweep and negative twist.
        vsp.SetParmValUpdate(tip_chord_pid, params["tip_chord"] * 2.5)
        vsp.SetParmValUpdate(sweep_pid, params["sweep"] + 40.0)
        vsp.SetParmValUpdate(dihedral_pid, params["cant"] + 70.0)

    elif category == "wing_grid":
        # Wing grid (lattice): approximate by creating several narrow winglets spaced across a small arc.
        n_grid = max(2, int(params.get("extra", 3)))
        spacing = params["span"] / (n_grid + 1)
        for i in range(n_grid):
            g = vsp.AddGeom("WING", wing_parent_id)
            g_span_pid = vsp.GetParm(g, "Span", "Plan")
            g_root_pid = vsp.GetParm(g, "Root_Chord", "Plan")
            g_tip_pid = vsp.GetParm(g, "Tip_Chord", "Plan")
            g_dihedral_pid = vsp.GetParm(g, "Dihedral", "XForm")
            vsp.SetParmValUpdate(g_span_pid, spacing * 0.8)
            vsp.SetParmValUpdate(g_root_pid, params["root_chord"] * 0.25)
            vsp.SetParmValUpdate(g_tip_pid, params["tip_chord"] * 0.2)
            # stagger dihedral to make grid-like appearance
            vsp.SetParmValUpdate(g_dihedral_pid, params["cant"] + (i - n_grid/2.0) * 10.0)

    # Final update for this winglet and any children geoms added
    vsp.Update()
    return wl_id

def save_model(name):
    """Write the current VSP model to a .vsp3 file with name in OUT_DIR."""
    fname = os.path.join(OUT_DIR, f"{name}.vsp3")
    # Set VSP file name and write
    vsp.SetVSP3FileName(fname)
    # WriteVSPFile takes (filename, set) — use SET_ALL constant from API
    try:
        vsp.WriteVSPFile(vsp.GetVSPFileName(), vsp.SET_ALL)
    except Exception:
        # older API variants might accept WriteVSPFile(fname, vsp.SET_ALL)
        vsp.WriteVSPFile(fname, vsp.SET_ALL)
    print("Saved:", fname)

def generate_variants():
    base_wing = create_base_wing()

    categories = [
        "traditional", "blended", "fences", "sharklet",
        "split_scimitar", "raked", "spiroid", "wing_grid"
    ]

    # Variation ranges (coarse grid); adjust density if you want more files
    span_variants = [0.6, 0.9, 1.2]    # multipliers relative to mini winglet reference
    tip_chord_variants = [0.15, 0.25]  # meters
    root_chord_variants = [0.6, 0.9]   # meters
    sweep_variants = [10.0, 20.0, 30.0]  # degrees
    cant_variants = [60.0, 75.0, 90.0]    # degrees (90 = vertical)
    twist_variants = [-2.0, 0.0, 2.0]    # degrees
    extra_param = [2, 3, 4]              # used by some categories (grid count, etc)

    # Use a small sample per category to keep total file count manageable
    variant_index = 0
    for cat in categories:
        # For each category pick a few combinations
        for span_mult, tip_c, root_c, sweep, cant, twist, extra in itertools.islice(
                itertools.product(span_variants, tip_chord_variants, root_chord_variants,
                                  sweep_variants, cant_variants, twist_variants, extra_param), 0, 9):
            # Clear model to base wing only (so variants don't accumulate)
            vsp.ClearVSPModel()
            base_wing = create_base_wing()

            # Build a nominal winglet parameter set (meters)
            base_span_ref = 2.0  # baseline winglet span (m)
            params = {
                "span": base_span_ref * span_mult,
                "tip_chord": tip_c,
                "root_chord": root_c,
                "sweep": sweep,
                "cant": cant,
                "twist": twist,
                "extra": extra
            }

            # Create the winglet variant
            wl_id = make_winglet(base_wing, cat, params)

            # Save the model
            fname = f"wing_{cat}_v{variant_index:03d}_s{params['span']:.2f}_tc{params['tip_chord']:.2f}"
            save_model(fname)
            variant_index += 1

    print(f"Generated {variant_index} variant files in '{OUT_DIR}'")

if __name__ == "__main__":
    generate_variants()