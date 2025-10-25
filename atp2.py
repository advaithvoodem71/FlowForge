"""
generate_winglets_stl.py

Generates a single base wing and exports multiple variants each with a different
winglet type. Outputs ASCII STL files in ./output_stls/

Dependencies: numpy
Run: python generate_winglets_stl.py
"""

import os
import math
import numpy as np

OUTPUT_DIR = "output_stls"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Geometry helper functions

def naca4_symmetric(thickness=0.12, n_points=40):
    """
    Returns x, y coordinates (upper and lower stacked) for a symmetric NACA-like
    thickness distribution from 0..1 chord.
    Upper: (x, +y_t), Lower: (x, -y_t). x from 0..1 (leading edge to trailing edge)
    """
    x = np.linspace(0.0, 1.0, n_points)
    t = thickness
    # thickness distribution (standard 4-digit formula)
    y_t = 5 * t * (0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)
    # create coordinates: leading edge to trailing edge on upper, trailing to leading on lower for closed loop
    xu = x
    yu = +y_t
    xl = x[::-1]
    yl = -y_t[::-1]
    xs = np.concatenate([xu, xl[1:]])  # remove duplicate trailing edge
    ys = np.concatenate([yu, yl[1:]])
    pts = np.vstack([xs, ys]).T  # shape (2n-1, 2)
    return pts

def transform_section(section_pts, chord, x_le, y_le, z_le, sweep_deg=0.0, dihedral_deg=0.0, twist_deg=0.0):
    """
    Place a 2D airfoil cross-section in 3D.
      - section_pts: Nx2 (x/c, y/c)
      - chord: physical chord length
      - x_le, y_le, z_le: leading-edge location in 3D
      - sweep_deg: sweep angle in degrees rotates the section around Z so the leading edge shifts in x
      - dihedral_deg: cant angle rotates section around local X axis producing vertical offset
      - twist_deg: rotate around spanwise axis (local Y)
    Returns Nx3 points.
    """
    # scale by chord
    arr = np.array(section_pts) * chord
    # start at 2D coords (x_local, z_local) where z is thickness coordinate (vertical)
    x2 = arr[:,0]
    z2 = arr[:,1]
    y2 = np.zeros_like(x2)

    # apply twist around chord (rotate in x-z plane about leading edge)
    th_twist = math.radians(twist_deg)
    cos_t = math.cos(th_twist); sin_t = math.sin(th_twist)
    x2_t = x2 * cos_t - z2 * sin_t
    z2_t = x2 * sin_t + z2 * cos_t

    # apply sweep: shift leading-edge x offset if desired (we handle sweep by shifting x_le per station)
    # apply dihedral (rotate about x axis)
    th_dih = math.radians(dihedral_deg)
    cos_d = math.cos(th_dih); sin_d = math.sin(th_dih)
    y3 = y2 * cos_d - z2_t * sin_d
    z3 = y2 * sin_d + z2_t * cos_d
    x3 = x2_t

    # assemble and translate to (x_le, y_le, z_le)
    pts3d = np.vstack([x3 + x_le, y3 + y_le, z3 + z_le]).T
    return pts3d

def quad_to_tris(v00, v10, v11, v01):
    """Split quad (v00,v10,v11,v01) into two triangles oriented consistently"""
    return [ (v00, v10, v11), (v00, v11, v01) ]

def write_ascii_stl(filename, triangles, name="winglet_model"):
    """Write triangles (list of (p1,p2,p3) where p are arrays or lists of 3 floats) to ASCII STL"""
    with open(filename, "w") as f:
        f.write(f"solid {name}\n")
        for tri in triangles:
            p1, p2, p3 = [np.array(p) for p in tri]
            # compute normal
            normal = np.cross(p2 - p1, p3 - p1)
            norm_len = np.linalg.norm(normal)
            if norm_len == 0:
                normal = np.array([0.0,0.0,0.0])
            else:
                normal = normal / norm_len
            f.write(f"  facet normal {normal[0]:.6e} {normal[1]:.6e} {normal[2]:.6e}\n")
            f.write("    outer loop\n")
            for p in (p1,p2,p3):
                f.write(f"      vertex {p[0]:.6e} {p[1]:.6e} {p[2]:.6e}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write(f"endsolid {name}\n")

# Wing lofting functions

def loft_sections_to_mesh(sections):
    """
    Given a list of sections, each is an Nx3 array (airfoil contour, closed loop),
    returns triangle list connecting successive sections as a closed surface (no caps).
    sections: [sec0_pts (M0x3), sec1_pts (M1x3), ...]
    Implementation: We require same number of points per section; if not, resample.
    """
    # resample all sections to same number of points = max length
    npt = max(len(s) for s in sections)
    resampled = []
    for s in sections:
        # compute parameterization along 2D contour by cumulative chordwise distance
        arr = np.array(s)
        # compute distances
        d = np.sqrt(((np.roll(arr, -1, axis=0) - arr)**2).sum(axis=1))
        cum = np.concatenate([[0.0], np.cumsum(d[:-1])])
        total = cum[-1] + d[-1]
        t = cum / total
        # new t grid
        tnew = np.linspace(0.0, 1.0, npt)
        # interp x,y,z separately with periodic wrap
        xs = np.interp(tnew, t, arr[:,0])
        ys = np.interp(tnew, t, arr[:,1])
        zs = np.interp(tnew, t, arr[:,2])
        resampled.append(np.vstack([xs,ys,zs]).T)
    # build quads between successive sections
    triangles = []
    for i in range(len(resampled)-1):
        A = resampled[i]
        B = resampled[i+1]
        m = len(A)
        for j in range(m):
            jn = (j+1) % m
            v00 = A[j]
            v10 = B[j]
            v11 = B[jn]
            v01 = A[jn]
            triangles.extend(quad_to_tris(v00, v10, v11, v01))
    return triangles

# Wing & winglet assembly

def build_base_wing(span=20.0, root_chord=3.5, tip_chord=1.0, n_span_sections=6, airfoil_thickness=0.12, n_af_points=48):
    """
    Build base wing as list of sections (each section is Nx3 points).
    We only build one half-wing (right side, y >= 0). The wing root at y=0.
    """
    af2d = naca4_symmetric(thickness=airfoil_thickness, n_points=n_af_points)
    sections = []
    # place sections across half-span (0..span/2)
    half_span = span / 2.0
    for i in range(n_span_sections):
        eta = i / (n_span_sections - 1)  # 0..1
        y = eta * half_span
        chord = root_chord + (tip_chord - root_chord) * eta
        # sweep linear: choose moderate leading-edge sweep
        sweep_le = 20.0  # degrees baseline
        x_le = math.tan(math.radians(sweep_le)) * y
        # no dihedral on main wing baseline
        z_le = 0.0
        # transform section
        sec3d = transform_section(af2d, chord, x_le, y, z_le, sweep_deg=0.0, dihedral_deg=0.0, twist_deg=0.0)
        sections.append(sec3d)
    return sections

def attach_winglet_sections(wing_sections, winglet_params, category):
    """
    Given base wing sections (list of Nx3 arrays), create extra sections at the tip
    to represent the winglet. winglet_params: dict with span, chord_scale, cant, sweep, twist, extra.
    category determines special constructions.
    Returns combined sections (wing + winglet sections).
    """
    # copy existing sections
    sections = list(wing_sections)
    tip_sec = sections[-1]  # the tip section (last of base wing)
    # determine wingtip location (leading edge of last section)
    # find leading edge point as min x value index
    tip_coords = np.array(tip_sec)
    le_idx = np.argmin(tip_coords[:,0])
    le_point = tip_coords[le_idx]  # (x,y,z)
    # winglet base parameters
    wl_span = winglet_params.get("span", 2.0)      # meters (outboard from tip)
    chord_scale = winglet_params.get("chord_scale", 0.4)
    cant = winglet_params.get("cant_deg", 75.0)    # degrees (0 flat, 90 vertical)
    sweep_deg = winglet_params.get("sweep_deg", 20.0)
    twist_deg = winglet_params.get("twist_deg", 0.0)
    n_wl_sections = winglet_params.get("n_sections", 3)
    # for fences / grid / multiple small devices, we will return different section lists
    if category == "wingtip_fence":
        # add 2 small vertical fence plates offset slightly outboard and aft
        plates = []
        for f in [0.0, 0.2]:
            # small rectangle-like airfoil (thin)
            rect = np.array([[0.0,0.0,0.0],
                             [0.2,0.0,0.0],
                             [0.2,0.0,0.8],
                             [0.0,0.0,0.8]])
            # shift outward along y and x
            x_le = le_point[0] + 0.1 + f*0.05
            y_le = le_point[1] + 0.03 + f*0.02
            z_le = le_point[2]
            # extrude as tiny wing-like section sequence
            secA = rect + np.array([x_le, y_le, z_le])
            secB = secA + np.array([0.0, wl_span*0.25, 0.0])
            plates.append(secA)
            plates.append(secB)
        # append plates as extra sections and return
        sections.extend(plates)
        return sections

    if category == "wing_grid":
        # multiple narrow vertical slats along a small arc at tip
        n_grid = winglet_params.get("extra", 4)
        for i in range(n_grid):
            frac = (i+1)/(n_grid+1)
            # create a small narrow plate
            width = 0.05
            height = 0.8 * frac
            rect = np.array([[0.0,0.0,0.0],
                             [width,0.0,0.0],
                             [width,0.0,height],
                             [0.0,0.0,height]])
            # position around tip in y direction
            angle = -10.0 + frac*20.0
            x_le = le_point[0] + math.tan(math.radians(angle)) * (wl_span*frac)
            y_le = le_point[1] + wl_span * frac
            z_le = le_point[2]
            secA = rect + np.array([x_le, y_le, z_le])
            secB = rect + np.array([x_le + 0.0, y_le + 0.02, z_le])
            sections.extend([secA, secB])
        return sections

    if category == "spiroid":
        # approximate a loop by adding a semicircular set of sections that bend up and back
        n_loop = 12
        R = wl_span * 0.6
        for k in range(n_loop):
            theta = math.pi * (k / (n_loop - 1))  # 0..pi
            # ring center at tip_le
            x = le_point[0] + R * math.cos(theta) * 0.4
            y = le_point[1] + wl_span * 0.5 * math.sin(theta)
            z = le_point[2] + R * math.sin(theta) * 0.6
            # use a small circular section
            circle = []
            for a in np.linspace(0, 2*math.pi, 12, endpoint=False):
                circle.append([0.02*math.cos(a), 0.0, 0.02*math.sin(a)])
            circle = np.array(circle) + np.array([x,y,z])
            sections.append(circle)
        return sections

    if category == "split_scimitar":
        # emulate split scimitar by having a main up-winglet + a small lower scimitar trailing piece
        # first the main winglet upward sections
        for i in range(n_wl_sections):
            frac = (i+1)/n_wl_sections
            # chord reduces outboard
            chord = (tip_sec[:,0].max() - tip_sec[:,0].min()) * chord_scale * (1.0 - 0.4*frac)
            x_le = le_point[0] + math.tan(math.radians(sweep_deg)) * (wl_span * frac * 0.5)
            y_le = le_point[1] + wl_span * frac
            z_le = le_point[2] + math.tan(math.radians(cant)) * (wl_span * frac)
            # use scaled airfoil
            af2d = naca4_symmetric(thickness=0.08, n_points=32)
            sec = transform_section(af2d, chord, x_le, y_le, z_le, sweep_deg=sweep_deg, dihedral_deg=0.0, twist_deg=twist_deg)
            sections.append(sec)
        # add a small lower downward scimitar
        for i in range(2):
            frac = 0.3 + 0.2*i
            chord = (tip_sec[:,0].max() - tip_sec[:,0].min()) * 0.25
            x_le = le_point[0] + math.tan(math.radians(sweep_deg+10)) * (wl_span * frac)
            y_le = le_point[1] + wl_span * frac
            z_le = le_point[2] - 0.2 * (i+1)
            af2d = naca4_symmetric(thickness=0.06, n_points=28)
            sec = transform_section(af2d, chord, x_le, y_le, z_le, sweep_deg=sweep_deg+10, dihedral_deg=0.0, twist_deg=-5.0)
            sections.append(sec)
        return sections

    # default: create n_wl_sections along cant direction
    for i in range(n_wl_sections):
        frac = (i+1)/n_wl_sections
        chord = (tip_sec[:,0].max() - tip_sec[:,0].min()) * chord_scale * (1.0 - 0.3*frac)
        # winglet leading edge location moves outboard and upward according to cant
        x_le = le_point[0] + math.tan(math.radians(sweep_deg)) * (wl_span * frac * 0.6)
        y_le = le_point[1] + wl_span * frac
        z_le = le_point[2] + math.tan(math.radians(cant)) * (wl_span * frac)
        af2d = naca4_symmetric(thickness=0.08, n_points=32)
        sec = transform_section(af2d, chord, x_le, y_le, z_le, sweep_deg=sweep_deg, dihedral_deg=0.0, twist_deg=twist_deg)
        sections.append(sec)

    # for blended type, add an extra close-in transition section that smooths the junction
    if category == "blended":
        # add a small transition ring that lies between tip and first winglet sec
        # create radial interpolation points and insert at index -2
        pass

    if category == "raked":
        # raked wingtips are basically long continuation; ensure sections extend with low cant
        pass

    if category == "sharklet":
        # make slightly taller, more swept
        pass

    return sections

# Variant generation

def generate_variants_and_write():
    # base parameters
    base_span = 40.0   # total span for full wing (we build half wing of 20 m)
    root_chord = 4.0
    tip_chord = 1.4
    base_sections = build_base_wing(span=base_span, root_chord=root_chord, tip_chord=tip_chord,
                                    n_span_sections=8, airfoil_thickness=0.12, n_af_points=64)

    categories = [
        "traditional", "blended", "wingtip_fence", "sharklet",
        "split_scimitar", "raked", "spiroid", "wing_grid"
    ]

    # small parameter grid per category (keeps file count moderate)
    variant_count = 0
    for cat in categories:
        # define param ranges tuned per category
        if cat == "traditional":
            spans = [1.6, 2.0, 2.6]          # m
            chord_scales = [0.35, 0.45]
            cants = [70, 75, 82]
            sweeps = [10, 20]
        elif cat == "blended":
            spans = [2.0, 2.6]
            chord_scales = [0.4, 0.55]
            cants = [68, 75]
            sweeps = [5, 18]
        elif cat == "wingtip_fence":
            spans = [0.8, 1.2]
            chord_scales = [0.15]
            cants = [90]
            sweeps = [0]
        elif cat == "sharklet":
            spans = [2.0, 2.4]
            chord_scales = [0.3, 0.4]
            cants = [78, 85]
            sweeps = [18, 28]
        elif cat == "split_scimitar":
            spans = [2.2, 2.8]
            chord_scales = [0.25, 0.35]
            cants = [78, 83]
            sweeps = [15, 25]
        elif cat == "raked":
            spans = [3.2, 3.8]
            chord_scales = [0.5, 0.7]
            cants = [15, 28]
            sweeps = [30, 45]
        elif cat == "spiroid":
            spans = [2.6, 3.2]
            chord_scales = [0.12]
            cants = [100, 180]
            sweeps = [0]
        elif cat == "wing_grid":
            spans = [1.6, 2.2]
            chord_scales = [0.2]
            cants = [90]
            sweeps = [0]
        else:
            spans = [2.0]; chord_scales=[0.4]; cants=[75]; sweeps=[20]

        # iterate small grid and make files
        for span in spans:
            for cs in chord_scales:
                for cant in cants:
                    for sw in sweeps:
                        # build sections copy
                        wing_secs = list(base_sections)
                        params = {
                            "span": span,
                            "chord_scale": cs,
                            "cant_deg": cant,
                            "sweep_deg": sw,
                            "twist_deg": 0.0,
                            "n_sections": 4,
                            "extra": 4
                        }
                        assembled_sections = attach_winglet_sections(wing_secs, params, cat)
                        # produce triangles (this connects sequential sections into a surface)
                        triangles = loft_sections_to_mesh(assembled_sections)
                        # Also cap root by triangulating root contour to close volume (optional)
                        # We'll just export the shell (wing half). For full model combine mirrored half if needed.
                        fname = f"{cat}_s{span:.2f}_cs{cs:.2f}_cant{cant}_sw{sw}.stl"
                        fname = fname.replace(" ", "_")
                        path = os.path.join(OUTPUT_DIR, fname)
                        write_ascii_stl(path, triangles, name=f"{cat}")
                        print("Wrote:", path)
                        variant_count += 1
    print(f"Done. Generated {variant_count} STL files in '{OUTPUT_DIR}'")

if __name__ == "__main__":
    generate_variants_and_write()