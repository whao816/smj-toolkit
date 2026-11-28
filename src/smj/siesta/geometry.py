"""Periodic Au(111) geometries in Å. Atom numbers in serialised records start at one."""
from __future__ import annotations

from itertools import product
import math

import numpy as np

from ..junction.contacts import rotation_matrix_to
from ..junction.geometry import Atom, strip_thiol_h
from ..phonon.device import fcc111_in_plane_vectors


def _settings(mapping, section, keys):
    if not isinstance(mapping, dict):
        raise ValueError(f"siesta.{section} must be a mapping")
    unknown = mapping.keys() - set(keys)
    if unknown:
        raise ValueError(f"unknown siesta.{section} field(s): {', '.join(map(str, unknown))}")
    return mapping


def _integer(value, name, minimum):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _positive(value, name):
    if isinstance(value, bool) or not math.isfinite(float(value)) or float(value) <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return float(value)


def _planes(lead, first, count, z0):
    """Order each matching block identically to the periodic electrode, even on the right."""
    u, v = np.array(lead["primitive"])
    nx, ny = lead["in_plane"]
    plane = np.array([i*u + j*v for i, j in product(range(nx), range(ny))])
    planes = []
    for k in range(count):
        xyz = plane + ((first + k) % 3) * (u + v)/3
        xyz[:, 2] = z0 + k * lead["spacing_A"]
        planes.append(xyz)
    return np.vstack(planes)


def build_electrode(siesta):
    """Build a periodic ABC principal layer with its repeat vector parallel to z."""
    cfg = _settings(siesta.get("electrode", {}), "electrode",
                    ("a_lattice", "in_plane", "principal_layers", "slab_layers", "movable_layers"))
    a = _positive(cfg.get("a_lattice", 4.07823), "a_lattice")
    nx, ny = cfg.get("in_plane", [3, 3])
    nx, ny = [_integer(n, "in_plane", 1) for n in (nx, ny)]
    n = _integer(cfg.get("principal_layers", 3), "principal_layers", 1)
    if n % 3:
        raise ValueError("principal_layers must be a multiple of 3 for a normal ABC repeat")
    u, v = np.array(fcc111_in_plane_vectors(a, (1, 1)))
    dz = a/math.sqrt(3)
    lead = {"primitive": [u.tolist(), v.tolist()], "in_plane": [nx, ny], "spacing_A": dz,
            "principal_layers": n, "cell": [(nx*u).tolist(), (ny*v).tolist(), [0, 0, n*dz]],
            "elements": ["Au"] * (n*nx*ny)}
    lead["positions"] = _planes(lead, 0, n, 0).tolist()
    return lead


def _flat_sites(u, v, sign):
    hollow = sign*(u+v)/3
    return [("top", 0*u), ("bridge", u/2), ("bridge", v/2), ("bridge", (u+v)/2),
            ("fcc", hollow), ("hcp", -hollow)]


def _contact(cfg, lead, phase, sign):
    _settings(cfg, "contact.left" if sign == 1 else "contact.right",
              ("motif", "site", "height", "movable_tip_atoms"))
    u, v = np.array(lead["primitive"])
    nx, ny = lead["in_plane"]
    # Choose a central lattice point. Positive/negative signs continue ABC towards the gap.
    surface = (nx//2)*u + (ny//2)*v + (phase % 3)*(u+v)/3
    motif = cfg.get("motif", "flat")
    if motif not in ("flat", "adatom", "trimer", "pyramid"):
        raise ValueError(f"unsupported contact motif {motif!r}")
    height = _positive(cfg.get("height", 2.35), "anchor height")
    tips = []
    if motif == "flat":
        site = cfg.get("site", "top")
        offsets = [p for name, p in _flat_sites(u, v, sign) if name == site]
        if not offsets:
            raise ValueError("flat contact requires site top, bridge, fcc or hcp")
        target = surface + offsets[0]
    else:
        site = "center" if motif == "trimer" else "apex"
        if cfg.get("site", site) != site:
            raise ValueError(f"{motif} contact requires local site {site}")
        if motif in ("trimer", "pyramid") and min(nx, ny) < 2:
            raise ValueError(f"{motif} footprint overlaps its periodic image in this in_plane cell")
        step = sign*(u+v)/3
        first = surface + step
        first[2] = sign*lead["spacing_A"]
        tips.append(first)
        if motif != "adatom":
            tips.extend((first + sign*u, first + sign*v))
        target = np.mean(tips, axis=0)
        if motif == "pyramid":
            target[2] += sign*lead["spacing_A"]
            tips.append(target.copy())
    target = target.copy()
    target[2] += sign*height
    return {"motif": motif, "site": site, "height_A": height}, target, np.array(tips).reshape(-1, 3)


def _nearest_registry_shift(left, right, primitive):
    # Exact for the equal-length 60-degree primitive basis.
    base = np.rint((left-right)[:2] @ np.linalg.inv(primitive[:, :2])).astype(int)
    candidates = [(base + np.array(ij)) @ primitive for ij in product((-1, 0, 1), repeat=2)]
    return min(candidates, key=lambda s: np.linalg.norm((right+s-left)[:2]))


def build_device(molecule: list[Atom], anchors, electrode, siesta, *, strip_h=False):
    """Rigidly fit both anchors; close the far ends on the undisturbed ABC lattice."""
    expected = build_electrode(siesta)
    if electrode != expected:
        raise ValueError("electrode geometry does not match the resolved periodic Au input")
    n = electrode["principal_layers"]
    lead_cfg = siesta.get("electrode", {})
    layers = _integer(lead_cfg.get("slab_layers", n+2), "slab_layers", n)
    movable = _integer(lead_cfg.get("movable_layers", 1), "movable_layers", 0)
    if layers-movable < n:
        raise ValueError("movable_layers would relax the electrode matching block")
    left_id, right_id = anchors
    if left_id == right_id or not all(1 <= i <= len(molecule) for i in anchors):
        raise ValueError("invalid molecular anchor ids")
    orig = molecule
    removed = []
    if strip_h:
        molecule, removal = strip_thiol_h(molecule, anchors)
        removed = removal["removed_h_indices"]
        if len(removed) != len(set(removed)):
            raise ValueError("ambiguous hydrogen shared by both molecular anchors")
    retained = [i for i in range(1, len(orig)+1) if i not in removed]
    xyz = np.array([a.coord() for a in molecule])
    left, right = xyz[[retained.index(i) for i in anchors]]
    length = np.linalg.norm(right-left)
    if not np.isfinite(xyz).all() or length < 1.e-8:
        raise ValueError("invalid molecular coordinates or coincident anchors")

    cfg = _settings(siesta.get("contact", {}), "contact", ("left", "right", "gap"))
    lmeta, lt, ltips = _contact(cfg.get("left", {}), electrode, layers-1, 1)
    rmeta, rt, rtips = _contact(cfg.get("right", {}), electrode, -layers, -1)
    shift = _nearest_registry_shift(lt, rt, np.array(electrode["primitive"]))
    rt += shift
    rtips += shift
    lateral = np.linalg.norm((rt-lt)[:2])
    if lateral >= length:
        raise ValueError("contact registry separation exceeds the molecular anchor span")
    gap = math.sqrt(length**2 - lateral**2) + lt[2] - rt[2]
    if cfg.get("gap") is not None:
        desired = _positive(cfg["gap"], "gap")
        free_heights = [(meta, target, sign) for side, meta, target, sign in
                        (("left", lmeta, lt, 1), ("right", rmeta, rt, -1))
                        if "height" not in cfg.get(side, {})]
        if not free_heights and abs(desired-gap) > 1.e-5:
            raise ValueError(f"specified gap conflicts with anchor heights and span (requires {gap:g} Å)")
        for meta, target, sign in free_heights:
            change = (desired-gap)/len(free_heights)
            meta["height_A"] = _positive(meta["height_A"]+change, "height implied by gap")
            target[2] += sign*change
        if free_heights:
            gap = desired
    dz = electrode["spacing_A"]
    lz = (layers-1)*dz
    rz = lz + gap
    lt[2] += lz
    rt[2] += rz
    ltips[:, 2] += lz
    rtips[:, 2] += rz
    rot = np.array(rotation_matrix_to(tuple(right-left), tuple(rt-lt)))
    molecule_xyz = (xyz-left) @ rot.T + lt

    block_l = _planes(electrode, 0, layers, 0)
    block_r = _planes(electrode, -layers, layers, rz)
    positions = np.vstack((block_l, ltips, molecule_xyz, rtips, block_r))
    nm = len(molecule)
    start_mol = len(block_l) + len(ltips) + 1
    start_r = start_mol + nm + len(rtips)
    nat = len(positions)
    nmlead = len(electrode["positions"])
    nx, ny = electrode["in_plane"]
    per_plane = nx*ny
    nfixed = per_plane*(layers-movable)
    fixed = list(range(1, nfixed+1)) + list(range(nat-nfixed+1, nat+1))
    mapping = {str(atom): start_mol+i for i, atom in enumerate(retained)}
    for side, meta, anchor, tiprange, start in (
            ("left", lmeta, left_id, range(len(block_l)+1, start_mol), 1),
            ("right", rmeta, right_id, range(start_mol+nm, start_r), start_r)):
        planes = [list(range(i, i+per_plane)) for i in range(start, start+layers*per_plane, per_plane)]
        selected = cfg.get(side, {}).get("movable_tip_atoms", list(range(1, len(tiprange)+1)))
        if (not isinstance(selected, list)
                or any(type(i) is not int or not 1 <= i <= len(tiprange) for i in selected)
                or len(set(selected)) != len(selected)):
            raise ValueError(f"{side} movable_tip_atoms must list distinct local tip ids in 1..{len(tiprange)}")
        fixed.extend(atom for i, atom in enumerate(tiprange, 1) if i not in selected)
        anchor_h = [i for i, atom in enumerate(orig, 1) if atom.element == "H"
                    and math.dist(orig[anchor-1].coord(), atom.coord()) < 1.70]
        meta.update(anchor_id=mapping[str(anchor)], original_anchor_id=anchor,
                    anchor_element=orig[anchor-1].element, tip_ids=list(tiprange),
                    layer_ids=planes, surface_ids=planes[-1 if side == "left" else 0],
                    original_h_neighbor_ids=anchor_h, h_neighbor_cutoff_A=1.70,
                    removed_original_h_ids=[i for i in anchor_h if i in removed])
    cell = [*electrode["cell"][:2], [0, 0, rz+layers*dz]]
    device = {"cell": cell, "positions": positions.tolist(),
              "orientation": "rigid shortest rotation between original and target anchor axes",
              "constraints": "fixed_ids fixed in all Cartesian directions; movable_ids free in all directions",
              "elements": ["Au"]*(start_mol-1) + [a.element for a in molecule] + ["Au"]*(nat-start_mol-nm+1),
              "gap_A": gap, "fixed_ids": sorted(fixed),
              "movable_ids": sorted(set(range(1, nat+1)).difference(fixed)),
              "mapping": {"original_to_device": mapping, "removed_h_ids": removed, "strip_h": strip_h},
              "contacts": {"left": lmeta, "right": rmeta},
              "electrode_ids": {"left": list(range(1, nmlead+1)), "right": list(range(nat-nmlead+1, nat+1))},
              "semi_infinite_direction": {"left": "-a3", "right": "+a3"}}
    return device


def _neighbors(point, ids, xyz, cell, radius):
    """Every in-plane image in a ball. Small cells can contribute several images of one id."""
    ids = np.asarray(ids, dtype=int)
    delta = xyz[ids-1] - point
    inv = np.linalg.inv(cell[:2, :2])
    frac = delta[:, :2] @ inv
    nearest = -np.rint(frac).astype(int)
    # Reciprocal-vector bound includes the full search ball even for oblique unit cells.
    limits = np.floor(radius*np.linalg.norm(inv, axis=0) + .5).astype(int)
    found = []
    for i, j in product(range(-limits[0], limits[0]+1), range(-limits[1], limits[1]+1)):
        image = nearest + (i, j)
        distances = np.linalg.norm(delta + image @ cell[:2], axis=1)
        for k in np.flatnonzero(distances <= radius):
            found.append({"id": int(ids[k]), "image": image[k].tolist(),
                          "distance_A": float(distances[k])})
    return sorted(found, key=lambda pair: (pair["distance_A"], pair["id"]))


def _contact_state(meta, sign, xyz, lead, cell):
    anchor = xyz[meta["anchor_id"]-1]
    layers = meta["layer_ids"]
    ids = [k for layer in layers for k in layer] + meta["tip_ids"]
    # Mechanical contact diagnostic, in Angstrom, applied to both reference and final coordinates.
    partners = _neighbors(anchor, ids, xyz, cell, 3.35)
    dz = np.array([xyz[np.asarray(group)-1, 2].mean() for group in layers])
    site = "detached"
    primitive = np.asarray(lead["primitive"])
    if partners:
        u, v = primitive
        # Describe registration to the mean surface lattice; individual Au coordination is
        # measured from actual Cartesian positions separately.
        nx, ny = lead["in_plane"]
        reference = xyz[np.asarray(meta["surface_ids"])-1].mean(axis=0) - ((nx-1)*u + (ny-1)*v)/2
        candidates = [(name, reference+off) for name, off in _flat_sites(u, v, sign)]
        def lateral_distance(item):
            target = item[1]
            shift = _nearest_registry_shift(anchor, target, primitive)
            return np.linalg.norm((target+shift-anchor)[:2])
        nearest = min(candidates, key=lateral_distance)
        site = nearest[0] if lateral_distance(nearest) <= .12*np.linalg.norm(u) else "off-site"
        tips = meta["tip_ids"]
        if tips and partners[0]["id"] in tips:
            if meta["motif"] == "trimer":
                # The bonded triad can straddle a lateral boundary in the output; xyz is
                # unwrapped to the assembled reference before this arithmetic.
                centre = xyz[np.array(tips)-1].mean(axis=0)
                shift = _nearest_registry_shift(anchor, centre, cell[:2])
                lateral = np.linalg.norm((centre + shift-anchor)[:2])
                site = "center" if lateral <= .12*np.linalg.norm(u) else "off-center"
            else:
                closest = partners[0]
                local = xyz[closest["id"]-1] + np.asarray(closest["image"]) @ cell[:2]
                if closest["id"] != tips[-1]:
                    site = "base"
                else:
                    site = "apex" if np.linalg.norm((local-anchor)[:2]) <= .12*np.linalg.norm(u) else "off-apex"
    return {"final_site": site, "coordination": len(partners), "neighbors": partners,
            "layer_spacings_A": np.diff(dz).tolist()}


def align_device_positions(device, atoms, cell):
    """Choose XV images nearest the assembled structure after one common translation.

    TranSIESTA compares the ordered block's *Cartesian differences* directly.
    Save these aligned coordinates as its input. The untouched XV retains the raw positions.
    """
    cell = np.asarray(cell)
    xyz = np.array([a.coord() for a in atoms])
    ref = np.asarray(device["positions"])
    delta = xyz-ref-(xyz[0]-ref[0])
    return (ref + delta - np.rint(delta @ np.linalg.inv(cell)) @ cell).tolist()


def inspect_device(device, electrode, final_atoms: list[Atom], final_cell):
    """Audit the two independently ordered electrode blocks and the relaxed interfaces."""
    issues = []
    report = {"passed": False, "issues": issues, "contacts": {}, "electrode_match": {},
              "coordination_cutoff_A": 3.35, "overlap_cutoff_A": 1.3,
              "matching_tolerance_A": 1.e-4}
    if [a.element for a in final_atoms] != device["elements"]:
        issues.append("atom count, species or order changed")
        return report
    cell = np.asarray(final_cell, dtype=float)
    if cell.shape != (3, 3) or not np.allclose(cell, device["cell"], rtol=0, atol=1.e-5):
        issues.append("relaxed cell does not match fixed device cell")
        return report
    if not np.allclose(cell[:2], electrode["cell"][:2], rtol=0, atol=1.e-5):
        issues.append("transverse cell does not match electrode")
        return report
    repeat = [0, 0, electrode["spacing_A"]*electrode["principal_layers"]]
    if not np.allclose(electrode["cell"][2], repeat, rtol=0, atol=1.e-5):
        issues.append("electrode repeat does not close its ABC stack")
        return report
    xyz = np.array([a.coord() for a in final_atoms])
    ref = np.array(device["positions"])
    if not np.isfinite(xyz).all():
        issues.append("nonfinite final coordinates")
        return report
    xyz = np.asarray(align_device_positions(device, final_atoms, cell))
    tol = report["matching_tolerance_A"]
    if np.max(np.linalg.norm((xyz-ref)[np.asarray(device["fixed_ids"])-1], axis=1)) > tol:
        issues.append("fixed atoms moved relative to common periodic reference")
    lead_xyz = np.array(electrode["positions"])
    lead_count = len(lead_xyz)
    plane_size = math.prod(electrode["in_plane"])
    actual_step = xyz[:plane_size] + cell[2] - xyz[-plane_size:]
    bulk_step = lead_xyz[:plane_size] + repeat - lead_xyz[-plane_size:]
    closure = actual_step - bulk_step
    closure[:, :2] -= np.rint(closure[:, :2] @ np.linalg.inv(cell[:2, :2])) @ cell[:2, :2]
    report["closure_error_A"] = float(np.max(np.linalg.norm(closure, axis=1)))
    if report["closure_error_A"] > tol:
        issues.append("periodic boundary closure differs from bulk Au stacking or spacing")
    metal = []
    for side, sign in (("left", 1), ("right", -1)):
        ids = device["electrode_ids"][side]
        start = 1 if side == "left" else len(xyz)-lead_count+1
        contiguous = ids == list(range(start, start+lead_count))
        direction = device["semi_infinite_direction"][side]
        matched = contiguous and direction == ("-a3" if side == "left" else "+a3")
        if matched:
            block = xyz[np.asarray(ids)-1]
            error = np.max(np.linalg.norm(block-block[0] - (lead_xyz-lead_xyz[0]), axis=1))
            matched = bool(error <= tol and [device["elements"][i-1] for i in ids] == electrode["elements"])
        report["electrode_match"][side] = {"matched": matched, "direction": direction}
        if not matched:
            issues.append(f"{side} electrode block does not match")
        meta = device["contacts"][side]
        initial = _contact_state(meta, sign, ref, electrode, cell)
        final = _contact_state(meta, sign, xyz, electrode, cell)
        def neighbor_ids(record):
            return sorted((p["id"], p["image"]) for p in record["neighbors"])
        final.update(initial_site=meta["site"], motif=meta["motif"], anchor_id=meta["anchor_id"],
                     anchor_element=meta["anchor_element"], initial_neighbors=initial["neighbors"],
                     migrated=(final["final_site"] != meta["site"] or neighbor_ids(initial) != neighbor_ids(final)))
        report["contacts"][side] = final
        if not final["neighbors"]:
            issues.append(f"{side} anchor detached")
        slab = [i for group in meta["layer_ids"] for i in group]
        tips = meta["tip_ids"]
        # Each extra layer must remain connected to its backing layer by an Au-Au bond.
        for tip in tips:
            lower = tips[:-1] if meta["motif"] == "pyramid" and tip == tips[-1] else slab
            bonded = _neighbors(xyz[tip-1], lower, xyz, cell, 1.25*np.linalg.norm(electrode["primitive"][0]))
            if not bonded:
                issues.append(f"{side} tip atom {tip} detached from its backing Au layer")
        metal.extend(slab + tips)
    for atom_id in set(metal).intersection(device["movable_ids"]):
        overlap = _neighbors(xyz[atom_id-1], metal, xyz, cell, .65*np.linalg.norm(electrode["primitive"][0]))
        others = [p for p in overlap if p["id"] != atom_id or p["image"] != [0, 0]]
        if others:
            issues.append(f"Au {atom_id} overlaps Au {others[0]['id']}")
    for atom_id in device["mapping"]["original_to_device"].values():
        overlap = _neighbors(xyz[atom_id-1], metal, xyz, cell, report["overlap_cutoff_A"])
        if overlap:
            issues.append(f"molecular atom {atom_id} overlaps Au {overlap[0]['id']}")
    report["passed"] = not issues
    return report
