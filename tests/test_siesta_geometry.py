import json

import numpy as np
import pytest

from smj.junction.geometry import Atom
from smj.siesta.geometry import build_device, build_electrode, inspect_device


def protocol(left="top", right="top", **extra):
    flat = {"top", "bridge", "fcc", "hcp"}

    def contact(value):
        return dict(height=2.35, **({"site": value} if value in flat else {"motif": value}))

    return {"electrode": {"a_lattice": 4.08, "in_plane": [3, 3],
                          "principal_layers": 3, "slab_layers": 5, "movable_layers": 1},
            "contact": {"left": contact(left), "right": contact(right), **extra}}


def junction(cfg):
    # Terminal H are first and last to exercise original -> relaxed device numbering.
    mol = [Atom("H", -1.3, 0, 0), Atom("S", 0, 0, 0), Atom("C", 3, 1, .2),
           Atom("S", 6, 0, 0), Atom("H", 7.3, 0, 0)]
    lead = build_electrode(cfg)
    device = build_device(mol, (2, 4), lead, cfg, strip_h=True)
    return lead, json.loads(json.dumps(device))


def atoms(d):
    return [Atom(el, *xyz) for el, xyz in zip(d["elements"], d["positions"])]


def test_bulk_lead_and_closed_boundaries():
    cfg = protocol()
    lead, d = junction(cfg)
    p = np.asarray(lead["positions"]).reshape(3, 9, 3)
    a1, a2, c = np.asarray(lead["cell"])
    np.testing.assert_allclose(p[1, 0, :2] - p[0, 0, :2], ((a1 + a2)/9)[:2])
    np.testing.assert_allclose(p[2, 0, :2] - p[1, 0, :2], ((a1 + a2)/9)[:2])
    assert c[2] == pytest.approx(3*4.08/np.sqrt(3))
    pos = np.asarray(d["positions"])
    nlead = len(lead["elements"])
    for side, part in (("left", pos[:nlead]), ("right", pos[-nlead:])):
        np.testing.assert_allclose(part-part[0], p.reshape(-1, 3)-p[0, 0], atol=1e-10)
        assert set(d["electrode_ids"][side]) <= set(d["fixed_ids"])
    assert pos[0, 2] + d["cell"][2][2] - pos[-1, 2] == pytest.approx(c[2]/3)
    assert len(d["movable_ids"]) == 18 + 3
    assert sorted(d["fixed_ids"] + d["movable_ids"]) == list(range(1, len(pos)+1))
    assert d["mapping"]["removed_h_ids"] == [1, 5]
    assert d["contacts"]["left"]["anchor_id"] == d["mapping"]["original_to_device"]["2"]
    assert d["contacts"]["right"]["anchor_id"] == d["mapping"]["original_to_device"]["4"]


@pytest.mark.parametrize("left,right,counts", [
    ("top", "bridge", (0, 0)), ("fcc", "hcp", (0, 0)),
    ("adatom", "trimer", (1, 3)), ("pyramid", "top", (4, 0)),
])
def test_contact_templates_and_coordination(left, right, counts):
    lead, d = junction(protocol(left, right))
    report = inspect_device(d, lead, atoms(d), d["cell"])
    assert report["passed"], report["issues"]
    local_sites = {"adatom": "apex", "pyramid": "apex", "trimer": "center"}
    shell = {"top": 1, "bridge": 2, "fcc": 3, "hcp": 3, "apex": 1, "center": 3}
    for side, name, count in zip(("left", "right"), (left, right), counts):
        assert len(d["contacts"][side]["tip_ids"]) == count
        final = report["contacts"][side]
        assert final["initial_site"] == final["final_site"] == local_sites.get(name, name)
        assert final["coordination"] == shell[final["final_site"]]
        assert not final["migrated"]


def test_match_allows_common_translation_but_detects_bad_boundary():
    lead, d = junction(protocol())
    final = [Atom(a.element, a.x+5, a.y-3, a.z+1) for a in atoms(d)]
    assert inspect_device(d, lead, final, d["cell"])["passed"]
    final[-1].z += .12
    assert not inspect_device(d, lead, final, d["cell"])["passed"]
    different = np.array(d["cell"])
    different[0, 0] += .1
    assert not inspect_device(d, lead, atoms(d), different)["passed"]


def test_site_migration_and_detachment():
    lead, d = junction(protocol())
    final = atoms(d)
    anchor = final[d["contacts"]["left"]["anchor_id"]-1]
    anchor.x += lead["cell"][0][0] / 6  # from top to adjacent bridge
    report = inspect_device(d, lead, final, d["cell"])
    assert report["contacts"]["left"]["final_site"] == "bridge"
    assert report["contacts"]["left"]["migrated"]
    anchor.z += 3
    assert not inspect_device(d, lead, final, d["cell"])["passed"]


def test_geometrical_conflicts_fail():
    cfg = protocol(gap=4)
    with pytest.raises(ValueError, match="gap"):
        junction(cfg)
    cfg = protocol()
    cfg["electrode"]["principal_layers"] = 2
    with pytest.raises(ValueError, match="ABC"):
        build_electrode(cfg)
    cfg = protocol()
    cfg["electrode"]["movable_layers"] = 3
    with pytest.raises(ValueError, match="match"):
        junction(cfg)


def test_explicit_gap_resolves_unspecified_heights_without_straining_molecule():
    cfg = protocol("bridge", "pyramid", gap=18)
    del cfg["contact"]["right"]["height"]
    lead, d = junction(cfg)
    lm, rm = (d["contacts"][side] for side in ("left", "right"))
    assert d["gap_A"] == pytest.approx(18)
    assert lm["height_A"] == 2.35
    p = np.array(d["positions"])
    assert np.linalg.norm(p[lm["anchor_id"]-1]-p[rm["anchor_id"]-1]) == pytest.approx(6)
    assert (p[rm["tip_ids"][-1]-1, 2] - p[rm["anchor_id"]-1, 2]) == pytest.approx(rm["height_A"])
    del cfg["contact"]["left"]["height"]
    _, d = junction(cfg)
    assert d["contacts"]["left"]["height_A"] == d["contacts"]["right"]["height_A"]
    cfg["contact"]["gap"] = 5
    with pytest.raises(ValueError, match="gap"):
        junction(cfg)


def test_reject_inconsistent_lead_repeat_and_colliding_tip():
    lead, d = junction(protocol("trimer", "pyramid"))
    lead["cell"][2][2] += .3
    assert not inspect_device(d, lead, atoms(d), d["cell"])["passed"]
    lead["cell"][2][2] -= .3
    final = atoms(d)
    tip1, tip2, _ = d["contacts"]["left"]["tip_ids"]
    final[tip2-1] = final[tip1-1]
    assert not inspect_device(d, lead, final, d["cell"])["passed"]


def test_fcc_hollow_is_above_third_layer_and_hcp_above_second():
    lead, d = junction(protocol("fcc", "hcp"))
    p = np.array(d["positions"])
    cell = np.array(d["cell"])[:2, :2]
    for side, delta in (("left", -3), ("right", 1)):
        meta = d["contacts"][side]
        registry = p[np.array(meta["layer_ids"][delta])-1, :2]
        frac = (p[meta["anchor_id"]-1, :2] - registry) @ np.linalg.inv(cell)
        assert np.min(np.linalg.norm(frac-np.rint(frac), axis=1)) < 1.e-10


def test_tip_constraints_and_removed_hydrogens_are_traceable_by_side():
    cfg = protocol("trimer", "pyramid")
    cfg["contact"]["right"]["movable_tip_atoms"] = [4]
    _, d = junction(cfg)
    lm, rm = (d["contacts"][side] for side in ("left", "right"))
    assert lm["removed_original_h_ids"] == [1] and rm["removed_original_h_ids"] == [5]
    assert set(rm["tip_ids"][:3]) <= set(d["fixed_ids"])
    assert rm["tip_ids"][3] in d["movable_ids"]
    assert set(lm["tip_ids"]) <= set(d["movable_ids"])


def test_independent_periodic_closure_check():
    lead, d = junction(protocol())
    d["cell"][2][2] += .3
    result = inspect_device(d, lead, atoms(d), d["cell"])
    assert not result["passed"]
    assert any("closure" in issue for issue in result["issues"])


def test_default_initialization_and_reject_misspelled_geometry_keys():
    _, d = junction({})
    for meta in d["contacts"].values():
        assert (meta["motif"], meta["site"], meta["tip_ids"]) == ("flat", "top", [])
    for cfg in ({"electrode": {"principal_layer": 6}}, {"contact": {"gapp": 7}},
                {"contact": {"left": {"sites": "fcc"}}}):
        with pytest.raises(ValueError, match="unknown"):
            junction(cfg)


def test_reconstruct_transiesta_electrode_coordinates_after_xv_image_wrapping():
    from smj.siesta.geometry import align_device_positions

    lead, d = junction(protocol())
    final = atoms(d)
    final[-1].z += d["cell"][2][2]
    aligned = align_device_positions(d, final, d["cell"])
    np.testing.assert_allclose(aligned, d["positions"], atol=1e-7)
    assert inspect_device(d, lead, final, d["cell"])["passed"]
