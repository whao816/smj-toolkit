"""Read the text and sequential Fortran records written by SIESTA 5.4.2."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re
import struct

import numpy as np

from ..checks import ATOMIC_NUMBERS
from ..errors import StageError
from ..junction.geometry import Atom

BOHR = 0.529177210903
RY_EV = 13.605693122994
VERSION = "5.4.2"
NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dump(path, data):
    path = Path(path)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")
    return path


def read_json(path):
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            raise ValueError("expected object")
        return data
    except (OSError, ValueError) as exc:
        raise StageError(f"{path}: {exc}") from exc


def verify_files(directory, manifest):
    files = manifest.get("files")
    if not isinstance(files, dict) or not files:
        raise StageError(f"{directory}: missing file digests")
    for name, expected in files.items():
        if Path(name).name != name or not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
            raise StageError(f"invalid manifest filename {name!r}")
        path = Path(directory) / name
        if not path.is_file() or sha(path) != expected:
            raise StageError(f"missing or changed file (sha256): {path}")


def convergence(log, *, open_system=False):
    version = re.search(r"(?im)^\s*(?:Siesta\s+)?Version\s*:\s*(\S+)", log)
    if not version or version[1] != VERSION:
        raise StageError(f"expected SIESTA {VERSION} version in output")
    if not re.search(r"(?im)^\s*>>\s*End of run:", log):
        raise StageError("SIESTA output missing End of run")
    if re.search(r"SCF_NOT_CONV|GEOM_NOT_CONV|OUT_OF_TIME|^\s*(?:FATAL|ERROR):", log, re.M | re.I):
        raise StageError("SIESTA reported failure or nonconvergence")
    key = "ts-scf" if open_system else "scf"
    iterations = list(re.finditer(rf"(?im)^\s*{key}:\s*(\d+)\s+[^\n]+", log))
    if not iterations:
        raise StageError(f"missing {key} iteration")
    tail = log[iterations[-1].end():]
    if not re.search(r"(?im)^\s*SCF cycle converged after\s+\d+\s+iterations?", tail):
        raise StageError(f"{key} final cycle did not converge")
    moment = re.findall(rf"(?im)^\s*spin moment:[^\n]*\}}\s*({NUM})", log)
    return {"version": version[1], "scf_converged": True, "scf_iterations": int(iterations[-1][1]),
            "final_moment_muB": float(moment[-1]) if moment else None}


def open_charge(log, spin):
    """In ts_charge.F90 the collinear output is (up, down, sum)."""
    pattern = (rf"(?im)^\s*transiesta:\s*Charge distribution, target\s*=\s*({NUM})"
               r"\s+Total charge\s*\[Q\]\s*:\s*([^\r\n]+)")
    matches = re.findall(pattern, log)
    if not matches:
        raise StageError("transiesta_scf: missing final charge distribution")
    target, row = matches[-1]
    fields = row.split()
    if (len(fields) != (3 if spin == "collinear" else 1)
            or not all(re.fullmatch(NUM, v) for v in fields)):
        raise StageError("transiesta_scf: invalid spin-resolved charge distribution")
    return float(target), float(fields[-1])


def mulliken_populations(log):
    """Return the last complete ``Mulliken Atomic Populations`` table.

    SIESTA prints the table once per geometry/end event; the last one belongs to
    the converged charge density.  Rows are ``(index, charge_e, valence_e,
    species)`` with one-based indices in device atom order.
    """
    marker = "Mulliken Atomic Populations:"
    parsed = None
    for block in log.split(marker)[1:]:
        rows = []
        started = False
        for line in block.splitlines():
            if not started:
                started = bool(re.match(r"\s*Atom #\s+charge \[q\]", line))
                continue
            if re.match(r"^-{10,}\s*$", line):
                if rows:
                    parsed = rows
                break
            match = re.fullmatch(rf"\s*(\d+)\s+({NUM})\s+({NUM})\s+([A-Za-z]{{1,3}})\s*", line)
            if match:
                rows.append((int(match[1]), float(match[2]), float(match[3]), match[4]))
    if parsed is None:
        raise StageError("SIESTA output missing a complete Mulliken Atomic Populations table")
    return parsed


def relaxed_geometry(path, expected):
    """XV rows contain species and atomic number, with cell/coordinates in Bohr."""
    try:
        lines = Path(path).read_text(encoding="ascii").splitlines()
        cell = np.array([list(map(float, line.split()[:3])) for line in lines[:3]]) * BOHR
        count = int(lines[3].strip())
        if count != len(expected["elements"]) or len(lines[4:]) != count:
            raise ValueError("XV atom count")
        if not np.allclose(cell, expected["cell"], rtol=2e-6, atol=1e-5):
            raise ValueError("XV cell changed")
        species = sorted(set(expected["elements"]))
        atoms = []
        for symbol, line in zip(expected["elements"], lines[4:]):
            sid, z, *fields = line.split()
            if int(z) != ATOMIC_NUMBERS[symbol] or int(sid) != species.index(symbol) + 1:
                raise ValueError("XV atom ordering")
            xyz = np.array(list(map(float, fields[:3]))) * BOHR
            if xyz.shape != (3,) or not np.isfinite(xyz).all():
                raise ValueError("XV coordinates")
            atoms.append(Atom(symbol, *xyz))
        return atoms, cell.tolist()
    except (OSError, IndexError, ValueError) as exc:
        raise StageError(f"{path}: {exc}") from exc


def _records(stream, size, first_length):
    endian = "<" if stream.read(4) == struct.pack("<i", first_length) else ">"
    stream.seek(0)

    def record(expected=None):
        n, = struct.unpack(endian + "i", stream.read(4))
        if n <= 0 or n > size - stream.tell() - 4 or (expected is not None and n != expected):
            raise ValueError(f"invalid record length {n} (expected {expected})")
        data = stream.read(n)
        if struct.unpack(endian + "i", stream.read(4))[0] != n:
            raise ValueError("incomplete Fortran record")
        return data
    return endian, record


def inspect_restart(path, hsx, spin):
    """Check the DM/energy-DM records written by Src/m_ts_iodm.F90."""
    try:
        with path.open("rb") as stream:
            endian, record = _records(stream, path.stat().st_size, 20)
            no, ns, *images = struct.unpack(endian + "5i", record(20))
            if (no != hsx["orbitals"] or ns != (2 if spin == "collinear" else 1)
                    or min(images) < 1):
                raise ValueError("incompatible dimensions")
            counts = np.frombuffer(record(no * 4), dtype=endian + "i4")
            maxcol = no * math.prod(images)
            if (counts < 1).any() or (counts > maxcol).any():
                raise ValueError("invalid row lengths")
            for n in counts:
                cols = np.frombuffer(record(int(n) * 4), dtype=endian + "i4")
                if (cols < 1).any() or (cols > maxcol).any():
                    raise ValueError("invalid orbital indices")
            for _ in range(2 * ns):
                for n in counts:
                    row = np.frombuffer(record(int(n) * 8), dtype=endian + "f8")
                    if not np.isfinite(row).all():
                        raise ValueError("nonfinite density matrix")
            ef, = struct.unpack(endian + "d", record(8))
            if not math.isfinite(ef) or stream.read(1):
                raise ValueError("invalid final TSDE record")
    except (OSError, ValueError, struct.error) as exc:
        raise StageError(f"invalid TSDE {path}: {exc}") from exc


def inspect_hsx(path, geometry, spin, *, electrode=False):
    """Verify v2 records against Src/io_hsx.F90 and Src/io_sparse.F90.

    Read sparsity, Hamiltonian and overlap one orbital row at a time; never
    materialise a complete matrix. Cell translations by a lattice vector are
    equivalent for coordinates returned by SIESTA.
    """
    path = Path(path)
    try:
        with path.open("rb") as f:
            endian, record = _records(f, path.stat().st_size, 4)

            def ints(data):
                return np.frombuffer(data, dtype=endian + "i4")

            def doubles(data):
                return np.frombuffer(data, dtype=endian + "f8")

            if ints(record(4))[0] != 2:
                raise ValueError("expected v2 HSX")
            real_width = 8 if ints(record(4))[0] else 4
            na, no, ns, nsp, *nsc = map(int, ints(record(28)))
            symbols = geometry["elements"]
            if (na != len(symbols) or no < na or ns != (2 if spin == "collinear" else 1)
                    or nsp != len(set(symbols)) or min(nsc) < 1):
                raise ValueError("incompatible atom/orbital/species/spin count")
            meta = doubles(record(96))
            cell = meta[:9].reshape(3, 3) * BOHR
            if not np.isfinite(meta).all() or not np.allclose(cell, geometry["cell"], atol=1e-5, rtol=2e-6):
                raise ValueError("incompatible unit cell or reference energy")
            nimages = math.prod(nsc)
            block = record(nimages * 12 + na * 32)
            offsets = ints(block[:nimages*12]).reshape(-1, 3)
            boundary = nimages*12 + na*24
            positions = doubles(block[nimages*12:boundary]).reshape(-1, 3) * BOHR
            ids = ints(block[boundary:boundary+na*4])
            delta = (positions - np.asarray(geometry["positions"])) @ np.linalg.inv(cell)
            if not np.allclose((delta - delta.round()) @ cell, 0., atol=2e-4):
                raise ValueError("HSX coordinates differ from input")
            labels = [label.decode("ascii").strip() for label, _, _ in
                      struct.iter_unpack(endian + "20sdi", record(nsp * 32))]
            if any(i <= 0 or i > nsp for i in ids) or [labels[i-1] for i in ids] != symbols:
                raise ValueError("HSX atom labels/order differ from input")
            for _ in range(nsp):
                record()
            record(60)  # k-cell and displacements
            counts = ints(record(no * 4))
            if (counts <= 0).any() or (counts > no * nimages).any():
                raise ValueError("invalid sparse row lengths")
            nearest_cells = set()
            for n in counts:
                cols = ints(record(4 * int(n))).astype(np.int64) - 1
                if (cols < 0).any() or (cols >= no * nimages).any():
                    raise ValueError("invalid orbital index")
                if electrode:
                    z = offsets[cols // no, 2]
                    if (np.abs(z) > 1).any():
                        raise ValueError("couplings extend beyond neighbouring principal layers")
                    nearest_cells.update(np.unique(z))
            if electrode and not {-1, 1}.issubset(nearest_cells):
                raise ValueError("uncoupled electrode principal layer along transport axis")
            for _ in range(ns + 1):
                for n in counts:
                    row = np.frombuffer(record(real_width * int(n)), dtype=endian + f"f{real_width}")
                    if not np.isfinite(row).all():
                        raise ValueError("nonfinite H/S matrix")
            if f.read(1):
                raise ValueError("unexpected trailing HSX data")
            return {"fermi_ev": float(meta[9] * RY_EV), "orbitals": no, "hsx_version": 2}
    except (ValueError, OSError, struct.error, np.linalg.LinAlgError) as exc:
        raise StageError(f"invalid HSX {path}: {exc}") from exc
