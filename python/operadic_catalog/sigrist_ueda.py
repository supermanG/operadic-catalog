"""Sigrist-Ueda channel table and spacegroup-to-point-group mapping.

Data source: Sigrist and Ueda, Rev. Mod. Phys. 63, 239 (1991).
Adapted from v14 (run_v14.py) and v6 (run_v6.py) of the RTSC pipeline.

Channel labels:
  's':     s-wave (trivial irrep, always present)
  'p':     p-wave (odd-parity, requires inversion)
  'd':     d-wave (any non-trivial 2D irrep or B-type irrep)
  'd_x2y2': d_{x^2-y^2} specifically (B_1g of D_4h-type)
  'd_xy':  d_xy specifically (B_2g of D_4h-type)
  'f':     f-wave (3D irreps in cubic, T_1u / T_2u)
  'z3':    Z/3 trion (E or non-trivial 1D complex irrep of C_3)
  'z4':    Z/4 (B-irreps of C_4 / D_4 / D_4h-type)
  'z6':    Z/6 (E-irreps of C_6)
"""

from __future__ import annotations

POINT_GROUP_CHANNELS: dict[str, frozenset[str]] = {
    "C1":  frozenset({"s"}),
    "Ci":  frozenset({"s", "p"}),
    "C2":  frozenset({"s", "d"}),
    "Cs":  frozenset({"s", "d"}),
    "C2h": frozenset({"s", "p", "d"}),
    "D2":  frozenset({"s", "d"}),
    "C2v": frozenset({"s", "d"}),
    "D2h": frozenset({"s", "p", "d"}),
    "C3":  frozenset({"s", "z3"}),
    "S6":  frozenset({"s", "p", "z3"}),
    "D3":  frozenset({"s", "z3"}),
    "C3v": frozenset({"s", "z3"}),
    "D3d": frozenset({"s", "p", "z3"}),
    "C4":  frozenset({"s", "d", "z4"}),
    "S4":  frozenset({"s", "d", "z4"}),
    "C4h": frozenset({"s", "p", "d", "z4"}),
    "D4":  frozenset({"s", "d", "d_x2y2", "d_xy", "z4"}),
    "C4v": frozenset({"s", "d", "d_x2y2", "d_xy", "z4"}),
    "D2d": frozenset({"s", "d", "d_x2y2", "d_xy", "z4"}),
    "D4h": frozenset({"s", "p", "d", "d_x2y2", "d_xy", "z4"}),
    "C6":  frozenset({"s", "d", "z3", "z6"}),
    "C3h": frozenset({"s", "d", "z3", "z6"}),
    "C6h": frozenset({"s", "p", "d", "z3", "z6"}),
    "D6":  frozenset({"s", "d", "z3", "z6"}),
    "C6v": frozenset({"s", "d", "z3", "z6"}),
    "D3h": frozenset({"s", "d", "z3"}),
    "D6h": frozenset({"s", "p", "d", "z3", "z6"}),
    "T":   frozenset({"s", "d", "z3", "f"}),
    "Th":  frozenset({"s", "p", "d", "z3", "f"}),
    "O":   frozenset({"s", "d", "z3", "f"}),
    "Td":  frozenset({"s", "d", "z3", "f"}),
    "Oh":  frozenset({"s", "p", "d", "z3", "f"}),
}

ALL_CHANNELS = frozenset({"s", "p", "d", "d_x2y2", "d_xy", "f", "z3", "z4", "z6"})
ALL_POINT_GROUPS = frozenset(POINT_GROUP_CHANNELS.keys())

CHANNEL_CYCLIC_ORDER: dict[str, int] = {
    "s": 1,
    "p": 2,
    "d": 2,
    "d_x2y2": 2,
    "d_xy": 2,
    "f": 3,
    "z3": 3,
    "z4": 4,
    "z6": 6,
}

SCHUR_MULTIPLIER_ORDER: dict[str, int] = {
    "C1": 1, "C2": 1, "C3": 1, "C4": 1, "C6": 1,
    "D2": 2, "D3": 1, "D4": 2, "D6": 2,
    "T":  2, "O":  2,
    "Ci":  1, "Cs":  1, "C2h": 1, "C4h": 1, "C6h": 1,
    "C3h": 1, "S4":  1, "S6":  1,
    "C2v": 2, "C4v": 2, "C6v": 2, "C3v": 1,
    "D2d": 2, "D3d": 2, "D3h": 2,
    "D2h": 8, "D4h": 8, "D6h": 8,
    "Th":  2, "Td":  2, "Oh":  4,
}


def spg_to_point_group(spg: int | None) -> str | None:
    """Map spacegroup number (1-230) to one of the 32 crystallographic
    point groups (Schoenflies notation)."""
    if spg is None:
        return None
    s = int(spg)
    if s == 1:               return "C1"
    if s == 2:               return "Ci"
    if 3 <= s <= 5:          return "C2"
    if 6 <= s <= 9:          return "Cs"
    if 10 <= s <= 15:        return "C2h"
    if 16 <= s <= 24:        return "D2"
    if 25 <= s <= 46:        return "C2v"
    if 47 <= s <= 74:        return "D2h"
    if 75 <= s <= 80:        return "C4"
    if 81 <= s <= 82:        return "S4"
    if 83 <= s <= 88:        return "C4h"
    if 89 <= s <= 98:        return "D4"
    if 99 <= s <= 110:       return "C4v"
    if 111 <= s <= 122:      return "D2d"
    if 123 <= s <= 142:      return "D4h"
    if 143 <= s <= 146:      return "C3"
    if 147 <= s <= 148:      return "S6"
    if 149 <= s <= 155:      return "D3"
    if 156 <= s <= 161:      return "C3v"
    if 162 <= s <= 167:      return "D3d"
    if 168 <= s <= 173:      return "C6"
    if s == 174:             return "C3h"
    if 175 <= s <= 176:      return "C6h"
    if 177 <= s <= 182:      return "D6"
    if 183 <= s <= 186:      return "C6v"
    if 187 <= s <= 190:      return "D3h"
    if 191 <= s <= 194:      return "D6h"
    if 195 <= s <= 199:      return "T"
    if 200 <= s <= 206:      return "Th"
    if 207 <= s <= 214:      return "O"
    if 215 <= s <= 220:      return "Td"
    if 221 <= s <= 230:      return "Oh"
    return None


def allowed_channels(point_group: str | None) -> frozenset[str]:
    """Return channels allowed by symmetry for the given point group."""
    if point_group is None:
        return frozenset()
    return POINT_GROUP_CHANNELS.get(point_group, frozenset({"s"}))


def h2_anomaly_class(point_group: str | None) -> str:
    """Return the H^2 anomaly class: 'trivial', 'Z2', 'Z2^2', or 'Z2^3'."""
    if point_group is None:
        return "unknown"
    m = SCHUR_MULTIPLIER_ORDER.get(point_group, 1)
    if m == 1:
        return "trivial"
    if m == 2:
        return "Z2"
    if m == 4:
        return "Z2^2"
    if m == 8:
        return "Z2^3"
    return f"order_{m}"
