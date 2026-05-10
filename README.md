# Operadic Mechanism Catalog

Recasting the Sigrist-Ueda symmetry-imposed pairing table as a colored
operad of pairing mechanisms, where each operadic operation is a pairing
channel and each material defines a representation (algebra over the
operad).

## What this is

Sigrist and Ueda (Rev. Mod. Phys., 1991) tabulated, per crystallographic
point group, the irreps available for the superconducting order parameter.
That table is discrete: no composition, no family relationships, no
unified structure across the 32 point groups.

This project re-casts that table as an **operad**, capturing:

- **Compositional structure**: channel composition via colored-operad
  insertion.
- **Functorial materials**: each material is a functor (algebra) from the
  operad to vector spaces.
- **Parametric stability**: families related by symmetry reduction give
  compatible operad representations.

## Structure

```
python/
  operadic_catalog/        Python implementation
    operad.py              colored operad definition + axiom checks
    generators.py          generators per crystallographic point group
    algebra.py             material -> O-algebra functor
    equivalence.py         operadic equivalence classes
    morphisms.py           symmetry-reduction morphisms
    sigrist_ueda.py        v14 channel data + spacegroup mapping
    __init__.py
  tests/                   48 tests, all passing
    test_operad_axioms.py
    test_generators.py
    test_algebra.py
    test_equivalence.py

lean/                      Lean 4 formal verification
  OperadicCatalog.lean     all theorems verified with Mathlib
  lakefile.lean            Lake project config
  lean-toolchain           leanprover/lean4:v4.30.0-rc2

latex/                     Paper draft
  main.tex                 main paper
  supplementary.tex        supplementary material
  refs.bib                 bibliography
  figs/                    figures

data/                      data files (if needed)
```

## Verified theorems (Lean 4)

All theorems are formally verified in `lean/OperadicCatalog.lean`:

| Theorem | Statement |
|---------|-----------|
| `s_wave_universal` | s-wave is allowed in all 32 point groups |
| `unit_right` | f circ id = f |
| `unit_left` | id circ f = id |
| `comp_assoc` | (f circ g) circ h = f circ (g circ h) |
| `restrict_unit` | restriction preserves identity |
| `restrict_comp` | restriction preserves composition (equivariance) |
| `sigrist_ueda_is_free_algebra` | Sigrist-Ueda = free O-algebra |
| `card_eq` | exactly 32 crystallographic point groups |
| `d4h_channel_count` | D4h admits exactly 6 channels |
| `oh_channel_count` | Oh admits exactly 5 channels |
| `c1_channel_count` | C1 admits exactly 1 channel |
| `maxSub_irrefl` | maximal-subgroup relation is irreflexive |
| `isSubgroup_antisymm` | subgroup relation is antisymmetric |

Axiom footprint: `{propext, Classical.choice, Quot.sound}` only.

## Context

Task T1.3 in the RTSC post-v14 roadmap.

## Dependencies

**Python**: Python 3.10+, numpy, networkx

**Lean**: Lean 4.30.0-rc2, Mathlib (rev 7bdf4031)

## References

- Sigrist, Ueda. "Phenomenological theory of unconventional
  superconductivity." Rev. Mod. Phys. 63, 239 (1991).
- Loday, Vallette. "Algebraic Operads." Springer (2012).
- Markl, Shnider, Stasheff. "Operads in Algebra, Topology and Physics."
  AMS (2002).
