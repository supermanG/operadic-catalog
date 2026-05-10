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
src/
  operadic_catalog/
    operad.py          -- colored operad definition + axiom checks
    generators.py      -- generators per crystallographic point group
    algebra.py         -- material -> O-algebra functor
    equivalence.py     -- operadic equivalence classes
    morphisms.py       -- symmetry-reduction morphisms
    sigrist_ueda.py    -- v14 channel data + spacegroup mapping
    __init__.py
tests/
  test_operad_axioms.py
  test_generators.py
  test_algebra.py
  test_equivalence.py
paper/
  (LaTeX draft, later)
```

## Context

Task T1.3 in the RTSC post-v14 roadmap. See
`roadmap/HANDOFF_operadic_mechanism_catalog.md` in the parent project
for the full mathematical specification.

## Dependencies

- Python 3.10+
- numpy
- networkx (for morphism graphs)

## References

- Sigrist, Ueda. "Phenomenological theory of unconventional
  superconductivity." Rev. Mod. Phys. 63, 239 (1991).
- Loday, Vallette. "Algebraic Operads." Springer (2012).
- Markl, Shnider, Stasheff. "Operads in Algebra, Topology and Physics."
  AMS (2002).
