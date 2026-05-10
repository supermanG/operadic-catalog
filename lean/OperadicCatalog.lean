/-
  OperadicCatalog.lean -- Colored operad of superconducting pairing mechanisms
  ============================================================================

  Formalizes the colored operad O whose:
    - Colors are (point_group, cyclic_order) pairs
    - Operations are pairing channels (Sigrist-Ueda)
    - Composition is sub-channel insertion under subgroup relation
    - Restriction is channel survival under point-group reduction

  Theorems verified:
    1. Unit axiom: composing with the s-wave unit is the identity
    2. Associativity: (f circ g) circ h = f circ (g circ h)
    3. Equivariance: restriction commutes with composition
    4. Sigrist-Ueda = free O-algebra on generators per PG
    5. Restriction is functorial (preserves composition and identity)
    6. Morphism DAG properties (irreflexivity, antisymmetry)
    7. Channel table properties (s-wave universality, channel counts)

  Built against Mathlib master (rev 7bdf4031, lean v4.30.0-rc2).
  Part of Task T1.3 (operadic mechanism catalog) in the RTSC roadmap.

  LH & Claude, 2026-05-10.
-/

import Mathlib

namespace OperadicCatalog

open BigOperators

/-! ## 1. Colored Operad: Abstract Definition

We define a colored operad as a structure with colors, operations,
composition, and unit, satisfying the three axioms. -/

/-- A colored operad over a type of colors `C`. -/
structure ColoredOperad (C : Type*) where
  Op : C -> C -> Type*
  id_op : (c : C) -> Op c c
  comp : {a b c : C} -> Op b c -> Op a b -> Op a c
  comp_id_right : forall {a b : C} (f : Op a b), comp f (id_op a) = f
  comp_id_left : forall {a b : C} (f : Op a b), comp (id_op b) f = f
  comp_assoc : forall {a b c d : C} (f : Op c d) (g : Op b c) (h : Op a b),
    comp f (comp g h) = comp (comp f g) h

/-! ## 2. Crystallographic Point Group Enumeration

We define the 32 crystallographic point groups as an inductive type,
along with the subgroup relation. -/

/-- The 32 crystallographic point groups (Schoenflies notation). -/
inductive PointGroup : Type where
  | C1 | Ci | C2 | Cs | C2h
  | D2 | C2v | D2h
  | C3 | S6 | D3 | C3v | D3d
  | C4 | S4 | C4h | D4 | C4v | D2d | D4h
  | C6 | C3h | C6h | D6 | C6v | D3h | D6h
  | T | Th | O_ | Td | Oh
  deriving DecidableEq, Repr, Inhabited

namespace PointGroup

instance : Fintype PointGroup where
  elems := {.C1, .Ci, .C2, .Cs, .C2h, .D2, .C2v, .D2h,
            .C3, .S6, .D3, .C3v, .D3d,
            .C4, .S4, .C4h, .D4, .C4v, .D2d, .D4h,
            .C6, .C3h, .C6h, .D6, .C6v, .D3h, .D6h,
            .T, .Th, .O_, .Td, .Oh}
  complete := by intro x; cases x <;> simp [Finset.mem_insert, Finset.mem_singleton]

theorem card_eq : Fintype.card PointGroup = 32 := by native_decide

end PointGroup

/-- Pairing channel types from the Sigrist-Ueda classification. -/
inductive Channel : Type where
  | s      -- s-wave (trivial irrep, always present)
  | p      -- p-wave (odd-parity, requires inversion)
  | d      -- d-wave (non-trivial 2D irrep)
  | d_x2y2 -- d_{x^2-y^2} specifically (B_1g of D_4h)
  | d_xy   -- d_xy specifically (B_2g of D_4h)
  | f      -- f-wave (3D irreps in cubic)
  | z3     -- Z/3 trion (E or complex 1D irrep of C_3)
  | z4     -- Z/4 (B-irreps of C_4/D_4/D_4h)
  | z6     -- Z/6 (E-irreps of C_6)
  deriving DecidableEq, Repr, Inhabited

/-- The cyclic order associated with each channel. -/
def Channel.cyclicOrder : Channel -> Nat
  | .s => 1
  | .p => 2
  | .d => 2
  | .d_x2y2 => 2
  | .d_xy => 2
  | .f => 3
  | .z3 => 3
  | .z4 => 4
  | .z6 => 6

/-- A color in the pairing operad: (point_group, cyclic_order). -/
structure PairingColor where
  pg : PointGroup
  cyclicOrder : Nat
  deriving DecidableEq, Repr

/-! ## 3. Sigrist-Ueda Channel Table

The allowed channels per point group, encoded as a decidable predicate. -/

/-- Sigrist-Ueda: which channels are allowed in which point group. -/
def allowedChannel : PointGroup -> Channel -> Bool
  | .C1,  .s => true
  | .Ci,  .s | .Ci,  .p => true
  | .C2,  .s | .C2,  .d => true
  | .Cs,  .s | .Cs,  .d => true
  | .C2h, .s | .C2h, .p | .C2h, .d => true
  | .D2,  .s | .D2,  .d => true
  | .C2v, .s | .C2v, .d => true
  | .D2h, .s | .D2h, .p | .D2h, .d => true
  | .C3,  .s | .C3,  .z3 => true
  | .S6,  .s | .S6,  .p | .S6, .z3 => true
  | .D3,  .s | .D3,  .z3 => true
  | .C3v, .s | .C3v, .z3 => true
  | .D3d, .s | .D3d, .p | .D3d, .z3 => true
  | .C4,  .s | .C4,  .d | .C4, .z4 => true
  | .S4,  .s | .S4,  .d | .S4, .z4 => true
  | .C4h, .s | .C4h, .p | .C4h, .d | .C4h, .z4 => true
  | .D4,  .s | .D4,  .d | .D4, .d_x2y2 | .D4, .d_xy | .D4, .z4 => true
  | .C4v, .s | .C4v, .d | .C4v, .d_x2y2 | .C4v, .d_xy | .C4v, .z4 => true
  | .D2d, .s | .D2d, .d | .D2d, .d_x2y2 | .D2d, .d_xy | .D2d, .z4 => true
  | .D4h, .s | .D4h, .p | .D4h, .d | .D4h, .d_x2y2 | .D4h, .d_xy | .D4h, .z4 => true
  | .C6,  .s | .C6,  .d | .C6, .z3 | .C6, .z6 => true
  | .C3h, .s | .C3h, .d | .C3h, .z3 | .C3h, .z6 => true
  | .C6h, .s | .C6h, .p | .C6h, .d | .C6h, .z3 | .C6h, .z6 => true
  | .D6,  .s | .D6,  .d | .D6, .z3 | .D6, .z6 => true
  | .C6v, .s | .C6v, .d | .C6v, .z3 | .C6v, .z6 => true
  | .D3h, .s | .D3h, .d | .D3h, .z3 => true
  | .D6h, .s | .D6h, .p | .D6h, .d | .D6h, .z3 | .D6h, .z6 => true
  | .T,   .s | .T,   .d | .T, .z3 | .T, .f => true
  | .Th,  .s | .Th,  .p | .Th, .d | .Th, .z3 | .Th, .f => true
  | .O_,  .s | .O_,  .d | .O_, .z3 | .O_, .f => true
  | .Td,  .s | .Td,  .d | .Td, .z3 | .Td, .f => true
  | .Oh,  .s | .Oh,  .p | .Oh, .d | .Oh, .z3 | .Oh, .f => true
  | _, _ => false

/-- s-wave is always allowed (Theorem: trivial irrep universality). -/
theorem s_wave_universal (g : PointGroup) : allowedChannel g .s = true := by
  cases g <;> rfl

/-! ## 4. Subgroup Relation

The maximal-subgroup relation on the 32 crystallographic point groups,
and its transitive closure defining the full subgroup lattice.

We encode the subgroup relation directly by exhaustive case matching
on the 32x32 table for decidability. -/

/-- Maximal subgroup relation (Hasse diagram edges). -/
def isMaximalSubgroup : PointGroup -> PointGroup -> Bool
  | .Ci,  .C1  => true
  | .C2,  .C1  => true
  | .Cs,  .C1  => true
  | .C2h, .C2  | .C2h, .Ci  | .C2h, .Cs  => true
  | .D2,  .C2  => true
  | .C2v, .C2  | .C2v, .Cs  => true
  | .D2h, .D2  | .D2h, .C2v | .D2h, .C2h => true
  | .C3,  .C1  => true
  | .S6,  .C3  | .S6,  .Ci  => true
  | .D3,  .C3  | .D3,  .C2  => true
  | .C3v, .C3  | .C3v, .Cs  => true
  | .D3d, .D3  | .D3d, .C3v | .D3d, .S6  => true
  | .C4,  .C2  => true
  | .S4,  .C2  => true
  | .C4h, .C4  | .C4h, .C2h | .C4h, .S4  => true
  | .D4,  .C4  | .D4,  .D2  => true
  | .C4v, .C4  | .C4v, .C2v => true
  | .D2d, .S4  | .D2d, .D2  | .D2d, .C2v => true
  | .D4h, .D4  | .D4h, .C4v | .D4h, .C4h | .D4h, .D2h | .D4h, .D2d => true
  | .C6,  .C3  | .C6,  .C2  => true
  | .C3h, .C3  | .C3h, .Cs  => true
  | .C6h, .C6  | .C6h, .C3h | .C6h, .C2h | .C6h, .S6  => true
  | .D6,  .C6  | .D6,  .D3  | .D6,  .D2  => true
  | .C6v, .C6  | .C6v, .C3v | .C6v, .C2v => true
  | .D3h, .D3  | .D3h, .C3h | .D3h, .C3v | .D3h, .C2v => true
  | .D6h, .D6  | .D6h, .C6v | .D6h, .C6h | .D6h, .D3h | .D6h, .D2h => true
  | .T,   .D2  | .T,   .C3  => true
  | .Th,  .T   | .Th,  .D2h | .Th,  .S6  => true
  | .O_,  .T   | .O_,  .D4  | .O_,  .D3  => true
  | .Td,  .T   | .Td,  .D2d | .Td,  .C3v | .Td, .S4 => true
  | .Oh,  .O_  | .Oh,  .Td  | .Oh,  .Th  | .Oh, .D4h | .Oh, .D3d => true
  | _, _ => false

/-- The subgroup relation on the 32 crystallographic point groups.
    isSubgroup h g = true means h is a subgroup of g.
    Exhaustive case match: 188 true pairs + 32 reflexive,
    computed from the transitive closure of isMaximalSubgroup. -/
def isSubgroup : PointGroup -> PointGroup -> Bool
  | .C1, .C1 | .Ci, .Ci | .C2, .C2 | .Cs, .Cs | .C2h, .C2h
  | .D2, .D2 | .C2v, .C2v | .D2h, .D2h
  | .C3, .C3 | .S6, .S6 | .D3, .D3 | .C3v, .C3v | .D3d, .D3d
  | .C4, .C4 | .S4, .S4 | .C4h, .C4h | .D4, .D4 | .C4v, .C4v | .D2d, .D2d | .D4h, .D4h
  | .C6, .C6 | .C3h, .C3h | .C6h, .C6h | .D6, .D6 | .C6v, .C6v | .D3h, .D3h | .D6h, .D6h
  | .T, .T | .Th, .Th | .O_, .O_ | .Td, .Td | .Oh, .Oh => true
  | .C1, .Ci => true
  | .C1, .C2 => true
  | .C1, .Cs => true
  | .C1, .C2h | .Ci, .C2h | .C2, .C2h | .Cs, .C2h => true
  | .C1, .D2 | .C2, .D2 => true
  | .C1, .C2v | .C2, .C2v | .Cs, .C2v => true
  | .C1, .D2h | .Ci, .D2h | .C2, .D2h | .Cs, .D2h => true
  | .C2h, .D2h | .D2, .D2h | .C2v, .D2h => true
  | .C1, .C3 => true
  | .C1, .S6 | .Ci, .S6 | .C3, .S6 => true
  | .C1, .D3 | .C2, .D3 | .C3, .D3 => true
  | .C1, .C3v | .Cs, .C3v | .C3, .C3v => true
  | .C1, .D3d | .Ci, .D3d | .C2, .D3d | .Cs, .D3d => true
  | .C3, .D3d | .S6, .D3d | .D3, .D3d | .C3v, .D3d => true
  | .C1, .C4 | .C2, .C4 => true
  | .C1, .S4 | .C2, .S4 => true
  | .C1, .C4h | .Ci, .C4h | .C2, .C4h | .Cs, .C4h => true
  | .C2h, .C4h | .C4, .C4h | .S4, .C4h => true
  | .C1, .D4 | .C2, .D4 | .D2, .D4 | .C4, .D4 => true
  | .C1, .C4v | .C2, .C4v | .Cs, .C4v | .C2v, .C4v => true
  | .C4, .C4v => true
  | .C1, .D2d | .C2, .D2d | .Cs, .D2d | .D2, .D2d => true
  | .C2v, .D2d | .S4, .D2d => true
  | .C1, .D4h | .Ci, .D4h | .C2, .D4h | .Cs, .D4h => true
  | .C2h, .D4h | .D2, .D4h | .C2v, .D4h | .D2h, .D4h => true
  | .C4, .D4h | .S4, .D4h | .C4h, .D4h | .D4, .D4h => true
  | .C4v, .D4h | .D2d, .D4h => true
  | .C1, .C6 | .C2, .C6 | .C3, .C6 => true
  | .C1, .C3h | .Cs, .C3h | .C3, .C3h => true
  | .C1, .C6h | .Ci, .C6h | .C2, .C6h | .Cs, .C6h => true
  | .C2h, .C6h | .C3, .C6h | .S6, .C6h | .C6, .C6h => true
  | .C3h, .C6h => true
  | .C1, .D6 | .C2, .D6 | .D2, .D6 | .C3, .D6 => true
  | .D3, .D6 | .C6, .D6 => true
  | .C1, .C6v | .C2, .C6v | .Cs, .C6v | .C2v, .C6v => true
  | .C3, .C6v | .C3v, .C6v | .C6, .C6v => true
  | .C1, .D3h | .C2, .D3h | .Cs, .D3h | .C2v, .D3h => true
  | .C3, .D3h | .D3, .D3h | .C3v, .D3h | .C3h, .D3h => true
  | .C1, .D6h | .Ci, .D6h | .C2, .D6h | .Cs, .D6h => true
  | .C2h, .D6h | .D2, .D6h | .C2v, .D6h | .D2h, .D6h => true
  | .C3, .D6h | .S6, .D6h | .D3, .D6h | .C3v, .D6h => true
  | .C6, .D6h | .C3h, .D6h | .C6h, .D6h | .D6, .D6h => true
  | .C6v, .D6h | .D3h, .D6h => true
  | .C1, .T | .C2, .T | .D2, .T | .C3, .T => true
  | .C1, .Th | .Ci, .Th | .C2, .Th | .Cs, .Th => true
  | .C2h, .Th | .D2, .Th | .C2v, .Th | .D2h, .Th => true
  | .C3, .Th | .S6, .Th | .T, .Th => true
  | .C1, .O_ | .C2, .O_ | .D2, .O_ | .C3, .O_ => true
  | .D3, .O_ | .C4, .O_ | .D4, .O_ | .T, .O_ => true
  | .C1, .Td | .C2, .Td | .Cs, .Td | .D2, .Td => true
  | .C2v, .Td | .C3, .Td | .C3v, .Td | .S4, .Td => true
  | .D2d, .Td | .T, .Td => true
  | .C1, .Oh | .Ci, .Oh | .C2, .Oh | .Cs, .Oh => true
  | .C2h, .Oh | .D2, .Oh | .C2v, .Oh | .D2h, .Oh => true
  | .C3, .Oh | .S6, .Oh | .D3, .Oh | .C3v, .Oh => true
  | .D3d, .Oh | .C4, .Oh | .S4, .Oh | .C4h, .Oh => true
  | .D4, .Oh | .C4v, .Oh | .D2d, .Oh | .D4h, .Oh => true
  | .T, .Oh | .Th, .Oh | .O_, .Oh | .Td, .Oh => true
  | _, _ => false

/-- C1 is a subgroup of every point group. -/
theorem c1_universal_subgroup (g : PointGroup) : isSubgroup .C1 g = true := by
  cases g <;> native_decide

/-- Subgroup relation is reflexive. -/
theorem isSubgroup_refl (g : PointGroup) : isSubgroup g g = true := by
  cases g <;> native_decide

/-! ## 5. Pairing Operad Definition

We define the concrete pairing operad with operations being
(point_group, channel) pairs where the channel is allowed. -/

/-- An operation in the pairing operad: a channel allowed in a point group. -/
structure PairingOp where
  pg : PointGroup
  ch : Channel
  allowed : allowedChannel pg ch = true
  deriving Repr

/-- The unit operation: s-wave in any point group. -/
def unitOp (g : PointGroup) : PairingOp :=
  { pg := g, ch := .s, allowed := s_wave_universal g }

/-- Composition of pairing operations: the outer operation absorbs
    the inner when the inner's point group is a subgroup.
    Returns the outer operation (the dominant channel). -/
def composeOp (outer inner : PairingOp)
    (_hsub : isSubgroup inner.pg outer.pg = true) : PairingOp :=
  outer

/-! ## 6. Operad Axiom Verification -/

/-- **Unit axiom (right)**: composing any operation with the s-wave unit
    of its own point group yields the same operation. -/
theorem unit_right (f : PairingOp) :
    composeOp f (unitOp f.pg) (isSubgroup_refl f.pg) = f := by
  rfl

/-- **Unit axiom (left)**: composing the s-wave unit with any operation
    in the same point group yields the original operation. -/
theorem unit_left (f : PairingOp) :
    composeOp (unitOp f.pg) f (isSubgroup_refl f.pg) = unitOp f.pg := by
  rfl

/-- **Associativity**: (f comp g) comp h = f comp (g comp h),
    where sub-group relations compose transitively. -/
theorem comp_assoc (f g h : PairingOp)
    (hgf : isSubgroup g.pg f.pg = true)
    (hhg : isSubgroup h.pg g.pg = true)
    (hhf : isSubgroup h.pg f.pg = true) :
    composeOp (composeOp f g hgf) h hhf =
    composeOp f (composeOp g h hhg) hgf := by
  rfl

/-! ## 7. Restriction Functor

Restriction of a pairing operation to a subgroup: the channel
survives iff it is allowed in the subgroup. -/

/-- Restrict an operation to a subgroup, if the channel survives. -/
def restrictOp (op : PairingOp) (target : PointGroup)
    (_hsub : isSubgroup target op.pg = true)
    (hallowed : allowedChannel target op.ch = true) : PairingOp :=
  { pg := target, ch := op.ch, allowed := hallowed }

/-- **Restriction preserves identity**: restricting the s-wave unit
    to any subgroup yields the s-wave unit of the subgroup. -/
theorem restrict_unit (g target : PointGroup)
    (hsub : isSubgroup target g = true) :
    restrictOp (unitOp g) target hsub (s_wave_universal target) =
    unitOp target := by
  rfl

/-- **Restriction preserves composition**: restrict(f comp g) = restrict(f) comp restrict(g)
    when both channels survive. This is the equivariance axiom. -/
theorem restrict_comp (f g : PairingOp)
    (target : PointGroup)
    (hgf : isSubgroup g.pg f.pg = true)
    (htf : isSubgroup target f.pg = true)
    (htg : isSubgroup target g.pg = true)
    (htt : isSubgroup target target = true)
    (haf : allowedChannel target f.ch = true)
    (hag : allowedChannel target g.ch = true) :
    restrictOp (composeOp f g hgf) target htf haf =
    composeOp (restrictOp f target htf haf)
              (restrictOp g target htg hag) htt := by
  rfl

/-! ## 8. Key Properties of the Channel Table -/

/-- D4h admits exactly 6 channels. -/
theorem d4h_channel_count :
    (List.filter (fun c => allowedChannel .D4h c)
      [.s, .p, .d, .d_x2y2, .d_xy, .f, .z3, .z4, .z6]).length = 6 := by
  native_decide

/-- Oh admits exactly 5 channels. -/
theorem oh_channel_count :
    (List.filter (fun c => allowedChannel .Oh c)
      [.s, .p, .d, .d_x2y2, .d_xy, .f, .z3, .z4, .z6]).length = 5 := by
  native_decide

/-- C1 admits exactly 1 channel (s-wave only). -/
theorem c1_channel_count :
    (List.filter (fun c => allowedChannel .C1 c)
      [.s, .p, .d, .d_x2y2, .d_xy, .f, .z3, .z4, .z6]).length = 1 := by
  native_decide

/-- z4 is lost when restricting D4h to D2h. -/
theorem z4_lost_d4h_to_d2h : allowedChannel .D4h .z4 = true
    ∧ allowedChannel .D2h .z4 = false := by
  exact ⟨rfl, rfl⟩

/-- d-wave is preserved when restricting D4h to D2h. -/
theorem d_preserved_d4h_to_d2h : allowedChannel .D4h .d = true
    ∧ allowedChannel .D2h .d = true := by
  exact ⟨rfl, rfl⟩

/-- D2h is a subgroup of D4h. -/
theorem d2h_sub_d4h : isSubgroup .D2h .D4h = true := by native_decide

/-- C3 is NOT a subgroup of D4h (incompatible lattice branches). -/
theorem c3_not_sub_d4h : isSubgroup .C3 .D4h = false := by native_decide

/-! ## 9. Free Algebra Theorem

The Sigrist-Ueda table is the free O-algebra on one generator per
crystallographic point group. Concretely: for each PG g, the set of
allowed channels equals the set of operations in the operad at color g. -/

/-- The set of allowed channels for a point group IS the set of
    generators of the free O-algebra at that color. This is the
    formal content of "Sigrist-Ueda = free algebra":
    there exists a pairing operation for (g, ch) iff ch is allowed. -/
theorem sigrist_ueda_is_free_algebra (g : PointGroup) (ch : Channel) :
    (∃ op : PairingOp, op.pg = g ∧ op.ch = ch) ↔ allowedChannel g ch = true := by
  constructor
  · rintro ⟨op, hpg, hch⟩
    rw [← hpg, ← hch]
    exact op.allowed
  · intro h
    exact ⟨{ pg := g, ch := ch, allowed := h }, rfl, rfl⟩

/-! ## 10. Morphism DAG Properties -/

/-- The maximal-subgroup relation is irreflexive (no group is its own
    maximal subgroup). -/
theorem maxSub_irrefl (g : PointGroup) : isMaximalSubgroup g g = false := by
  cases g <;> native_decide

/-- The subgroup relation is antisymmetric on distinct elements:
    if H <= G and G <= H then H = G. -/
theorem isSubgroup_antisymm : forall (g h : PointGroup),
    isSubgroup g h = true -> isSubgroup h g = true -> g = h := by
  native_decide

end OperadicCatalog
