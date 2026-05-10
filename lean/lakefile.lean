import Lake
open Lake DSL

package «operadic» where
  leanOptions := #[
    { name := `autoImplicit, value := .ofBool false }
  ]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "master"

@[default_target]
lean_lib «OperadicCatalog» where
  roots := #[`OperadicCatalog]
