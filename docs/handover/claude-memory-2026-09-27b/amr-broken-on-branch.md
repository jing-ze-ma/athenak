---
name: amr-broken-on-branch
description: AMR segfaulted on the cubed-sphere branch because Mesh::ptree was never constructed; fixed in 1e432b81
metadata: 
  node_type: memory
  type: project
  originSessionId: 75f1abe2-02a4-4f90-9678-ee140865bae8
  modified: 2026-08-08T15:05:56.652Z
---

**FIXED in `1e432b81`** (2026-08-08). Kept because the shape of the bug says something
about the cubed-sphere refactor that is still half-finished elsewhere.

`9e37cfb7` ("Begin constructing cubed-sphere: 6 panels") replaced the single
`MeshBlockTree` with one tree per panel and **commented out the construction of
`Mesh::ptree`** in `BuildTreeFromScratch`, without updating the code that still reads it.
`ptree` was a null `unique_ptr` for every from-scratch run, and
`MeshRefinement::UpdateMeshBlockTree` dereferenced it on the first refinement check.

Debugging note worth reusing: the fault address **0x18 is the byte offset of `lloc_.level`**
in `MeshBlockTree` (`pleaf_`@0, `gid_`@8, `lloc_`@12 so `level`@24). A segfault at 0x18 in
that class means a null node, not a corrupt tree. Guarding the `FindMeshBlock` return did
NOT fire, which is what proved the null was `this` rather than the result.

The migration was half-done in the other direction too: `BuildTreeFromRestart` built only
`ptree` and never populated `panel_trees`, which `MeshBlock`'s neighbour search indexes
directly — a restart read past the end of an empty vector.

Fix: `panel_trees` owns the nodes in both paths; `ptree` is a non-owning pointer to panel 0.
AMR with `npanels > 1` is now refused (it would refine one sixth of the mesh).

**Suspect the same half-migration elsewhere.** `shearing_box.cpp:120` uses `pm->ptree->`
and was equally broken before this fix; anything else added since that assumes one global
tree is worth checking against the panel decomposition.

See [[general-eos-project]].
