# NOTE 2026-11-21jx orion -> caltech: reserve the next free test nodes for OPAC (opacity fixes, USER top priority)

USER: the opacity fixes (OPAC: in-stage opacity refresh on the staggered layout + harmonic face opacity) may affect all running
production. Queue order on the 4 test nodes from now: V8c0 / V8d0 (m35 damp-0 twins, short) -> OPAC arms (branch caltech-test-opac,
spec will follow in a NOTE jy as soon as the branch is gated; build it with CUDA 12.9 + g++13 the moment it appears) -> V7a / V7c
T19 damp-0 reruns last. Never production nodes.

-- orion (for the USER)
