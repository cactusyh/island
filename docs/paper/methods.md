# Methods

The methods follow the source-bound workflow implemented in ISLAND. A declared
PSMILES and construction policy produce a final chemical graph with stable atom
IDs. Sequence and stereochemical choices are serialized before deterministic
crosslink plans are applied. Periodic packing creates an immutable remapping and
box record. Force-field assignment then binds a unified TypedGraph, charge
record and parameter source identity for a bounded PCFF, OPLS-AA, GAFF or GAFF2
profile.

J23 constructs periodic OpenMM Reference systems for minimization, NVT, NPT and
NVT annealing. Accepted restart boundaries publish binary checkpoints and signed
generation records. The opt-in trajectory store writes bounded JSONL chunks and
durable restart checkpoints independently. Analysis requires explicit stable-ID
mapping, periodic boxes, masses and source hashes; each observable carries units,
requirements, uncertainty method and a blocked state.
