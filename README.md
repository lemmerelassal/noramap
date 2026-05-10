# NoraMap

**Normalised-rank Optimal Reaction Atom Mapper**

NoraMap is a command-line tool for atom-to-atom mapping of chemical reactions. It reads and writes V2000 RXN files and requires no training data or internet connection beyond the optional KEGG download mode.

## How it works

Each atom in a molecule receives a vector fingerprint computed by breadth-first traversal to depth 6. Each component of the vector is the product of prime numbers assigned to the atoms at that shell distance, weighted by bond order, ring membership, formal charge, and stereo configuration. After the initial BFS pass, the fingerprints are refined iteratively in the style of the Morgan algorithm until the relative ranking stabilises. The raw floating-point values are then replaced by their per-component ranks across all atoms in the molecule, making the fingerprint reproducible across platforms and independent of atom ordering in the input file.

Mapping is solved per element group using the Hungarian algorithm, which gives a globally optimal bijection under the rank-distance cost metric. Atoms that cannot be matched in the primary assignment (for example, hydrogens that leave as H2) are assigned in a second pass by nearest-neighbour search. After assignment, the unchanged region of the molecule is verified by checking that every bond between two mapped atoms is preserved in the product with the same bond order. The adaptive threshold for classifying atoms as changed or unchanged is based on the median absolute deviation of the assignment distances, which makes it robust to reactions where a single atom changes dramatically.

Explicit hydrogens are added via OpenBabel before fingerprinting, so KEGG mol files (which store hydrogens implicitly) are handled correctly.

## Installation

Requires Go 1.21 or later and OpenBabel.

```
git clone https://github.com/placeholder/noramap
cd noramap
go build -o noramap .
```

OpenBabel must be on your PATH. On Windows, download the installer from openbabel.org and add the installation directory to PATH. On Linux:

```
apt install openbabel
```

## Usage

Map a single reaction:

```
noramap input.rxn output.rxn
```

Download and map all KEGG reactions:

```
noramap kegg
```

Resume an interrupted KEGG run:

```
noramap kegg resume:R01234
```

Input and output files must be V2000 RXN format. The output file uses CRLF line endings and is compatible with ACD/ChemSketch and similar software.

KEGG reactions are written to `input/` and mapped outputs to `output/`, both relative to the working directory. Reactions already present in both directories are skipped automatically on resume.

## Output

The console output shows a mapping table and a verification result for each reaction. For example:

```
Reactants: 1  Products: 2

RMol RAtom  Elem   PMol PAtom  Elem   R.FP       P.FP       Dist      MapNum
--------------------------------------------------------------------------------
R1   A1     C      P1   A1     C      11.58      8.66       4.36      1
R1   A2     C      P1   A2     C      11.75      8.78       3.00      2
R1   A8     H      P2   A1     H      9.00       0.00       9.00      7 ~
...

=== Subgraph Isomorphism Verification ===
Total assigned atoms : 9
Unchanged region     : 7 atoms (dist <= median+2.5*MAD, threshold=8.27)
Reaction centre atoms: [2 8 9]
Result               : VERIFIED - unchanged subgraph is isomorphic
```

The `~` marker indicates a second-pass (fallback) assignment where the atom's nearest neighbour was already claimed. The reaction centre atoms listed are those whose fingerprint distance exceeds the MAD threshold and therefore changed environment between reactant and product.

## Supported elements

All common organic elements plus biologically relevant metals: Fe, Mg, Zn, Ca, Mn, Co, Cu, Ni, Mo, V, W, Cr, and alkali/alkaline earth metals. Unknown elements fall back to a default prime value.

## Limitations

- Stereochemistry encoding is present but sensitivity is lower in small molecules where few atoms have distinguishable rank positions.
- Large symmetric systems (porphyrins, nucleotide repeats) may produce tied fingerprints for genuinely equivalent atoms. The assignment within a tied group is arbitrary but chemically valid.
- The accuracy has not been benchmarked against a standard dataset such as USPTO-50K. Estimates in the accompanying paper are based on algorithmic analysis.
- OpenBabel must be installed separately.

## Code structure

The codebase follows SOLID principles with one file per responsibility:

| File | Responsibility |
|---|---|
| `molecule.go` | Core data types (Atom, Molecule, Reaction, Mapping) |
| `elements.go` | Element-to-prime assignment table |
| `molparser.go` | V2000 mol block parsing |
| `rxnparser.go` | RXN file parsing, dependency injection |
| `hydrogen.go` | Explicit hydrogen addition via OpenBabel |
| `fingerprint.go` | Vector fingerprint computation and rank canonicalisation |
| `mapper.go` | Hungarian algorithm assignment grouped by element |
| `hungarian.go` | Hungarian algorithm implementation (Kuhn-Munkres) |
| `verifier.go` | Subgraph isomorphism verification with MAD threshold |
| `rxnwriter.go` | V2000 RXN file writing with CRLF line endings |
| `pipeline.go` | Single-reaction orchestration with injected dependencies |
| `kegg.go` | KEGG download pipeline with parallel worker pools |
| `main.go` | CLI entry point |

## Reference

El Assal Lemmer Reiad Petro. "A Prime-Product Vector Fingerprinting Algorithm with Iterative Morgan Refinement, Stereochemistry and Ionic State Encoding, Rank Canonicalisation, Optimal Hungarian Assignment, Adaptive Reaction Centre Detection, and Subgraph Isomorphism Verification for Atom-to-Atom Mapping in Biochemical Reactions." Preprint, 2026.

## License

MIT
