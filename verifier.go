package main

import (
	"fmt"
	"math"
	"sort"
)

const madScaleFactor = 2.5

// VerificationResult holds the outcome of subgraph isomorphism verification.
type VerificationResult struct {
	Verified       bool     // true if unchanged region is isomorphic
	UnchangedAtoms int      // number of atoms below the MAD threshold
	TotalAtoms     int      // total mapped atoms
	FailedPairs    [][2]int // reactant atom pairs whose edge is missing in product
	Threshold      float64  // adaptive MAD threshold used
	ReactionCentre []int    // 1-based reactant atom indices classified as changed
}

// Verifier checks that a Mapping is consistent with the molecular graphs.
type Verifier interface {
	Verify(rxn *Reaction, m Mapping, rFPs, pFPs [][][]float64) VerificationResult
}

// SubgraphVerifier implements Verifier using MAD-based adaptive thresholding
// and bond-neighbourhood isomorphism checking.
type SubgraphVerifier struct{}

// Verify implements Verifier.
func (v *SubgraphVerifier) Verify(
	rxn *Reaction,
	m Mapping,
	rFPs, pFPs [][][]float64,
) VerificationResult {
	pairs := v.collectPairs(m, rFPs, pFPs)
	if len(pairs) == 0 {
		return VerificationResult{}
	}

	threshold := v.madThreshold(pairs)
	unchanged, reactionCentre, rToP := v.classify(pairs, threshold)
	failedPairs, additionalRC := v.checkEdges(rxn, unchanged, rToP)
	reactionCentre = dedupSorted(append(reactionCentre, additionalRC...))

	return VerificationResult{
		Verified:       len(failedPairs) == 0,
		UnchangedAtoms: len(unchanged),
		TotalAtoms:     len(pairs),
		FailedPairs:    failedPairs,
		Threshold:      threshold,
		ReactionCentre: reactionCentre,
	}
}

type assignedPair struct {
	rMol, rAtom int
	pMol, pAtom int
	dist        float64
}

func (v *SubgraphVerifier) collectPairs(m Mapping, rFPs, pFPs [][][]float64) []assignedPair {
	mapNumToEnc := make(map[int]int)
	for enc, mn := range m.PMap {
		mapNumToEnc[mn] = enc
	}

	var pairs []assignedPair
	for rKey, mn := range m.RMap {
		if mn == 0 {
			continue
		}
		enc, ok := mapNumToEnc[mn]
		if !ok {
			continue
		}
		pMol, pAtom := decodePA(enc)
		rFP := rFPs[rKey[0]][rKey[1]]
		pFP := make([]float64, maxLevel)
		if pMol < len(pFPs) && pAtom < len(pFPs[pMol]) {
			pFP = pFPs[pMol][pAtom]
		}
		pairs = append(pairs, assignedPair{
			rMol: rKey[0], rAtom: rKey[1],
			pMol: pMol, pAtom: pAtom,
			dist: vecDist(rFP, pFP),
		})
	}
	return pairs
}

func (v *SubgraphVerifier) madThreshold(pairs []assignedPair) float64 {
	dists := make([]float64, len(pairs))
	for i, p := range pairs {
		dists[i] = p.dist
	}

	sorted := make([]float64, len(dists))
	copy(sorted, dists)
	sort.Float64s(sorted)
	median := medianSorted(sorted)

	absDevs := make([]float64, len(dists))
	for i, d := range dists {
		absDevs[i] = math.Abs(d - median)
	}
	sort.Float64s(absDevs)
	mad := medianSorted(absDevs)

	if mad < 1e-9 {
		mad = 1e-9
	}
	return median + madScaleFactor*mad
}

func medianSorted(sorted []float64) float64 {
	n := len(sorted)
	if n == 0 {
		return 0
	}
	if n%2 == 0 {
		return (sorted[n/2-1] + sorted[n/2]) / 2.0
	}
	return sorted[n/2]
}

func (v *SubgraphVerifier) classify(
	pairs []assignedPair,
	threshold float64,
) (unchanged map[[2]int]bool, reactionCentre []int, rToP map[[2]int][2]int) {
	unchanged = make(map[[2]int]bool)
	rToP = make(map[[2]int][2]int)

	for _, p := range pairs {
		rKey := [2]int{p.rMol, p.rAtom}
		rToP[rKey] = [2]int{p.pMol, p.pAtom}
		if p.dist <= threshold {
			unchanged[rKey] = true
		} else {
			reactionCentre = append(reactionCentre, p.rAtom+1)
		}
	}
	return unchanged, reactionCentre, rToP
}

func (v *SubgraphVerifier) checkEdges(
	rxn *Reaction,
	unchanged map[[2]int]bool,
	rToP map[[2]int][2]int,
) (failedPairs [][2]int, additionalRC []int) {
	seen := map[[2]int]bool{}

	for rKey := range unchanged {
		rMol := rxn.Reactants[rKey[0]]
		rAtom := rKey[1]

		for _, nb := range rMol.Adjacency[rAtom] {
			nbKey := [2]int{rKey[0], nb}
			if !unchanged[nbKey] {
				continue
			}

			// Canonical edge key to avoid duplicates
			a, b := rAtom+1, nb+1
			if a > b {
				a, b = b, a
			}
			edgeKey := [2]int{a, b}
			if seen[edgeKey] {
				continue
			}
			seen[edgeKey] = true

			rBO := rMol.BondOrder[[2]int{rAtom, nb}]
			if rBO == 0 {
				rBO = 1
			}

			pKey := rToP[rKey]
			pNbKey := rToP[nbKey]

			if pKey[0] != pNbKey[0] {
				// Cross-molecule bond: leaving group, not a mapping error
				additionalRC = append(additionalRC, rAtom+1, nb+1)
				continue
			}

			pMol := rxn.Products[pKey[0]]
			pBO := v.lookupBondOrder(pMol, pKey[1], pNbKey[1])

			if pBO == 0 {
				failedPairs = append(failedPairs, edgeKey)
			} else if rBO != pBO {
				// Bond order changed at this edge — reaction centre
				additionalRC = append(additionalRC, rAtom+1)
			}
		}
	}
	return failedPairs, additionalRC
}

func (v *SubgraphVerifier) lookupBondOrder(mol *Molecule, a1, a2 int) int {
	if bo, ok := mol.BondOrder[[2]int{a1, a2}]; ok && bo > 0 {
		return bo
	}
	// Fallback: check adjacency list
	for _, nb := range mol.Adjacency[a1] {
		if nb == a2 {
			bo := mol.BondOrder[[2]int{a1, nb}]
			if bo == 0 {
				return 1
			}
			return bo
		}
	}
	return 0
}

func dedupSorted(in []int) []int {
	seen := map[int]bool{}
	var out []int
	for _, v := range in {
		if !seen[v] {
			seen[v] = true
			out = append(out, v)
		}
	}
	sort.Ints(out)
	return out
}

// PrintVerification prints a human-readable summary to stdout.
func PrintVerification(res VerificationResult) {
	fmt.Printf("\n=== Subgraph Isomorphism Verification ===\n")
	fmt.Printf("Total assigned atoms : %d\n", res.TotalAtoms)
	fmt.Printf("Unchanged region     : %d atoms (dist <= median+%.1f*MAD, threshold=%.4f)\n",
		res.UnchangedAtoms, madScaleFactor, res.Threshold)
	fmt.Printf("Reaction centre atoms: %v\n", res.ReactionCentre)

	if res.Verified {
		fmt.Println("Result               : VERIFIED — unchanged subgraph is isomorphic")
	} else {
		fmt.Printf("Result               : FAILED — %d inconsistent edge(s)\n", len(res.FailedPairs))
		for _, fp := range res.FailedPairs {
			fmt.Printf("  Edge missing in product: R-atom %d -- R-atom %d\n", fp[0], fp[1])
		}
		fmt.Println("  → Mapping may be incorrect; consider manual review")
	}
}
