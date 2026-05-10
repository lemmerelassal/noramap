package main

import (
	"math"
	"sort"
)

// AtomMapper is the interface for computing a Mapping from a Reaction.
// Separating the interface allows alternative mapping strategies (greedy,
// ILP, ML-based) to be swapped in without modifying callers (O/D).
type AtomMapper interface {
	Map(rxn *Reaction) Mapping
}

// FlatAtom is an atom with its pre-computed fingerprint vector, used
// internally by HungarianMapper.
type FlatAtom struct {
	MolIdx  int
	AtomIdx int
	FP      []float64
	Symbol  string
}

// HungarianMapper implements AtomMapper using the Hungarian algorithm grouped
// by element, with a second-pass fallback for unmatched atoms.
type HungarianMapper struct {
	Fingerprinter Fingerprinter
}

// NewHungarianMapper constructs a HungarianMapper with the standard
// MorganFingerprinter injected.
func NewHungarianMapper() *HungarianMapper {
	return &HungarianMapper{Fingerprinter: &MorganFingerprinter{}}
}

// Map implements AtomMapper.
func (m *HungarianMapper) Map(rxn *Reaction) Mapping {
	rFlat := m.flatten(rxn.Reactants)
	pFlat := m.flatten(rxn.Products)

	rMap := make(map[[2]int]int)
	pMap := make(map[int]int)
	usedP := make(map[int]bool)
	mapNum := 1

	for _, elem := range sortedElements(rFlat, pFlat) {
		rGroup := filterByElement(rFlat, elem)
		pGroup := filterByElement(pFlat, elem)
		if len(rGroup) == 0 || len(pGroup) == 0 {
			continue
		}

		cost := buildCostMatrix(rGroup, pGroup)
		assignment := hungarian(cost)

		for i, j := range assignment {
			if j < 0 || j >= len(pGroup) {
				continue
			}
			enc := encodePA(pGroup[j].MolIdx, pGroup[j].AtomIdx)
			rKey := [2]int{rGroup[i].MolIdx, rGroup[i].AtomIdx}
			rMap[rKey] = mapNum
			pMap[enc] = mapNum
			usedP[enc] = true
			mapNum++
		}
	}

	// Second pass: fallback nearest-neighbour for unmatched reactant atoms
	for _, ra := range rFlat {
		rKey := [2]int{ra.MolIdx, ra.AtomIdx}
		if rMap[rKey] != 0 {
			continue
		}
		enc, found := nearestUnused(ra, pFlat, usedP)
		if found {
			usedP[enc] = true
			rMap[rKey] = mapNum
			pMap[enc] = mapNum
			mapNum++
		}
	}

	return Mapping{RMap: rMap, PMap: pMap}
}

func (m *HungarianMapper) flatten(mols []*Molecule) []FlatAtom {
	var flat []FlatAtom
	for mi, mol := range mols {
		fps := m.Fingerprinter.Compute(mol)
		for ai, fp := range fps {
			flat = append(flat, FlatAtom{
				MolIdx:  mi,
				AtomIdx: ai,
				FP:      fp,
				Symbol:  mol.Atoms[ai].Symbol,
			})
		}
	}
	return flat
}

// FlatFingerprints returns pre-computed fingerprints as [molIdx][atomIdx][]float64,
// used by the verifier to avoid recomputing.
func (m *HungarianMapper) FlatFingerprints(mols []*Molecule) [][][]float64 {
	result := make([][][]float64, len(mols))
	for mi, mol := range mols {
		result[mi] = m.Fingerprinter.Compute(mol)
	}
	return result
}

func sortedElements(rFlat, pFlat []FlatAtom) []string {
	seen := map[string]bool{}
	for _, a := range rFlat {
		seen[a.Symbol] = true
	}
	for _, a := range pFlat {
		seen[a.Symbol] = true
	}
	elems := make([]string, 0, len(seen))
	for e := range seen {
		elems = append(elems, e)
	}
	sort.Strings(elems)
	return elems
}

func filterByElement(flat []FlatAtom, elem string) []FlatAtom {
	var result []FlatAtom
	for _, a := range flat {
		if a.Symbol == elem {
			result = append(result, a)
		}
	}
	return result
}

func buildCostMatrix(rGroup, pGroup []FlatAtom) [][]float64 {
	cost := make([][]float64, len(rGroup))
	for i, ra := range rGroup {
		cost[i] = make([]float64, len(pGroup))
		for j, pa := range pGroup {
			cost[i][j] = vecDist(ra.FP, pa.FP)
		}
	}
	return cost
}

func nearestUnused(ra FlatAtom, pFlat []FlatAtom, usedP map[int]bool) (int, bool) {
	bestDiff := math.MaxFloat64
	bestEnc := -1
	for _, pa := range pFlat {
		if pa.Symbol != ra.Symbol {
			continue
		}
		enc := encodePA(pa.MolIdx, pa.AtomIdx)
		if usedP[enc] {
			continue
		}
		d := vecDist(ra.FP, pa.FP)
		if d < bestDiff {
			bestDiff = d
			bestEnc = enc
		}
	}
	return bestEnc, bestEnc >= 0
}
