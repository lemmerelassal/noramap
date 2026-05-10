package main

import (
	"math"
	"sort"
)

const (
	maxLevel    = 6
	morganIters = 10

	// Stereo multipliers: incommensurate with primes to prevent cancellations.
	stereoNone     = 1.0
	stereoR        = 1.03
	stereoS        = 1.07
	stereoEither   = 1.01
	bondStereoUp   = 1.013
	bondStereoDown = 1.019
)

// Fingerprinter computes vector fingerprints for all atoms in a molecule.
// The interface allows alternative fingerprint strategies to be substituted
// without modifying the mapper (O/D principle).
type Fingerprinter interface {
	Compute(mol *Molecule) [][]float64
}

// MorganFingerprinter implements Fingerprinter using iterative Morgan-style
// refinement of prime-product BFS vectors, followed by rank canonicalisation.
type MorganFingerprinter struct{}

// Compute implements Fingerprinter.
func (mf *MorganFingerprinter) Compute(mol *Molecule) [][]float64 {
	n := len(mol.Atoms)
	fps := make([][]float64, n)
	for i := range mol.Atoms {
		fps[i] = mf.initialVector(mol, i)
	}
	fps = mf.refine(mol, fps)
	return rankComponents(fps)
}

// initialVector computes the base BFS prime-product vector for one root atom.
func (mf *MorganFingerprinter) initialVector(mol *Molecule, start int) []float64 {
	type entry struct{ idx, level, parent int }
	visited := make([]bool, len(mol.Atoms))
	queue := []entry{{start, 1, -1}}
	visited[start] = true
	byLevel := map[int][]entry{}

	for len(queue) > 0 {
		e := queue[0]
		queue = queue[1:]
		byLevel[e.level] = append(byLevel[e.level], e)
		if e.level < maxLevel {
			for _, nb := range mol.Adjacency[e.idx] {
				if !visited[nb] {
					visited[nb] = true
					queue = append(queue, entry{nb, e.level + 1, e.idx})
				}
			}
		}
	}

	vec := make([]float64, maxLevel)
	for lv := 1; lv <= maxLevel; lv++ {
		entries := byLevel[lv]
		if len(entries) == 0 {
			break
		}
		prod := 1.0
		for _, e := range entries {
			prod *= atomWeight(mol, e.idx, e.parent)
		}
		vec[lv-1] = prod
	}
	return vec
}

// refine applies iterative Morgan-style propagation until convergence or maxIter.
func (mf *MorganFingerprinter) refine(mol *Molecule, fps [][]float64) [][]float64 {
	n := len(fps)
	for iter := 0; iter < morganIters; iter++ {
		prev := make([][]float64, n)
		for i := range fps {
			prev[i] = make([]float64, maxLevel)
			copy(prev[i], fps[i])
		}

		for i := range mol.Atoms {
			next := make([]float64, maxLevel)
			copy(next, prev[i])
			for _, nb := range mol.Adjacency[i] {
				bo := float64(mol.BondOrder[[2]int{i, nb}])
				if bo == 0 {
					bo = 1
				}
				if bo == 4 {
					bo = 1.5 // aromatic
				}
				for k := range next {
					next[k] += prev[nb][k] * bo / 10.0
				}
			}
			fps[i] = next
		}

		if rankingUnchanged(prev, fps) {
			break
		}
	}
	return fps
}

// rankingUnchanged returns true if the lexicographic ordering of all vectors
// is identical between prev and curr — used as a convergence criterion.
func rankingUnchanged(prev, curr [][]float64) bool {
	for i := 0; i < len(prev)-1; i++ {
		if vecLess(prev[i], prev[i+1]) != vecLess(curr[i], curr[i+1]) {
			return false
		}
	}
	return true
}

// rankComponents replaces each component of every vector with its 0-based rank
// among all atoms for that component. Ties receive the same (min) rank.
// This eliminates floating-point magnitude sensitivity while preserving all
// relative structural distinctions.
func rankComponents(fps [][]float64) [][]float64 {
	n := len(fps)
	if n == 0 {
		return fps
	}
	ranked := make([][]float64, n)
	for i := range ranked {
		ranked[i] = make([]float64, maxLevel)
	}

	type pair struct {
		val float64
		idx int
	}

	for k := 0; k < maxLevel; k++ {
		pairs := make([]pair, n)
		for i := range fps {
			pairs[i] = pair{fps[i][k], i}
		}
		sort.Slice(pairs, func(a, b int) bool {
			return pairs[a].val < pairs[b].val
		})
		rank := 0.0
		for i := 0; i < len(pairs); {
			j := i
			for j < len(pairs) && pairs[j].val == pairs[i].val {
				j++
			}
			for _, p := range pairs[i:j] {
				ranked[p.idx][k] = rank
			}
			rank += float64(j - i)
			i = j
		}
	}
	return ranked
}

// atomWeight computes the composite prime-product weight for atom atomIdx
// with BFS parent parentIdx, encoding element, bond order, ring membership,
// formal charge, and stereochemistry.
func atomWeight(mol *Molecule, atomIdx, parentIdx int) float64 {
	prime := ElementPrime(mol.Atoms[atomIdx].Symbol)
	boMult := bondOrderMultiplier(mol, atomIdx, parentIdx)
	ringMult := ringMultiplier(mol, atomIdx)
	chargeMult := chargeMultiplier(mol.Atoms[atomIdx].FormalCharge)
	stereoMult := stereoMultiplier(mol.Atoms[atomIdx].StereoParity)
	bondStereoMult := bondStereoMultiplier(mol, atomIdx, parentIdx)
	return prime * boMult * ringMult * chargeMult * stereoMult * bondStereoMult
}

func bondOrderMultiplier(mol *Molecule, atomIdx, parentIdx int) float64 {
	if parentIdx < 0 {
		return 1.0
	}
	switch mol.BondOrder[[2]int{parentIdx, atomIdx}] {
	case 2:
		return 2.0
	case 3:
		return 3.0
	case 4:
		return 1.5 // aromatic
	default:
		return 1.0
	}
}

func ringMultiplier(mol *Molecule, atomIdx int) float64 {
	if isInRing(mol, atomIdx) {
		return 1.1
	}
	return 1.0
}

// isInRing returns true if atomIdx is part of any ring.
// Uses DFS: block atomIdx, check if one neighbour can reach another.
func isInRing(mol *Molecule, atomIdx int) bool {
	neighbors := mol.Adjacency[atomIdx]
	if len(neighbors) < 2 {
		return false
	}
	start, target := neighbors[0], neighbors[1]
	visited := make([]bool, len(mol.Atoms))
	visited[atomIdx] = true
	visited[start] = true
	stack := []int{start}
	for len(stack) > 0 {
		cur := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		if cur == target {
			return true
		}
		for _, nb := range mol.Adjacency[cur] {
			if !visited[nb] {
				visited[nb] = true
				stack = append(stack, nb)
			}
		}
	}
	return false
}

func chargeMultiplier(charge int) float64 {
	switch charge {
	case 3:
		return 1.093
	case 2:
		return 1.061
	case 1:
		return 1.031
	case -1:
		return 0.971
	case -2:
		return 0.943
	case -3:
		return 0.917
	default:
		return 1.0
	}
}

func stereoMultiplier(parity int) float64 {
	switch parity {
	case 1:
		return stereoR
	case 2:
		return stereoS
	case 3:
		return stereoEither
	default:
		return stereoNone
	}
}

func bondStereoMultiplier(mol *Molecule, atomIdx, parentIdx int) float64 {
	if parentIdx < 0 {
		return 1.0
	}
	switch mol.BondStereo[[2]int{parentIdx, atomIdx}] {
	case 1:
		return bondStereoUp
	case 6:
		return bondStereoDown
	default:
		return 1.0
	}
}

// vecDist returns the Euclidean distance between two fingerprint vectors.
func vecDist(a, b []float64) float64 {
	sum := 0.0
	for i := range a {
		d := a[i] - b[i]
		sum += d * d
	}
	return math.Sqrt(sum)
}

// vecLess returns true if a is lexicographically less than b.
func vecLess(a, b []float64) bool {
	for i := range a {
		if a[i] != b[i] {
			return a[i] < b[i]
		}
	}
	return false
}

// vecNorm returns the L2 norm of a vector.
func vecNorm(v []float64) float64 {
	s := 0.0
	for _, x := range v {
		s += x * x
	}
	return math.Sqrt(s)
}
