package main

// Atom represents a single atom in a molecule.
type Atom struct {
	Index        int
	Symbol       string
	X, Y, Z      float64
	StereoParity int // V2000 parity: 0=none, 1=odd(R), 2=even(S), 3=either
	FormalCharge int // decoded from V2000 charge code or M  CHG line
}

// Molecule holds all structural data for a single molecule.
type Molecule struct {
	Atoms      []Atom
	Adjacency  map[int][]int  // atom index -> neighbour indices
	BondOrder  map[[2]int]int // bond order: 1=single,2=double,3=triple,4=aromatic
	BondStereo map[[2]int]int // bond stereo: 0=none,1=up(wedge),6=down(dash)
	Lines      []string       // raw V2000 mol block lines for round-trip serialisation
}

// NewMolecule allocates an empty Molecule with initialised maps.
func NewMolecule(lines []string) *Molecule {
	return &Molecule{
		Adjacency:  make(map[int][]int),
		BondOrder:  make(map[[2]int]int),
		BondStereo: make(map[[2]int]int),
		Lines:      lines,
	}
}

// Reaction holds the parsed reactant and product molecules together with
// the raw RXN header lines needed for round-trip serialisation.
type Reaction struct {
	Reactants []  *Molecule
	Products  []*Molecule
	Header    []string // 5-line RXN header block
}

// Mapping holds the result of atom-to-atom mapping.
// RMap maps [molIdx, atomIdx] in reactants to a map number (1-based).
// PMap maps encodePA(molIdx, atomIdx) in products to the same map number.
type Mapping struct {
	RMap map[[2]int]int
	PMap map[int]int
}

// encodePA encodes a product molecule+atom index pair into a single int.
func encodePA(molIdx, atomIdx int) int { return molIdx*100000 + atomIdx }

// decodePA decodes an encoded product atom back into molecule and atom indices.
func decodePA(enc int) (molIdx, atomIdx int) { return enc / 100000, enc % 100000 }
