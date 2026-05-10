package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// ReactionReader is the interface for reading a Reaction from a source.
type ReactionReader interface {
	Read(path string) (*Reaction, error)
}

// HydrogenAdder is the interface for adding explicit hydrogens to a Molecule.
// Injected into RXNReader so the hydrogen strategy is swappable (O/D).
type HydrogenAdder interface {
	AddHydrogens(mol *Molecule)
}

// RXNReader reads V2000 RXN files, parsing each mol block via an injected
// MolParser and adding explicit hydrogens via an injected HydrogenAdder.
type RXNReader struct {
	Parser    MolParser
	HAdder    HydrogenAdder
}

// NewRXNReader constructs an RXNReader with the standard V2000 parser and
// the OpenBabel hydrogen adder.
func NewRXNReader() *RXNReader {
	return &RXNReader{
		Parser: &V2000Parser{},
		HAdder: &OpenBabelHAdder{},
	}
}

// Read implements ReactionReader.
func (r *RXNReader) Read(path string) (*Reaction, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var all []string
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		all = append(all, sc.Text())
	}

	rxnStart := -1
	for i, l := range all {
		if strings.HasPrefix(l, "$RXN") {
			rxnStart = i
			break
		}
	}
	if rxnStart < 0 {
		return nil, fmt.Errorf("no $RXN found in %s", path)
	}

	// Chemsketch RXN header: $RXN, blank, program/date, blank, counts (5 lines)
	if rxnStart+5 > len(all) {
		return nil, fmt.Errorf("RXN header truncated in %s", path)
	}
	header := make([]string, 5)
	copy(header, all[rxnStart:rxnStart+5])

	cf := strings.Fields(all[rxnStart+4])
	if len(cf) < 2 {
		return nil, fmt.Errorf("bad RXN counts line in %s", path)
	}
	nR, _ := strconv.Atoi(cf[0])
	nP, _ := strconv.Atoi(cf[1])

	cursor := rxnStart + 5
	reactants, err := r.collectMols(all, &cursor, nR)
	if err != nil {
		return nil, fmt.Errorf("reactants: %w", err)
	}
	products, err := r.collectMols(all, &cursor, nP)
	if err != nil {
		return nil, fmt.Errorf("products: %w", err)
	}

	return &Reaction{
		Reactants: reactants,
		Products:  products,
		Header:    header,
	}, nil
}

func (r *RXNReader) collectMols(all []string, cursor *int, n int) ([]*Molecule, error) {
	var mols []*Molecule
	for i := 0; i < n; i++ {
		// Advance to next $MOL marker
		for *cursor < len(all) && !strings.HasPrefix(all[*cursor], "$MOL") {
			*cursor++
		}
		if *cursor >= len(all) {
			return nil, fmt.Errorf("expected $MOL #%d, reached EOF", i+1)
		}
		*cursor++ // skip $MOL line

		start := *cursor
		for *cursor < len(all) &&
			!strings.HasPrefix(all[*cursor], "$MOL") &&
			!strings.HasPrefix(all[*cursor], "$RXN") &&
			!strings.HasPrefix(all[*cursor], "$$$$") {
			*cursor++
		}

		mol, err := r.Parser.Parse(all[start:*cursor])
		if err != nil {
			return nil, fmt.Errorf("mol %d: %w", i+1, err)
		}
		r.HAdder.AddHydrogens(mol)
		mols = append(mols, mol)
	}
	return mols, nil
}
