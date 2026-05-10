package main

import (
	"fmt"
	"strconv"
	"strings"
)

// MolParser is the interface for parsing V2000 mol block lines into a Molecule.
// Separating the interface from the implementation (O: open for extension)
// allows alternative parsers (e.g. V3000, SDF) to be substituted (L: Liskov).
type MolParser interface {
	Parse(lines []string) (*Molecule, error)
}

// V2000Parser parses MDL V2000 mol block lines.
type V2000Parser struct{}

// Parse implements MolParser.
func (p *V2000Parser) Parse(lines []string) (*Molecule, error) {
	if len(lines) < 4 {
		return nil, fmt.Errorf("mol block too short (%d lines)", len(lines))
	}

	counts := strings.Fields(lines[3])
	if len(counts) < 2 {
		return nil, fmt.Errorf("invalid counts line: %q", lines[3])
	}
	numAtoms, err := strconv.Atoi(counts[0])
	if err != nil {
		return nil, fmt.Errorf("bad atom count: %w", err)
	}
	numBonds, err := strconv.Atoi(counts[1])
	if err != nil {
		return nil, fmt.Errorf("bad bond count: %w", err)
	}

	mol := NewMolecule(lines)

	if err := p.parseAtoms(mol, lines, numAtoms); err != nil {
		return nil, err
	}
	if err := p.parseBonds(mol, lines, numAtoms, numBonds); err != nil {
		return nil, err
	}
	p.parseChargeLines(mol, lines)

	return mol, nil
}

func (p *V2000Parser) parseAtoms(mol *Molecule, lines []string, numAtoms int) error {
	// V2000 charge code → formal charge value
	chargeDecodeMap := map[int]int{0: 0, 1: 3, 2: 2, 3: 1, 4: 0, 5: -1, 6: -2, 7: -3}

	for i := 0; i < numAtoms; i++ {
		li := 4 + i
		if li >= len(lines) {
			return fmt.Errorf("unexpected end in atom block at line %d", li)
		}
		f := strings.Fields(lines[li])
		if len(f) < 4 {
			return fmt.Errorf("short atom line: %q", lines[li])
		}

		x, _ := strconv.ParseFloat(f[0], 64)
		y, _ := strconv.ParseFloat(f[1], 64)
		z, _ := strconv.ParseFloat(f[2], 64)

		stereo := 0
		if len(f) > 6 {
			stereo, _ = strconv.Atoi(f[6])
		}
		chargeCode := 0
		if len(f) > 5 {
			chargeCode, _ = strconv.Atoi(f[5])
		}

		mol.Atoms = append(mol.Atoms, Atom{
			Index:        i,
			Symbol:       f[3],
			X:            x,
			Y:            y,
			Z:            z,
			StereoParity: stereo,
			FormalCharge: chargeDecodeMap[chargeCode],
		})
	}
	return nil
}

func (p *V2000Parser) parseBonds(mol *Molecule, lines []string, numAtoms, numBonds int) error {
	for i := 0; i < numBonds; i++ {
		li := 4 + numAtoms + i
		if li >= len(lines) {
			break
		}
		f := strings.Fields(lines[li])
		if len(f) < 2 {
			continue
		}

		a1, _ := strconv.Atoi(f[0])
		a2, _ := strconv.Atoi(f[1])
		a1--
		a2-- // convert to 0-based

		mol.Adjacency[a1] = append(mol.Adjacency[a1], a2)
		mol.Adjacency[a2] = append(mol.Adjacency[a2], a1)

		bo := 1
		if len(f) > 2 {
			bo, _ = strconv.Atoi(f[2])
		}
		mol.BondOrder[[2]int{a1, a2}] = bo
		mol.BondOrder[[2]int{a2, a1}] = bo

		if len(f) > 3 {
			bs, _ := strconv.Atoi(f[3])
			if bs != 0 {
				// Directional: only a1→a2 (wedge direction is relative to bond source)
				mol.BondStereo[[2]int{a1, a2}] = bs
			}
		}
	}
	return nil
}

// parseChargeLines reads M  CHG property lines and overrides charge codes.
// M  CHG is authoritative over the atom-line charge code in modern mol files.
func (p *V2000Parser) parseChargeLines(mol *Molecule, lines []string) {
	for _, l := range lines {
		if !strings.HasPrefix(l, "M  CHG") {
			continue
		}
		fields := strings.Fields(l)
		if len(fields) < 4 {
			continue
		}
		n, _ := strconv.Atoi(fields[2])
		for i := 0; i < n && 3+2*i+1 < len(fields); i++ {
			atomIdx, _ := strconv.Atoi(fields[3+2*i])
			charge, _ := strconv.Atoi(fields[3+2*i+1])
			atomIdx-- // 1-based in mol file
			if atomIdx >= 0 && atomIdx < len(mol.Atoms) {
				mol.Atoms[atomIdx].FormalCharge = charge
			}
		}
	}
}
