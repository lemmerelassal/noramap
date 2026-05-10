package main

import (
	"bufio"
	"fmt"
	"os"
	"os/exec"
	"strings"
)

// OpenBabelHAdder adds explicit hydrogens by shelling out to obabel.
// Implements HydrogenAdder. Fails silently if obabel is unavailable,
// leaving the molecule unchanged (graceful degradation).
type OpenBabelHAdder struct{}

// AddHydrogens implements HydrogenAdder.
func (h *OpenBabelHAdder) AddHydrogens(mol *Molecule) {
	tmp, err := os.CreateTemp("", "mol_*.mol")
	if err != nil {
		return
	}
	tmpPath := tmp.Name()
	defer os.Remove(tmpPath)

	w := bufio.NewWriter(tmp)
	for _, l := range mol.Lines {
		fmt.Fprintln(w, l)
	}
	w.Flush()
	tmp.Close()

	out, err := exec.Command("obabel", tmpPath, "-omol", "-h").Output()
	if err != nil {
		return
	}

	molText := trimObabelOutput(string(out))
	if molText == "" {
		return
	}

	lines := strings.Split(strings.ReplaceAll(molText, "\r\n", "\n"), "\n")
	for len(lines) > 0 && strings.TrimSpace(lines[len(lines)-1]) == "" {
		lines = lines[:len(lines)-1]
	}

	parser := &V2000Parser{}
	newMol, err := parser.Parse(lines)
	if err != nil {
		return
	}

	// Replace mol data in-place (pointer receiver would mutate, but mol is *Molecule)
	mol.Atoms = newMol.Atoms
	mol.Adjacency = newMol.Adjacency
	mol.BondOrder = newMol.BondOrder
	mol.BondStereo = newMol.BondStereo
	mol.Lines = lines
}

// trimObabelOutput strips the trailing "N molecules converted" line that
// obabel appends to stdout, leaving only the mol block.
func trimObabelOutput(output string) string {
	endMarker := "\nM  END"
	idx := strings.LastIndex(output, endMarker)
	if idx < 0 {
		return ""
	}
	return output[:idx+len(endMarker)+1]
}

// NoOpHAdder is a HydrogenAdder that does nothing. Useful in tests or when
// the input mol files already contain explicit hydrogens.
type NoOpHAdder struct{}

// AddHydrogens implements HydrogenAdder (no-op).
func (h *NoOpHAdder) AddHydrogens(_ *Molecule) {}
