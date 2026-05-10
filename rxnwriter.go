package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// ReactionWriter writes a mapped Reaction to a file.
type ReactionWriter interface {
	Write(path string, rxn *Reaction, m Mapping) error
}

// ChemsketchRXNWriter writes V2000 RXN files with CRLF line endings,
// compatible with ACD/ChemSketch and other standard cheminformatics tools.
type ChemsketchRXNWriter struct{}

// Write implements ReactionWriter.
func (w *ChemsketchRXNWriter) Write(path string, rxn *Reaction, m Mapping) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()

	bw := bufio.NewWriter(f)
	defer bw.Flush()

	crlf := func(s string) { fmt.Fprintf(bw, "%s\r\n", s) }

	// Write RXN header (Chemsketch format: $RXN, blank, program, blank, counts)
	crlf("$RXN")
	crlf("")
	if len(rxn.Header) > 2 {
		crlf(rxn.Header[2])
	} else {
		crlf("")
	}
	crlf("")
	if len(rxn.Header) > 4 {
		crlf(rxn.Header[4])
	}

	for mi, mol := range rxn.Reactants {
		crlf("$MOL")
		amn := buildAtomMapNums(mol, mi, m.RMap)
		writeMolBlock(bw, mol, amn)
	}
	for mi, mol := range rxn.Products {
		crlf("$MOL")
		amn := buildProductMapNums(mol, mi, m.PMap)
		writeMolBlock(bw, mol, amn)
	}
	return nil
}

func buildAtomMapNums(mol *Molecule, molIdx int, rMap map[[2]int]int) map[int]int {
	amn := make(map[int]int)
	for ai := range mol.Atoms {
		if mn, ok := rMap[[2]int{molIdx, ai}]; ok {
			amn[ai] = mn
		}
	}
	return amn
}

func buildProductMapNums(mol *Molecule, molIdx int, pMap map[int]int) map[int]int {
	amn := make(map[int]int)
	for ai := range mol.Atoms {
		enc := encodePA(molIdx, ai)
		if mn, ok := pMap[enc]; ok {
			amn[ai] = mn
		}
	}
	return amn
}

func writeMolBlock(w *bufio.Writer, mol *Molecule, atomMapNums map[int]int) {
	crlf := func(s string) { fmt.Fprintf(w, "%s\r\n", s) }
	lines := mol.Lines
	if len(lines) < 4 {
		return
	}

	cf := strings.Fields(lines[3])
	numAtoms, _ := strconv.Atoi(cf[0])
	numBonds, _ := strconv.Atoi(cf[1])

	// Preserve original header lines (0–3)
	for i := 0; i < 4 && i < len(lines); i++ {
		crlf(lines[i])
	}

	for i := 0; i < numAtoms; i++ {
		li := 4 + i
		if li >= len(lines) {
			break
		}
		crlf(injectAtomMap(lines[li], atomMapNums[i]))
	}
	for i := 0; i < numBonds; i++ {
		li := 4 + numAtoms + i
		if li >= len(lines) {
			break
		}
		crlf(lines[li])
	}
	for i := 4 + numAtoms + numBonds; i < len(lines); i++ {
		crlf(lines[i])
		if strings.HasPrefix(lines[i], "M  END") {
			break
		}
	}
}

// injectAtomMap writes the map number into the correct fixed-width position
// in a V2000 atom line. The atom-atom mapping occupies chars 61–63 (1-based).
func injectAtomMap(line string, mapNum int) string {
	for len(line) < 69 {
		line += " "
	}
	return fmt.Sprintf("%s%3d%s", line[:60], mapNum, line[63:])
}
