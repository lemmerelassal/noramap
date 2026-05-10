package main

import (
	"fmt"
	"math"
	"os"
	"sort"
	"strings"
)

func main() {
	if len(os.Args) < 2 {
		printUsage()
		os.Exit(1)
	}

	switch os.Args[1] {
	case "kegg":
		runKEGG()
	default:
		runSingle()
	}
}

func runKEGG() {
	startFrom := ""
	if len(os.Args) >= 3 && strings.HasPrefix(os.Args[2], "resume:") {
		startFrom = strings.TrimPrefix(os.Args[2], "resume:")
		fmt.Printf("Resuming from %s\n", startFrom)
	}
	pipeline := NewMappingPipeline()
	if _, err := DownloadKEGG(pipeline, startFrom); err != nil {
		fmt.Fprintf(os.Stderr, "KEGG error: %v\n", err)
		os.Exit(1)
	}
}

func runSingle() {
	if len(os.Args) < 3 {
		fmt.Fprintln(os.Stderr, "Usage: atom_map <input.rxn> <output.rxn>")
		os.Exit(1)
	}

	pipeline := NewMappingPipeline()
	result, err := pipeline.Run(os.Args[1], os.Args[2])
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	printMappingTable(result)
	PrintVerification(result.Verification)
	fmt.Printf("\nWrote: %s\n", os.Args[2])
}

func printMappingTable(result *PipelineResult) {
	rxn := result.Rxn
	m := result.Mapping

	fmt.Printf("Reactants: %d  Products: %d\n", len(rxn.Reactants), len(rxn.Products))
	fmt.Printf("\n%-4s %-6s %-5s  %-4s %-6s %-5s  %-10s %-10s %-8s  %s\n",
		"RMol", "RAtom", "Elem", "PMol", "PAtom", "Elem", "R.FP", "P.FP", "Dist", "MapNum")
	fmt.Println(strings.Repeat("-", 80))

	rows := buildRows(rxn, m)
	sort.Slice(rows, func(i, j int) bool { return rows[i].mapNum < rows[j].mapNum })

	for _, r := range rows {
		suffix := ""
		if r.fallback {
			suffix = " ~"
		}
		mapStr := fmt.Sprintf("%d%s", r.mapNum, suffix)
		if r.mapNum == 0 {
			mapStr = "unmapped"
		}
		fmt.Printf("R%-3d A%-5d %-5s  P%-3d A%-5d %-5s  %-10.4f %-10.4f %-8.4f  %s\n",
			r.rmi+1, r.rai+1, r.sym,
			r.pmi+1, r.pai+1, r.sym,
			r.rNorm, r.pNorm, r.dist,
			mapStr)
	}
	fmt.Println("(~ = fallback: second-pass assignment, large environment change)")
}

type tableRow struct {
	rmi, rai int
	pmi, pai int
	sym      string
	rNorm    float64
	pNorm    float64
	dist     float64
	mapNum   int
	fallback bool
}

func buildRows(rxn *Reaction, m Mapping) []tableRow {
	hm := NewHungarianMapper()
	pFlatAll := hm.flatten(rxn.Products)

	mapNumToEnc := make(map[int]int)
	for enc, mn := range m.PMap {
		mapNumToEnc[mn] = enc
	}

	var rows []tableRow
	for mi, mol := range rxn.Reactants {
		fps := hm.Fingerprinter.Compute(mol)
		for ai, fp := range fps {
			key := [2]int{mi, ai}
			mapNum := m.RMap[key]

			pmi, pai := 0, 0
			var pFP []float64
			isFallback := false

			if assignedEnc, ok := mapNumToEnc[mapNum]; ok && mapNum != 0 {
				pmi, pai = decodePA(assignedEnc)
				if pmi < len(rxn.Products) && pai < len(rxn.Products[pmi].Atoms) {
					pFP = hm.Fingerprinter.Compute(rxn.Products[pmi])[pai]
				}
				isFallback = true
				for _, pa := range pFlatAll {
					if pa.Symbol == mol.Atoms[ai].Symbol && encodePA(pa.MolIdx, pa.AtomIdx) == assignedEnc {
						isFallback = false
						break
					}
				}
			} else {
				pFP = make([]float64, maxLevel)
			}

			dist := 0.0
			if len(pFP) > 0 {
				dist = vecDist(fp, pFP)
			}

			rows = append(rows, tableRow{
				rmi: mi, rai: ai,
				pmi: pmi, pai: pai,
				sym:      mol.Atoms[ai].Symbol,
				rNorm:    vecNorm(fp),
				pNorm:    vecNorm(pFP),
				dist:     dist,
				mapNum:   mapNum,
				fallback: isFallback,
			})
		}
	}
	return rows
}

func printUsage() {
	fmt.Fprintln(os.Stderr, "Usage:")
	fmt.Fprintln(os.Stderr, "  atom_map <input.rxn> <output.rxn>      map a single reaction")
	fmt.Fprintln(os.Stderr, "  atom_map kegg [resume:RXXXXX]          download and map all KEGG reactions")
}

// ensure math is used (vecNorm uses it via fingerprint.go, but keep import clean)
var _ = math.Sqrt
