package main

import (
	"bufio"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"sync"
	"sync/atomic"
	"time"
)

const (
	keggBase        = "https://rest.kegg.jp"
	keggDelay       = 400 * time.Millisecond
	inputDir        = "input"
	outputDir       = "output"
	downloadWorkers = 3
	mapWorkers      = 8
)

// KEGGClient handles all communication with the KEGG REST API.
// Isolated here so it can be mocked in tests (I principle).
type KEGGClient struct {
	rateTok <-chan struct{}
}

func (c *KEGGClient) get(path string) (string, error) {
	if c.rateTok != nil {
		<-c.rateTok
	}
	url := keggBase + path
	resp, err := http.Get(url)
	if err != nil {
		return "", fmt.Errorf("GET %s: %w", url, err)
	}
	defer resp.Body.Close()
	if resp.StatusCode != 200 {
		return "", fmt.Errorf("GET %s: status %d", url, resp.StatusCode)
	}
	b, err := io.ReadAll(resp.Body)
	return string(b), err
}

func (c *KEGGClient) listReactions() ([]string, error) {
	body, err := c.get("/list/reaction")
	if err != nil {
		return nil, err
	}
	var ids []string
	sc := bufio.NewScanner(strings.NewReader(body))
	for sc.Scan() {
		fields := strings.Fields(sc.Text())
		if len(fields) == 0 {
			continue
		}
		id := strings.TrimPrefix(fields[0], "rn:")
		if strings.HasPrefix(id, "R") {
			ids = append(ids, id)
		}
	}
	return ids, nil
}

func (c *KEGGClient) reactionEquation(rxnID string) (reactants, products []string, equation string, err error) {
	body, err := c.get("/get/" + rxnID)
	if err != nil {
		return nil, nil, "", err
	}

	sc := bufio.NewScanner(strings.NewReader(body))
	var eqLines []string
	inEq := false
	for sc.Scan() {
		line := sc.Text()
		if strings.HasPrefix(line, "EQUATION") {
			inEq = true
			eqLines = append(eqLines, strings.TrimPrefix(line, "EQUATION"))
		} else if inEq {
			if len(line) > 0 && line[0] == ' ' {
				eqLines = append(eqLines, line)
			} else {
				break
			}
		}
	}
	if len(eqLines) == 0 {
		return nil, nil, "", fmt.Errorf("no EQUATION in %s", rxnID)
	}

	equation = strings.TrimSpace(strings.Join(eqLines, " "))
	sides := strings.SplitN(equation, "<=>", 2)
	if len(sides) != 2 {
		return nil, nil, equation, fmt.Errorf("no <=> in equation: %s", equation)
	}

	return parseStoichiometry(sides[0]), parseStoichiometry(sides[1]), equation, nil
}

func (c *KEGGClient) fetchMol(compoundID string) (string, error) {
	body, err := c.get("/get/" + compoundID + "/mol")
	if err != nil {
		return "", err
	}
	if strings.TrimSpace(body) == "" || strings.Contains(body, "No such data") {
		return "", nil
	}
	return body, nil
}

// parseStoichiometry expands coefficients: "2 C00001 + C00002" → ["C00001","C00001","C00002"]
func parseStoichiometry(side string) []string {
	reCompound := regexp.MustCompile(`C\d{5}`)
	reCoeff := regexp.MustCompile(`(\d+)\s+(C\d{5})`)

	withCoeff := make(map[string]bool)
	var result []string
	for _, m := range reCoeff.FindAllStringSubmatch(side, -1) {
		coeff := 1
		fmt.Sscanf(m[1], "%d", &coeff)
		if coeff > 10 {
			coeff = 10
		}
		for i := 0; i < coeff; i++ {
			result = append(result, m[2])
		}
		withCoeff[m[2]] = true
	}
	for _, cid := range reCompound.FindAllString(side, -1) {
		if !withCoeff[cid] {
			result = append(result, cid)
		}
	}
	return result
}

// buildRXNContent assembles a CRLF RXN file string from raw mol blocks.
func buildRXNContent(rxnID string, rMols, pMols []string) string {
	crlf := "\r\n"
	var sb strings.Builder
	sb.WriteString("$RXN" + crlf)
	sb.WriteString(rxnID + crlf)
	sb.WriteString("  KEGG" + crlf)
	sb.WriteString(crlf)
	sb.WriteString(fmt.Sprintf("%3d%3d", len(rMols), len(pMols)) + crlf)

	writeMol := func(mol string) {
		sb.WriteString("$MOL" + crlf)
		raw := strings.Split(strings.ReplaceAll(mol, "\r\n", "\n"), "\n")
		var lines []string
		for _, l := range raw {
			lines = append(lines, strings.TrimRight(l, " \t\r"))
			if strings.HasPrefix(l, "M  END") {
				break
			}
		}
		for len(lines) < 4 {
			lines = append([]string{""}, lines...)
		}
		for _, l := range lines {
			sb.WriteString(l + crlf)
		}
	}
	for _, m := range rMols {
		writeMol(m)
	}
	for _, m := range pMols {
		writeMol(m)
	}
	return sb.String()
}

// newRateLimiter returns a token channel firing every keggDelay.
func newRateLimiter(stop <-chan struct{}) <-chan struct{} {
	ch := make(chan struct{}, downloadWorkers)
	go func() {
		ticker := time.NewTicker(keggDelay)
		defer ticker.Stop()
		for {
			select {
			case <-stop:
				return
			case <-ticker.C:
				select {
				case ch <- struct{}{}:
				default:
				}
			}
		}
	}()
	return ch
}

// MolCache is a thread-safe cache of compound ID → mol block string.
type MolCache struct {
	mu   sync.Mutex
	data map[string]string
}

func (mc *MolCache) Get(id string) (string, bool) {
	mc.mu.Lock()
	defer mc.mu.Unlock()
	v, ok := mc.data[id]
	return v, ok
}

func (mc *MolCache) Set(id, mol string) {
	mc.mu.Lock()
	defer mc.mu.Unlock()
	mc.data[id] = mol
}

// DownloadKEGG downloads and immediately maps all KEGG reactions.
func DownloadKEGG(pipeline *MappingPipeline, startFrom string) ([]string, error) {
	if err := os.MkdirAll(inputDir, 0755); err != nil {
		return nil, fmt.Errorf("mkdir %s: %w", inputDir, err)
	}
	if err := os.MkdirAll(outputDir, 0755); err != nil {
		return nil, fmt.Errorf("mkdir %s: %w", outputDir, err)
	}

	stopRL := make(chan struct{})
	rateTok := newRateLimiter(stopRL)
	client := &KEGGClient{rateTok: rateTok}
	cache := &MolCache{data: make(map[string]string)}

	fmt.Println("Fetching reaction list from KEGG...")
	ids, err := client.listReactions()
	if err != nil {
		return nil, fmt.Errorf("list reactions: %w", err)
	}
	fmt.Printf("Found %d reactions\n", len(ids))

	type job struct {
		idx, total int
		rxnID      string
	}

	jobCh := make(chan job, downloadWorkers*2)
	mapCh := make(chan string, mapWorkers*2)

	var downloaded, skipped, failed, mapped, mapFailed int64

	var dlWg sync.WaitGroup
	for w := 0; w < downloadWorkers; w++ {
		dlWg.Add(1)
		go func() {
			defer dlWg.Done()
			for j := range jobCh {
				outPath := filepath.Join(inputDir, j.rxnID+".rxn")
				if _, err := os.Stat(outPath); err == nil {
					mapCh <- j.rxnID
					atomic.AddInt64(&skipped, 1)
					continue
				}

				fmt.Printf("[%d/%d] %s downloading...\n", j.idx+1, j.total, j.rxnID)

				rIDs, pIDs, _, err := client.reactionEquation(j.rxnID)
				if err != nil || len(rIDs) == 0 || len(pIDs) == 0 {
					fmt.Printf("[%d/%d] %s skip (no equation)\n", j.idx+1, j.total, j.rxnID)
					atomic.AddInt64(&failed, 1)
					continue
				}

				fetchMols := func(cids []string) ([]string, bool) {
					var mols []string
					for _, cid := range cids {
						mol, ok := cache.Get(cid)
						if !ok {
							mol, err = client.fetchMol(cid)
							if err != nil || mol == "" {
								fmt.Printf("[%d/%d] %s skip (no mol %s)\n", j.idx+1, j.total, j.rxnID, cid)
								return nil, false
							}
							cache.Set(cid, mol)
						}
						mols = append(mols, mol)
					}
					return mols, true
				}

				rMols, ok := fetchMols(rIDs)
				if !ok {
					atomic.AddInt64(&failed, 1)
					continue
				}
				pMols, ok := fetchMols(pIDs)
				if !ok {
					atomic.AddInt64(&failed, 1)
					continue
				}

				content := buildRXNContent(j.rxnID, rMols, pMols)
				if err := os.WriteFile(outPath, []byte(content), 0644); err != nil {
					atomic.AddInt64(&failed, 1)
					continue
				}
				atomic.AddInt64(&downloaded, 1)
				mapCh <- j.rxnID
			}
		}()
	}

	var mapWg sync.WaitGroup
	var resultIDs struct {
		sync.Mutex
		ids []string
	}

	for w := 0; w < mapWorkers; w++ {
		mapWg.Add(1)
		go func() {
			defer mapWg.Done()
			for rxnID := range mapCh {
				inPath := filepath.Join(inputDir, rxnID+".rxn")
				outPath := filepath.Join(outputDir, rxnID+".rxn")

				if _, err := os.Stat(outPath); err == nil {
					resultIDs.Lock()
					resultIDs.ids = append(resultIDs.ids, rxnID)
					resultIDs.Unlock()
					continue
				}

				_, err := pipeline.Run(inPath, outPath)
				if err != nil {
					fmt.Printf("  map %s: %v\n", rxnID, err)
					atomic.AddInt64(&mapFailed, 1)
				} else {
					fmt.Printf("  mapped %s\n", rxnID)
					atomic.AddInt64(&mapped, 1)
				}

				resultIDs.Lock()
				resultIDs.ids = append(resultIDs.ids, rxnID)
				resultIDs.Unlock()
			}
		}()
	}

	skip := startFrom != ""
	total := len(ids)
	for i, rxnID := range ids {
		if skip {
			if rxnID == startFrom {
				skip = false
			} else {
				continue
			}
		}
		jobCh <- job{i, total, rxnID}
	}
	close(jobCh)
	dlWg.Wait()
	close(stopRL)
	close(mapCh)
	mapWg.Wait()

	fmt.Printf("\nDownloaded: %d  Skipped: %d  Failed: %d  Mapped: %d  MapFailed: %d\n",
		atomic.LoadInt64(&downloaded), atomic.LoadInt64(&skipped),
		atomic.LoadInt64(&failed), atomic.LoadInt64(&mapped),
		atomic.LoadInt64(&mapFailed))

	return resultIDs.ids, nil
}
