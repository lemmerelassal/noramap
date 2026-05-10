package main

import "fmt"

// MappingPipeline orchestrates the full atom-mapping workflow for a single
// reaction file: read → map → verify → write. All dependencies are injected
// as interfaces (D principle), making each stage independently testable.
type MappingPipeline struct {
	Reader   ReactionReader
	Mapper   AtomMapper
	Verifier Verifier
	Writer   ReactionWriter
}

// NewMappingPipeline constructs a MappingPipeline with the standard
// production implementations injected.
func NewMappingPipeline() *MappingPipeline {
	return &MappingPipeline{
		Reader:   NewRXNReader(),
		Mapper:   NewHungarianMapper(),
		Verifier: &SubgraphVerifier{},
		Writer:   &ChemsketchRXNWriter{},
	}
}

// PipelineResult holds the outcome of running the pipeline on one reaction.
type PipelineResult struct {
	Rxn          *Reaction
	Mapping      Mapping
	Verification VerificationResult
}

// Run reads inPath, maps the reaction, verifies, and writes to outPath.
// The mapped file is written regardless of verification outcome.
// Returns an error if reading or writing fails; verification failure is
// non-fatal and reported in PipelineResult.
func (p *MappingPipeline) Run(inPath, outPath string) (*PipelineResult, error) {
	rxn, err := p.Reader.Read(inPath)
	if err != nil {
		return nil, fmt.Errorf("read: %w", err)
	}

	for _, mol := range append(rxn.Reactants, rxn.Products...) {
		if len(mol.Atoms) == 0 {
			return nil, fmt.Errorf("empty mol block in %s", inPath)
		}
	}

	mapping := p.Mapper.Map(rxn)

	// Compute fingerprints for verification (mapper may cache these internally
	// in a production optimisation; for now we recompute via type assertion)
	hm, ok := p.Mapper.(*HungarianMapper)
	var verification VerificationResult
	if ok {
		rFPs := hm.FlatFingerprints(rxn.Reactants)
		pFPs := hm.FlatFingerprints(rxn.Products)
		verification = p.Verifier.Verify(rxn, mapping, rFPs, pFPs)
	}

	if err := p.Writer.Write(outPath, rxn, mapping); err != nil {
		return nil, fmt.Errorf("write: %w", err)
	}

	return &PipelineResult{
		Rxn:          rxn,
		Mapping:      mapping,
		Verification: verification,
	}, nil
}
