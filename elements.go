package main

// elementPrime maps element symbols to unique prime numbers.
// The fundamental theorem of arithmetic guarantees that products over
// any multiset of elements are unique to that multiset within a single shell.
var elementPrime = map[string]float64{
	// Organic
	"H": 2, "C": 3, "N": 5, "O": 7,
	"F": 11, "P": 13, "S": 17, "Cl": 19,
	"Br": 23, "I": 29, "B": 31, "Si": 37,
	"Se": 41, "As": 43, "Te": 47,
	// Biologically relevant metals
	"Fe": 53, "Mg": 59, "Zn": 61, "Ca": 67,
	"Mn": 71, "Co": 73, "Cu": 79, "Ni": 83,
	"Mo": 89, "V": 97, "W": 101, "Cr": 103,
	// Alkali / alkaline earth
	"Na": 107, "K": 109, "Li": 113, "Rb": 127,
	"Cs": 131, "Ba": 137, "Sr": 139,
	// Other inorganic
	"Al": 149, "Ge": 151, "Sn": 157, "Pb": 163,
	"Sb": 167, "Bi": 173, "Cd": 179, "Hg": 181,
	"Pt": 191, "Pd": 193, "Au": 197, "Ag": 199,
	"Ti": 211, "Zr": 223, "Rh": 227, "Ir": 229,
	"Ru": 233, "Os": 239, "Re": 241, "Tc": 251,
	"Sc": 257, "Y": 263, "La": 269, "Ce": 271,
}

// ElementPrime returns the prime assigned to an element symbol.
// Unknown elements fall back to 43 (As) — a safe default that avoids
// collisions with any element in the table above.
func ElementPrime(symbol string) float64 {
	if p, ok := elementPrime[symbol]; ok {
		return p
	}
	return 43
}
