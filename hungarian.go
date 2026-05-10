package main

// hungarian solves the rectangular assignment problem using the Hungarian
// algorithm (Kuhn-Munkres). costMatrix[i][j] is the cost of assigning row i
// to column j. Returns assignment[i] = j meaning row i is assigned to col j.
// Unassigned rows get -1. Minimises total cost.
// Complexity: O(n³) where n = max(rows, cols).
func hungarian(costMatrix [][]float64) []int {
	if len(costMatrix) == 0 {
		return nil
	}
	rows := len(costMatrix)
	cols := len(costMatrix[0])

	// Pad to square if needed
	n := rows
	if cols > n {
		n = cols
	}
	inf := 1e18

	// Work on a square padded matrix
	c := make([][]float64, n)
	for i := range c {
		c[i] = make([]float64, n)
		for j := range c[i] {
			if i < rows && j < cols {
				c[i][j] = costMatrix[i][j]
			} else {
				c[i][j] = 0 // padding rows/cols have zero cost
			}
		}
	}

	u := make([]float64, n+1) // potential for rows
	v := make([]float64, n+1) // potential for cols
	p := make([]int, n+1)     // p[j] = row assigned to col j (1-indexed, 0 = unassigned)
	way := make([]int, n+1)   // way[j] = previous col in augmenting path

	for i := 1; i <= n; i++ {
		p[0] = i
		j0 := 0
		minVal := make([]float64, n+1)
		used := make([]bool, n+1)
		for j := range minVal {
			minVal[j] = inf
		}

		for {
			used[j0] = true
			i0 := p[j0]
			delta := inf
			var j1 int

			for j := 1; j <= n; j++ {
				if !used[j] {
					cur := c[i0-1][j-1] - u[i0] - v[j]
					if cur < minVal[j] {
						minVal[j] = cur
						way[j] = j0
					}
					if minVal[j] < delta {
						delta = minVal[j]
						j1 = j
					}
				}
			}

			for j := 0; j <= n; j++ {
				if used[j] {
					u[p[j]] += delta
					v[j] -= delta
				} else {
					minVal[j] -= delta
				}
			}
			j0 = j1
			if p[j0] == 0 {
				break
			}
		}

		for j0 != 0 {
			p[j0] = p[way[j0]]
			j0 = way[j0]
		}
	}

	// Extract assignment: for each row i (1-indexed), find which col it's assigned to
	assignment := make([]int, rows)
	for i := range assignment {
		assignment[i] = -1
	}
	for j := 1; j <= n; j++ {
		if p[j] >= 1 && p[j] <= rows && j <= cols {
			assignment[p[j]-1] = j - 1
		}
	}
	return assignment
}
