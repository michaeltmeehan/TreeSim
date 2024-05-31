# TreeSim
TreeSim is a Julia package designed for simulating phylogenetic trees from transmission linelist data. The package provides tools to generate within-host trees, combine them into a single tree, and calculate probabilities and constraints associated with phylogenetic processes.

## Features
- Define and manipulate host structures using the Host type.
- Simulate within-host phylogenetic trees.
- Combine multiple within-host trees into a single phylogenetic tree.
- Calculate probabilities and enforce constraints for coalescent events.
- Normalize and relabel phylogenetic trees for analysis.

## Installation
To install TreeSim, use the following command:
```julia
using Pkg
Pkg.add("https://github.com/michaeltmeehan/TreeSim.jl")
```

## Usage
Here is an example of how to use TreeSim in conjunction with the EpiSim package to simulate a phylogenetic tree:
```julia
using TreeSim
using EpiSim
using DataFrames

# Define the parameters for the epidemiological simulation
params = CRBDParameters(N₀=1, λ=0.5, μ=0.2, ψ=0.1, ρ₀=0.1, r=0.1, t_max=100.0)

# Simulate the epidemiological process
out = simulate_outbreak(params, N_max = 50, S_max = 10)

# Generate the phylogenetic tree from the linelist
tree = simulate_phylogeny(out)

# Display the resulting tree
println(tree)
```

## Main Functions
### `simulate_phylogeny`
Simulate a phylogenetic tree from a linelist DataFrame.
```julia
simulate_phylogeny(linelist::DataFrame; Nₑ::Float64=1e-6, binarize=false) -> Phylogeny
```
#### Arguments
- `linelist::DataFrame`: A `DataFrame` representing the history of transmission of an outbreak. Each row should include columns for the infected individual's ID (child_id), the ID of the person that infected them (parent_id), the time of infection (t_birth), the time of sampling (t_sam), and additional columns like child_type and parent_type.
- `Nₑ::Float64`: The effective population size (default is 1e-6).

#### Returns
- `DataFrame`: A DataFrame representing the combined phylogenetic tree, containing columns for times (t), node IDs (id), left and right child IDs (left, right), leaf IDs (leaf_id), host IDs (host), and node types (type).

## Citation
If you use TreeSim in your research, please cite the following paper related to the bounded coalescent model theory and implementation:

Carson J, Ledda A, Ferretti L, Keeling MJ, Didelot X. The bounded coalescent model: Conditioning a genealogy on a minimum root date. J Theor Biol. 2022 Sep 7;548:111186. [DOI: 10.1016/j.jtbi.2022.111186](https://doi.org/10.1016/j.jtbi.2022.111186). Epub 2022 Jun 10. PMID: 35697144.

This paper describes a new algorithm to simulate from the bounded coalescent model directly, which is more computationally efficient than the rejection sampling approach.

## License
This package is licensed under the MIT License.

## Contributing
Contributions are welcome! Please feel free to submit a pull request or open an issue on GitHub.