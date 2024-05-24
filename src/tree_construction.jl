"""
    get_hosts(linelist::DataFrame) -> Dict{Int64, Host}

Construct a dictionary of `Host` objects from a linelist DataFrame, representing the transmission history of an outbreak.

# Arguments
- `linelist::DataFrame`: A DataFrame representing the history of transmission of an outbreak. Each row should include columns for the infected individual's ID (`child_id`), the ID of the person that infected them (`parent_id`), the time of infection (`t_birth`), and the time of sampling (`t_sam`), if they were sampled.

# Returns
- `Dict{Int64, Host}`: A dictionary where keys correspond to unique host IDs and values are `Host` objects representing each host in the transmission chain.
"""
function get_hosts(linelist::DataFrame)
	hosts = Dict{Int64, Host}()
	@eachrow! reverse(linelist) begin
		if :t_sam > 0. || haskey(hosts, :child_id)
			# Check if a record exists for the infector
			if haskey(hosts, :parent_id) # Append existing record
				pushfirst!(hosts[:parent_id].leaf_times, :t_birth)
				pushfirst!(hosts[:parent_id].leaf_ids, :child_id)
			else   # Create new record
				hosts[:parent_id] = Host(:parent_id, :parent_type, nothing, 0, NaN, [:t_birth], [:child_id])
			end

			if :t_sam > 0. && !haskey(hosts, :child_id)	# Sampled infectee
				hosts[:child_id] = Host(:child_id, :child_type, :parent_id, :parent_type, :t_birth, [:t_sam], [:child_id])	# Create a new record for the infectee
			
			elseif haskey(hosts, :child_id)	# Infector
				hosts[:child_id].root = :t_birth	# Update root
                hosts[:child_id].infector = :parent_id  # Update infector
                hosts[:child_id].type = :child_type # Update type
                hosts[:child_id].infector_type = :parent_type
				
				if :t_sam > 0.  # TODO: Generalize to sampling without removal (i.e., the ordered placement of leaves)
					push!(hosts[:child_id].leaf_times, :t_sam)	# Update leaves to include sampling
					push!(hosts[:child_id].leaf_ids, :child_id)
				end
			end
		end
	end
	return hosts
end


@forward Outbreak.linelist get_hosts


"""
    simulate_phylogeny(linelist::DataFrame; Nₑ::Float64=1e-6) -> DataFrame

Simulate a phylogenetic tree from a linelist DataFrame, generating within-host trees for each sampled host and combining them into a single tree.

# Arguments
- `linelist::DataFrame`: A DataFrame representing the history of transmission of an outbreak. Each row should include columns for the infected individual's ID (`child_id`), the ID of the person that infected them (`parent_id`), the time of infection (`t_birth`), the time of sampling (`t_sam`), and additional columns like `child_type` and `parent_type`.
- `Nₑ::Float64`: The effective population size (default is `1e-6`).

# Returns
- `nothing` : If linelist is empty or only contains a single individual, simulate_phylogeny returns `nothing`
- `DataFrame`: A DataFrame representing the combined phylogenetic tree, containing columns for times (`t`), node IDs (`id`), left and right child IDs (`left`, `right`), leaf IDs (`leaf_id`), host IDs (`host`), and node types (`type`).
"""
function simulate_phylogeny(linelist::DataFrame; Nₑ::Float64=1e-6, binarize=true)::Union{DataFrame, Nothing}
    nrow(linelist) ≤ 1 && return nothing
    hosts = get_hosts(linelist)
    sum(linelist[:, :parent_id] .== 0) != 1 && throw(ArgumentError("simulate_phylogeny cannot handle outbreaks with multiple seed cases"))
    # Generate within-host trees for each sampled host in linelist
    wtrees = [sample_wtree(host, Nₑ) for host in values(hosts) if host.id != 0]
    tree = join_wtrees(wtrees)
    if binarize
        tree = binarize_tree(tree)
    end
    return tree
end


@forward Outbreak.linelist simulate_phylogeny


function Base.convert(::Type{RecursiveTree}, df::DataFrame)
    # df_binary = binarize_tree(df)
    tree = RecursiveTree{OneRoot, String, Dict{String, Any}, Dict{String, Any}, PolytomousBranching, Float64,Dict{String, Any}}()
    @eachrow! df begin
        lab = string(:id)
        createnode!(tree, lab)
        setnodedata!(tree, lab, "id", :id)
        setnodedata!(tree, lab, "t", :t)
        setnodedata!(tree, lab, "host", :host)
        setnodedata!(tree, lab, "type", :type)
        :left != 0 && createbranch!(tree, lab, string(:left), df[:left, :].t - :t)
        :right != 0 && createbranch!(tree, lab, string(:right), df[:right, :].t - :t)
    end
    return tree
end


RecursiveTree(df::DataFrame) = convert(RecursiveTree, df)


function RecursiveTree(out::Outbreak; Nₑ=1e-6, binarize=true)
    tree = simulate_phylogeny(out, Nₑ=Nₑ, binarize=binarize)
    return convert(RecursiveTree, tree)
end


mutable struct Phylogeny
    tree::RecursiveTree
    Nₑ::Float64
end


function Phylogeny(out::Outbreak; Nₑ=1e-6, binarize=true)
    return Phylogeny(RecursiveTree(out, Nₑ=Nₑ, binarize=binarize), Nₑ)
end


@forward Phylogeny.tree RecursiveTree