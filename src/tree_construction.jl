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


function simulate(linelist::DataFrame, Nₑ::Float64)::DataFrame
    hosts = get_hosts(linelist)
    # Generate within-host trees for each sampled host in linelist
    wtrees = [sample_wtree(host, Nₑ) for host in values(hosts) if host.id != 0]
    return join_wtrees(wtrees)
end


# function simulate(proc::BDProcess, Nₑ::Float64)::DataFrame
#     return simulate(proc.linelist, Nₑ)
# end