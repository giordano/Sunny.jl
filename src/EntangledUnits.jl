# `forward` contains a list from sites of the original crystal to a site of the
# contracted crystal, including an extra index to keep track entangled units:
# `(contracted_crystal_site_index, intra_unit_site_index)`. If the first index
# refers to a site in the new crystal that does not contain multiple units, than
# the second index will always be 1. 
#
# `inverse` contains a list of length equal to the number of sites in the
# contracted crystal (corresponding to `contracted_crystal_site_index` above).
# Each element of this list is another list of tuples,
# `(site_of_original_crystal, position_offset)`. The position offset is applied
# to the position of the contracted crystal to recover the corresponding
# location in the original crystal. The length of these sublists corresponds to
# the number of sites within the entangled unit.
struct CrystalContractionInfo
    forward :: Vector{Tuple{Int64, Int64}}
    inverse :: Vector{Vector{Tuple{Int64, Vec3}}}
end

function contract_crystal(crystal, units)

    # Determine which sites are entangled and which are not.
    unentangled_sites = 1:natoms(crystal)
    entangled_sites = Int64[]
    for unit in units, site in unit
        push!(entangled_sites, site)
        unentangled_sites = filter(!=(site), unentangled_sites)
    end

    # Check that unit definitions are valid.
    @assert length(entangled_sites) + length(unentangled_sites) == natoms(crystal) "Invalid entangled unit specification."
    @assert isempty(intersect(entangled_sites, unentangled_sites)) "Invalid entangled unit specification." # Should be true by construction -- remove later

    # Prepare mapping dictionaries and initialize iteration variable.
    forward_map = Dict{Int64, Tuple{Int64, Int64}}()
    inverse_map = Dict{Tuple{Int64, Int64}, Tuple{Int64, Vec3}}()
    new_site_current = 1

    # Add sites that are *not* encapsulated in any unit as identical sites in
    # the new crystal.
    new_positions = []
    for site in unentangled_sites
        push!(new_positions, crystal.positions[site])
        new_pair = (new_site_current, 1)
        forward_map[site] = new_pair
        inverse_map[new_pair] = (site, Vec3(0, 0, 0))
        new_site_current += 1
    end

    # Assign entangled units to single site in new crystal and record mapping
    # information.
    for unit in units
        # Find new position by averaging location of entangled positions. 
        old_positions = [crystal.positions[i] for i in unit]
        new_position = sum(old_positions) / length(old_positions)
        push!(new_positions, new_position)

        # Record forward and inverse mapping information.
        for (j, site) in enumerate(unit)
            new_pair = (new_site_current, j)
            offset = crystal.positions[site] - new_position

            forward_map[site] = new_pair
            inverse_map[new_pair] = (site, offset)
        end

        new_site_current += 1
    end

    nsites_new = new_site_current - 1

    # Sorted list version of forward map
    forward_keys = collect(keys(forward_map))
    forward_vals = collect(values(forward_map))
    idcs = sortperm(forward_keys)
    forward_list = [forward_vals[i] for i in idcs]

    # Sorted list version of inverse map
    inverse_keys = collect(keys(inverse_map))
    inverse_vals = collect(values(inverse_map))
    idcs = sortperm(inverse_keys)
    inverse_keys = inverse_keys[idcs]
    inverse_vals = inverse_vals[idcs]

    inverse_list = [Tuple{Int64, Vec3}[] for _ in 1:nsites_new]
    for (n, key) in enumerate(inverse_keys)
        new_site, _ = key  # `key` is a tuple (global_site, local_site_in_unit)
        push!(inverse_list[new_site], inverse_vals[n])
    end

    # Generate a new contracted crystal and information to invert contraction.
    # NB: Perhaps set space group to 1 to allow appropriate interactions.
    # NB: Also add option to override to unconventional cell warnings and errors.
    new_crystal = Crystal(crystal.latvecs, new_positions)
    contraction_info = CrystalContractionInfo(forward_list, inverse_list)

    return new_crystal, contraction_info
end

function expand_crystal(contracted_crystal, contraction_info)
    (; forward, inverse) = contraction_info
    contracted_positions = contracted_crystal.positions
    nsites_expanded = length(forward) 
    expanded_positions = [Vec3(0, 0, 0) for _ in 1:nsites_expanded]
    for (contracted_idx, contracted_sites) in enumerate(inverse)
        for (original_site, offset) in contracted_sites
            expanded_positions[original_site] = contracted_positions[contracted_idx] + offset
        end
    end
    Crystal(contracted_crystal.latvecs, expanded_positions)
end


# unit_list should be an array of tuples
function contract_system(sys, units)
    # Construct contracted crystal
    contracted_cryst, contraction_info = contract_crystal(sys.crystal, units)

    # Iterate through interactions_union and construct new interactions

    
end