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
    forward :: Vector{Tuple{Int64, Int64}}          # Original site index -> full unit index (contracted crystal site index and unit subindex)
    inverse :: Vector{Vector{Tuple{Int64, Vec3}}}   # List ordered according to contracted crystal sites. Each element is itself a list containing original crystal site indices and corresponding offset information 
end

# 
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
    @assert isempty(intersect(entangled_sites, unentangled_sites)) "Invalid entangled unit specification." # Sanity check. Should be true by construction. Remove later

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
    # NB: Space group must be set to 1 to allow invalid anisotropies.
    # NB: Add option to override to unconventional cell warnings and errors? (Probably yes)
    new_crystal = Crystal(crystal.latvecs, new_positions, 1)
    contraction_info = CrystalContractionInfo(forward_list, inverse_list)

    return new_crystal, contraction_info
end

# Reconstruct original crystal from contracted Crystal and a CrystalContractionInfo
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

function contracted_Ns(sys, contracted_cryst, contraction_info)
    Ns = [Int64[] for _ in 1:natoms(contracted_cryst)] 
    for contracted_sites in contraction_info.inverse
        for (original_site, _) in contracted_sites
            push!(Ns[], sys.Ns[original_site])
        end
    end
    Ns
end

# Given a local operator, A, that lives within an entangled unit on local site
# i, construct I ⊗ … ⊗ I ⊗ A ⊗ I ⊗ … ⊗ I, where A is in the i position.
function local_op_to_unit_op(A, i, Ns)
    @assert size(A, 1) == Ns[i] "Given operator not consistent with dimension of local Hilbert space"
    nsites = length(Ns) # Number of sites in the unit.

    # If there is only one site in the unit, simply return the original operator.
    (nsites == 1) && (return A)

    # Otherwise generate the appropriate tensor product.
    ops = [unit_index == i ? A : I(Ns[unit_index]) for unit_index in 1:nsites]

    return reduce(kron, ops)
end


# Pull out original indices of sites in entangled unit
sites_in_unit(contraction_info, i) = [site[1] for site in contraction_info.inverse[i]] 

# List of all pair-wise bonds in a unit.
function bonds_in_unit(contraction_info, i)
    sites = sites_in_unit(contraction_info, i)
    nsites = length(sites)
    bonds = Bond[]
    for i in 1:nsites, j in i+1:nsites
        push!(bonds, Bond(i, j, [0, 0, 0]))
    end
    return bonds
end

function bond_is_in_unit(bond::Bond, ci::CrystalContractionInfo)
    (ci.forward[bond.i][1] == ci.forward[bond.j][1]) && (bond.n == [0, 0, 0])
end

function contract_system(sys, units)

    # Construct contracted crystal
    contracted_crystal, contraction_info = contract_crystal(sys.crystal, units)

    # Determine Ns for local Hilbert spaces (all must be equal). (TODO: Determine if alternative behavior preferable in mixed case.)
    Ns_local = contracted_Ns(sys, contracted_crystal, contraction_info)
    Ns_contracted = map(Ns -> prod(Ns), Ns_local)
    @assert allequal(Ns_contracted) "After contraction, the dimensions of the local Hilbert spaces on each generalized site must all be equal."

    # Construct empty contracted system
    dims = size(sys.dipoles)[1:3]
    spin_infos = [SpinInfo(i; S=(N-1)/2, g=1.0) for (i, N) in enumerate(Ns_contracted)]  # TODO: Decisions about g-factor 
    sys_contracted = System(contracted_crystal, dims, spin_infos, :SUN)

    # For each contracted site, scan original interactions and reconstruct as necessary.
    for (contracted_site, N) in zip(1:natoms(contracted_crystal), Ns_contracted)

        ## Onsite portion of interaction 
        relevant_sites = sites_in_unit(contraction_info, contracted_site)
        original_interactions = sys.interactions_union[relevant_sites] 
        onsite_contracted = zeros(ComplexF64, N, N)
        for (site, interaction) in zip(relevant_sites, original_interactions)
            onsite_original = interaction.onsite        # TODO: convert if Stevens expansion? I.e., Allow contracting dipole models?
            unit_index = contraction_info.forward[site][2]
            onsite_contracted += local_op_to_unit_op(onsite_original, unit_index, Ns_local[contracted_site])
        end

        # Sort all PairCouplings in couplings that will be within a unit and couplings that will be between units
        pcs_intra = PairCoupling[] 
        pcs_inter = PairCoupling[]
        for interaction in original_interactions
            for pc in interaction.pair
                (; isculled, bond) = pc
                if !isculled
                    if bond_is_in_unit(bond, contraction_info)
                        push!(pcs_intra, pc)
                    else
                        push!(pcs_inter, pc)
                    end
                end
            end
        end

        # Convert intra-unit PairCouplings to onsite couplings
        Ns = Ns_local[contracted_site]
        I_unit = reduce(kron, [I(N) for N in Ns])
        for pc in pcs_intra
            # Note -- can ignore `isculled` because `bonds_in_unit` returns unique bonds
            (; bond, scalar, bilin, biquad, general) = pc

            # Collect local Hilbert space information
            i, j = bond.i, bond.j
            @assert contraction_info.forward[i][1] == contracted_site "Sanity check -- remove later"
            @assert contraction_info.forward[j][1] == contracted_site "Sanity check -- remove later"
            i_unit = contraction_info.forward[i][2]
            j_unit = contraction_info.forward[j][2]
            Ni = sys.Ns[1, 1, 1, i] 
            Nj = sys.Ns[1, 1, 1, j] 

            # Add scalar part
            onsite_contracted += scalar*I_unit

            # Add bilinear part
            J = bilin isa Float64 ? bilin*I(3) : bilin
            Si = [local_op_to_unit_op(Sa, i_unit, Ns) for Sa in spin_matrices((Ni-1)/2)]
            Sj = [local_op_to_unit_op(Sb, j_unit, Ns) for Sb in spin_matrices((Nj-1)/2)]
            onsite_contracted += Si' * J * Sj

            # Add biquadratic part
            K = biquad isa Float64 ? diagm(biquad * Sunny.scalar_biquad_metric) : biquad
            Oi = [local_op_to_unit_op(Oa, i_unit, Ns) for Oa in stevens_matrices_of_dim(2; N=Ni)]
            Oj = [local_op_to_unit_op(Ob, j_unit, Ns) for Ob in stevens_matrices_of_dim(2; N=Nj)]
            onsite_contracted += Oi' * K * Oj

            # Add general part
            for (A, B) in general.data
                onsite_contracted += local_op_to_unit_op(A, i_unit, Ns) * local_op_to_unit_op(B, j_unit, Ns)
            end
        end
        set_onsite_coupling!(sys_contracted, onsite_contracted, contracted_site)
        

        ## Convert inter-unit PairCouplings into new pair couplings
        new_pair_data = Tuple{Bond, Matrix{ComplexF64}}[]
        for pc in pcs_inter
            (; bond, scalar, bilin, biquad, general) = pc
            (; i, j, n) = bond
            unit1, unitsub1 = contraction_info.forward[i]
            unit2, unitsub2 = contraction_info.forward[j]
            Ns1 = Ns_local[unit1]
            Ns2 = Ns_local[unit2]
            N1 = sys.Ns[1, 1, 1, i]
            N2 = sys.Ns[1, 1, 1, j]
            N = Ns_contracted[unit1] * Ns_contracted[unit2]
            newbond = Bond(unit1, unit2, n)

            bond_operator = zeros(ComplexF64, N, N)

            # Add scalar part
            bond_operator += scalar*I(N)

            # Add bilinear part
            J = bilin isa Float64 ? bilin*I(3) : bilin
            Si = [kron(local_op_to_unit_op(Sa, unitsub1, Ns1), I(Ns_contracted[unit1])) for Sa in spin_matrices((N1-1)/2)]
            Sj = [kron(I(Ns_contracted[unit2]), local_op_to_unit_op(Sa, unitsub2, Ns2)) for Sa in spin_matrices((N2-1)/2)]
            bond_operator += Si' * J * Sj

            # Add biquadratic part
            K = biquad isa Float64 ? diagm(biquad * Sunny.scalar_biquad_metric) : biquad
            Oi = [kron(local_op_to_unit_op(Oa, unitsub1, Ns1), I(Ns_contracted[unit1])) for Oa in stevens_matrices_of_dim(2; N=N1)]
            Oj = [kron(I(Ns_contracted[unit2]), local_op_to_unit_op(Ob, unitsub2, Ns2)) for Ob in stevens_matrices_of_dim(2; N=N2)]
            bond_operator += Oi' * K * Oj

            # Add general part
            for (A, B) in general.data
                bond_operator += kron(local_op_to_unit_op(A, unitsub1, Ns1), I(Ns_contracted[unit1])) * kron(I(Ns_contracted[unit2]), local_op_to_unit_op(B, unitsub2, Ns2))
            end

            push!(new_pair_data, (newbond, bond_operator))
        end

        # Consolidate interactions on same generalized bonds and make final pair Couplings
        unique_bonds = unique([data[1] for data in new_pair_data])
        for bond in unique_bonds
            bond_data = filter(data -> data[1] == bond, new_pair_data)
            bond_operator = sum(data[2] for data in bond_data)
            set_pair_coupling!(sys_contracted, bond_operator, bond)
        end
    end

    return sys_contracted, contraction_info
end