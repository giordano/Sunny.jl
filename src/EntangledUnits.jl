# Data for mapping one site inside a unit back to the site of the original
# system.
struct InverseData
    site   :: Int64
    offset :: Vec3
end

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
    inverse :: Vector{Vector{InverseData}}   # List ordered according to contracted crystal sites. Each element is itself a list containing original crystal site indices and corresponding offset information 
end

struct EntangledSystem{N}
    sys         :: System{N}
    sys_origin  :: System
    ci          :: CrystalContractionInfo
    Ns_internal :: Vector{Vector{Int64}}
end

function EntangledSystem(sys::System{N}, units) where N
    sys_contracted, contraction_info = contract_system(sys, units)
    Ns_internal = Ns_in_units(sys, contraction_info)
    EntangledSystem(sys_contracted, sys, contraction_info, Ns_internal)
end


# Takes a crystal and a list of integer tuples. The tuples indicate which sites
# in the original crystal are to be grouped, i.e., contracted into a single site
# of a new crystal.
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
    inverse_map = Dict{Tuple{Int64, Int64}, InverseData}()
    new_site_current = 1

    # Add sites that are *not* encapsulated in any unit as individual sites in
    # the new crystal with no associated displacement.
    new_positions = []
    for site in unentangled_sites
        push!(new_positions, crystal.positions[site])
        new_pair = (new_site_current, 1)
        forward_map[site] = new_pair
        inverse_map[new_pair] = InverseData(site, Vec3(0, 0, 0))
        new_site_current += 1
    end

    # Assign entangled units to single site in new crystal and record mapping
    # information.
    for unit in units
        # Find new position by averaging location of entangled positions. 
        old_positions = [crystal.positions[i] for i in unit]
        new_position = sum(old_positions) / length(old_positions)
        push!(new_positions, new_position)

        # Record forward and inverse mapping information, including
        # the displacement data from the new unit position.
        for (j, site) in enumerate(unit)
            new_pair = (new_site_current, j)
            offset = crystal.positions[site] - new_position

            forward_map[site] = new_pair
            inverse_map[new_pair] = InverseData(site, offset)
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

    inverse_list = [InverseData[] for _ in 1:nsites_new]
    for (n, key) in enumerate(inverse_keys)
        new_site, _ = key  # `key` is a tuple (global_site, local_site_in_unit)
        push!(inverse_list[new_site], inverse_vals[n])
    end

    # Generate a new contracted crystal and information to invert contraction.
    # Space group must be set to 1 to allow "invalid" anisotropies -- these are
    # generated by PairCouplings that become onsite couplings. NB: Add option to
    # override to unconventional cell warnings and errors? (Probably yes)
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
    for (contracted_idx, original_site_data) in enumerate(inverse)
        for (; site, offset) in original_site_data
            expanded_positions[site] = contracted_positions[contracted_idx] + offset
        end
    end
    Crystal(contracted_crystal.latvecs, expanded_positions)
end

# Returns a list of length equal to the number of "units" in the a contracted
# crystal. Each list element is itself a list of integers, each of which
# corresponds to the N of the corresponding site of the original system. The
# order is consistent with that given by the `inverse` field of a
# `CyrstalContractionInfo`.
function Ns_in_units(sys_original, contraction_info)
    Ns = [Int64[] for _ in 1:length(contraction_info.inverse)] 
    for (n, contracted_sites) in enumerate(contraction_info.inverse)
        for (; site) in contracted_sites
            push!(Ns[n], sys_original.Ns[site])
        end
    end
    Ns
end

# Given a local operator, A, that lives within an entangled unit on local site
# i, construct I ⊗ … ⊗ I ⊗ A ⊗ I ⊗ … ⊗ I, where A is in the i position.
# TODO: Incorporate this into `to_product_space` to unify code base.
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
sites_in_unit(contraction_info, i) = [inverse_data.site for inverse_data in contraction_info.inverse[i]] 

# List of all pair-wise bonds in a unit. The resulting bonds are to be
# interpreted in terms of the original crystal.
function bonds_in_unit(contraction_info, i)
    sites = sites_in_unit(contraction_info, i)
    nsites = length(sites)
    bonds = Bond[]
    for i in 1:nsites, j in i+1:nsites
        push!(bonds, Bond(i, j, [0, 0, 0]))
    end
    return bonds
end

# Checks whether a bond given in terms of the original crystal lies inside a
# single unit of the contracted system.
function bond_is_in_unit(bond::Bond, ci::CrystalContractionInfo)
    (ci.forward[bond.i][1] == ci.forward[bond.j][1]) && (bond.n == [0, 0, 0])
end

# Converts what was a pair coupling between two different sites in the original system into a single
# on-bond operator (an onsite operator in terms of the "units".)
function accum_pair_coupling_into_bond_operator_in_unit!(op, pc, sys, contracted_site, contraction_info)
    (; bond, scalar, bilin, biquad, general, isculled) = pc
    isculled && return

    Ns_all = Ns_in_units(sys, contraction_info)
    Ns_unit = Ns_all[contracted_site]
    I_unit = I(prod(Ns_unit))

    # Collect local Hilbert space information
    i, j = bond.i, bond.j
    @assert contraction_info.forward[i][1] == contracted_site "Sanity check -- remove later"
    @assert contraction_info.forward[j][1] == contracted_site "Sanity check -- remove later"
    i_unit = contraction_info.forward[i][2]
    j_unit = contraction_info.forward[j][2]
    Ni = sys.Ns[1, 1, 1, i] 
    Nj = sys.Ns[1, 1, 1, j] 

    # Add scalar part
    op .+= scalar*I_unit

    # Add bilinear part
    J = bilin isa Float64 ? bilin*I(3) : bilin
    Si = [local_op_to_unit_op(Sa, i_unit, Ns_unit) for Sa in spin_matrices((Ni-1)/2)]
    Sj = [local_op_to_unit_op(Sb, j_unit, Ns_unit) for Sb in spin_matrices((Nj-1)/2)]
    op .+= Si' * J * Sj

    # Add biquadratic part
    K = biquad isa Float64 ? diagm(biquad * Sunny.scalar_biquad_metric) : biquad
    Oi = [local_op_to_unit_op(Oa, i_unit, Ns_unit) for Oa in stevens_matrices_of_dim(2; N=Ni)]
    Oj = [local_op_to_unit_op(Ob, j_unit, Ns_unit) for Ob in stevens_matrices_of_dim(2; N=Nj)]
    op .+= Oi' * K * Oj

    # Add general part
    for (A, B) in general.data
        op .+= local_op_to_unit_op(A, i_unit, Ns_unit) * local_op_to_unit_op(B, j_unit, Ns_unit)
    end
end


# Converts a pair coupling from the original system into a pair coupling between
# units in the new system.
function pair_coupling_into_bond_operator_between_units(pc, sys, contraction_info)
    (; bond, scalar, bilin, biquad, general) = pc
    (; i, j, n) = bond
    unit1, unitsub1 = contraction_info.forward[i]
    unit2, unitsub2 = contraction_info.forward[j]

    Ns_local = Ns_in_units(sys, contraction_info)
    Ns_contracted = map(Ns -> prod(Ns), Ns_local)
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
        bond_operator += kron( local_op_to_unit_op(A, unitsub1, Ns1), I(Ns_contracted[unit1]) ) * 
                            kron( I(Ns_contracted[unit2]), local_op_to_unit_op(B, unitsub2, Ns2) )
    end

    return (; newbond, bond_operator)
end

function contract_system(::System{0}, _) 
    error("Cannot contract a dipole system.")
end

function contract_system(sys::System{M}, units) where M

    # Construct contracted crystal
    contracted_crystal, contraction_info = contract_crystal(sys.crystal, units)

    # Determine Ns for local Hilbert spaces (all must be equal). (TODO: Determine if alternative behavior preferable in mixed case.)
    Ns_local = Ns_in_units(sys, contraction_info)
    Ns_contracted = map(Ns -> prod(Ns), Ns_local)
    @assert allequal(Ns_contracted) "After contraction, the dimensions of the local Hilbert spaces on each generalized site must all be equal."

    # Construct empty contracted system
    dims = size(sys.dipoles)[1:3]
    spin_infos = [SpinInfo(i; S=(N-1)/2, g=1.0) for (i, N) in enumerate(Ns_contracted)]  # TODO: Decisions about g-factor 
    sys_contracted = System(contracted_crystal, dims, spin_infos, :SUN)

    # For each contracted site, scan original interactions and reconstruct as necessary.
    new_pair_data = Tuple{Bond, Matrix{ComplexF64}}[]
    for (contracted_site, N) in zip(1:natoms(contracted_crystal), Ns_contracted)


        ## Onsite portion of interaction 
        ## TODO: Add Zeeman term
        relevant_sites = sites_in_unit(contraction_info, contracted_site)
        original_interactions = sys.interactions_union[relevant_sites] 
        onsite_contracted = zeros(ComplexF64, N, N)
        for (site, interaction) in zip(relevant_sites, original_interactions)
            onsite_original = interaction.onsite
            unit_index = contraction_info.forward[site][2]
            onsite_contracted += local_op_to_unit_op(onsite_original, unit_index, Ns_local[contracted_site])
        end

        # Sort all PairCouplings in couplings that will be within a unit and couplings that will be between units
        pcs_intra = PairCoupling[] 
        pcs_inter = PairCoupling[]
        for interaction in original_interactions, pc in interaction.pair
            (; bond) = pc
            if bond_is_in_unit(bond, contraction_info)
                push!(pcs_intra, pc)
            else
                push!(pcs_inter, pc)
            end
        end

        # Convert intra-unit PairCouplings to onsite couplings
        for pc in pcs_intra
            accum_pair_coupling_into_bond_operator_in_unit!(onsite_contracted, pc, sys, contracted_site, contraction_info)
        end
        set_onsite_coupling!(sys_contracted, onsite_contracted, contracted_site)

        ## Convert inter-unit PairCouplings into new pair couplings
        for pc in pcs_inter
            (; newbond, bond_operator) = pair_coupling_into_bond_operator_between_units(pc, sys, contraction_info)
            push!(new_pair_data, (newbond, bond_operator))
        end
    end

    # Now have list of bonds and bond operators. First we must find individual
    # exemplars of each symmetry class of bonds in terms of the *new* crystal.
    all_bonds_with_interactions = [data[1] for data in new_pair_data]
    exemplars = Bond[]
    while length(all_bonds_with_interactions) > 0
        exemplar = all_bonds_with_interactions[1]
        all_bonds_with_interactions = filter(all_bonds_with_interactions) do bond
            !is_related_by_symmetry(contracted_crystal, bond, exemplar)
        end
        push!(exemplars, exemplar)
    end

    # We collected all bond_operators associated with a particular exemplar, sum
    # them, and set the interaction
    for bond in exemplars
        relevant_interactions = filter(data -> data[1] == bond, new_pair_data)
        bond_operator = sum(data[2] for data in relevant_interactions)
        set_pair_coupling!(sys_contracted, bond_operator, bond)
    end

    return (; sys_contracted, contraction_info)
end





# Make optimized version, also generalize to work on list of observables
function expand_contracted_system!(sys, sys_contracted, contraction_info)
    expectation(op, Z) = real(Z' * op * Z)

    for contracted_site in Sunny.eachsite(sys_contracted)
        i, j, k, n = contracted_site.I
        Ns_local = Ns_in_units(sys, contraction_info)
        Ns = Ns_local[n]

        # This iteration will be slow because it's on the last index...
        for (m, (; site)) in enumerate(contraction_info.inverse[n])
            S = spin_matrices((Ns[m]-1)/2)
            Sx = local_op_to_unit_op(S[1], m, Ns)
            Sy = local_op_to_unit_op(S[2], m, Ns)
            Sz = local_op_to_unit_op(S[3], m, Ns)
            dipole = Sunny.Vec3(
                expectation(Sx, sys_contracted.coherents[contracted_site]),
                expectation(Sy, sys_contracted.coherents[contracted_site]),
                expectation(Sz, sys_contracted.coherents[contracted_site]),
            )
            sys.dipoles[i, j, k, site] = dipole
        end
    end
end