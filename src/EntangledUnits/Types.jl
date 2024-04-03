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
    forward :: Vector{Tuple{Int64, Int64}}  # Original site index -> full unit index (contracted crystal site index and unit subindex)
    inverse :: Vector{Vector{InverseData}}  # List ordered according to contracted crystal sites. Each element is itself a list containing original crystal site indices and corresponding offset information 
end

# struct EntanglementData
#     contraction_info  :: CrystalContractionInfo
#     Ns_unit           :: Vector{Vector{Int64}}   # This eliminates need to carry original system around in many places
# end

mutable struct EntangledSystem
    const sys               :: System
    const sys_origin        :: System
    synced                  :: Bool
    const contraction_info  :: CrystalContractionInfo
    const Ns_unit           :: Vector{Vector{Int64}}   # This eliminates need to carry original system around in many places -- now that system is there, possibly eliminate
end

function EntangledSystem(sys, units)
    (; sys_entangled, contraction_info, Ns_unit) = entangle_system(sys, units)
    sys_origin = clone_system(sys)
    EntangledSystem(sys_entangled, sys_origin, false, contraction_info, Ns_unit)
end


################################################################################
# Aliasing and field forwarding
################################################################################

# Functions to access System fields of EntangledSystem
randomize_spins!(esys::EntangledSystem; kwargs...) = randomize_spins!(esys.sys; kwargs...)
minimize_energy!(esys::EntangledSystem; kwargs...) = minimize_energy!(esys.sys; kwargs...)
energy(esys::EntangledSystem; kwargs...) = energy(esys.sys; kwargs...)
set_coherent!(esys::EntangledSystem, coherent, site; kwargs...) = set_coherent!(esys.sys, coherent, site; kwargs...)
eachsite(esys::EntangledSystem) = eachsite(esys.sys) # Not sure that we want this

# Functions acting on original System of an EntangledSystem 
set_dipole!(esys::EntangledSystem, dipole, site; kwargs...) = error("Setting dipoles of an EntangledSystem not well defined.") # Could replicate behavior of normal SU(N) system
function magnetic_moment(esys::EntangledSystem, site; kwargs...) 
    magnetic_moment(esys.sys_origin, site; kwargs...)
end
function plot_spins(esys::EntangledSystem; kwargs...)
    sync_dipoles(esys)
    plot_spins(esys.sys_origin; kwargs...)
end

# Functions acting on both systems
function reshape_supercell(esys::EntangledSystem, shape)
    (; sys, sys_origin, contraction_info) = esys
    new_sys_origin = reshape_supercell(sys_origin, shape)

    new_crystal = sys_origin_new.crystal
    new_na = natoms(new_crystal)

    units = original_units(esys)
    new_units = []

    new_atoms = 1:new_na
    new_units = []
    while length(new_atoms) > 0
        # Pick any site from list of new sites
        new_atom = new_atoms[1]
        new_site = CartesianIndex(1, 1, 1, new_atom) # Just work with first unit cell

        # Find corresponding original atom number
        site = position_to_site(sys_origin, position(new_sys_origin, new_site))
        atom = site[4]

        # Find the unit to which this original atom belongs
        unit = findfirst(unit -> atom in unit, units)

        # Find positions of all atoms in the unit, find corresponding sites in reshape system, and define unit for reshaped system
        unit_positions = [position(sys_origin, CartesianIndex(1, 1, 1, atom)) for atom in unit]
        new_unit_sites = [position_to_site(new_sys_origin, position) for position in unit_positions]
        new_unit = Int64[]
        for new_site in new_unit_sites
            i, j, k, a = new_site
            if !(i == j == k == 1)
                error("Specified reshaping incomptable with specified entangled units. (Unit split between unit cells.)")
            end
            push!(new_unit, a)
        end
        push!(new_units, Tuple(new_unit))

        idcs = findall(atom -> atom in new_unit, new_atoms)
        deleteat!(new_atoms, idcs)
    end

    # Construct new ContractionInfo for reshaped system
    _, contraction_info = contract_crystal(new_sys_origin.crystal, new_units)

    # Reshape the entangled system as well
    new_sys = reshape_supercell(sys, shape)

    # new_esys = EntangledSystem()

end

function repeat_periodically(esys, counts)
    sync_dipoles(esys)
    sys_new = repeat_periodically(esys.sys, counts)
    sys_origin_new = repeat_periodically(esys.sys_origin, counts)
    return EntangledSystem(sys_new, sys_origin_new, true, esys.contraction_info, esys.Ns_unit)
end


# Forward field requests to internal system
function Base.getproperty(value::EntangledSystem, name::Symbol)
    if name == :coherents
        return getfield(value.sys, name)
    elseif name == :dipoles
        return getfield(value.sys_origin, name)
    end
    return getfield(value, name)
end
function Base.setproperty!(value::EntangledSystem, name::Symbol, x)
    if name == :coherents 
        return setfield!(value.sys, name, convert(fieldtype(System, name), x))
    elseif name == :dipoles
        error("Cannot set `dipoles` of EntangledSystem directly.")
    end
    return setfield!(value, name, convert(fieldtype(EntangledSystem, name), x))
end


# TODO
# set_external_field!(esys, B)