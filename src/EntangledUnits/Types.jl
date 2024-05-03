# Data for mapping one site inside a unit back to the site of the original
# system.
struct InverseData
    site   :: Int64  # Atom index of original, uncontracted crystal
    offset :: Vec3   # Position offset of original atom relative to center of unit
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

struct EntangledSystem
    sys               :: System                  # System containing entangled units
    sys_origin        :: System                  # Original "uncontracted" system
    contraction_info  :: CrystalContractionInfo  # Forward and inverse mapping data for sys <-> sys_origin
end

function EntangledSystem(sys, units)
    (; sys_entangled, contraction_info) = entangle_system(sys, units)
    sys_origin = clone_system(sys)
    esys = EntangledSystem(sys_entangled, sys_origin, contraction_info)
    set_expected_dipoles_of_entangled_system!(sys.dipoles, esys)
    return esys
end

function randomize_spins!(esys::EntangledSystem; kwargs...) 
    randomize_spins!(esys.sys; kwargs...)
    set_expected_dipoles_of_entangled_system!(esys.sys_origin.dipoles, esys)
end

function minimize_energy!(esys::EntangledSystem; kwargs...)
    optout = minimize_energy!(esys.sys; kwargs...)
    set_expected_dipoles_of_entangled_system!(esys.sys_origin.dipoles, esys)
    return optout
end

energy(esys::EntangledSystem; kwargs...) = energy(esys.sys; kwargs...)

set_dipole!(esys::EntangledSystem, dipole, site; kwargs...) = error("Setting dipoles of an EntangledSystem not well defined.") 

function set_coherent!(esys::EntangledSystem, coherent, site; kwargs...) 
    set_coherent!(esys.sys, coherent, site; kwargs...)
    set_expected_dipole_of_entangled_system!(esys.sys_origin.dipoles, esys, site)
end

eachsite(esys::EntangledSystem) = eachsite(esys.sys) # Not sure that we want this

# TODO: set_external_field!(esys, B)


function magnetic_moment(esys::EntangledSystem, site; kwargs...) 
    magnetic_moment(esys.sys_origin, site; kwargs...)
end

function plot_spins(esys::EntangledSystem; kwargs...)
    plot_spins(esys.sys_origin; kwargs...)
end