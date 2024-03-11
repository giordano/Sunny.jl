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

struct EntangledSystem
    sys               :: System
    crystal_origin    :: Crystal
    contraction_info  :: CrystalContractionInfo
    Ns_unit           :: Vector{Vector{Int64}}   # This eliminates need to carry original system around in many places
end

function EntangledSystem(sys, units)
    (; sys_entangled, contraction_info, Ns_unit) = entangle_system(sys, units)
    EntangledSystem(sys_entangled, orig_crystal(sys), contraction_info, Ns_unit)
end


################################################################################
# Aliases
################################################################################
randomize_spins!(esys::EntangledSystem; kwargs...) = randomize_spins!(esys.sys; kwargs...)
minimize_energy!(esys::EntangledSystem; kwargs...) = minimize_energy!(esys.sys; kwargs...)
plot_spins(esys::EntangledSystem; kwargs...) = plot_spins(esys.sys; kwargs...)