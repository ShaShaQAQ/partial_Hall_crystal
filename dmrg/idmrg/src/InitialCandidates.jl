struct InitialCandidate
    id::String
    kind::Symbol
    cdw_translation::Int
    dopant_sites::Vector{Int}
    occupied_sites::Vector{Int}
end

function Base.:(==)(first::InitialCandidate, second::InitialCandidate)
    return first.id == second.id &&
        first.kind == second.kind &&
        first.cdw_translation == second.cdw_translation &&
        first.dopant_sites == second.dopant_sites &&
        first.occupied_sites == second.occupied_sites
end

function paper_cdw_color(c::InfiniteCylinderConfig, site::Int)
    c.geometry == :paper_straight || throw(
        ArgumentError("paper CDW colors require paper_straight geometry")
    )
    triangular_x, triangular_y = paper_triangular_coordinates(c, site)
    return mod(triangular_x - triangular_y, 3)
end

function _translate_paper_site_y(
    c::InfiniteCylinderConfig,
    site::Int,
    translation::Int,
)
    x, y, orbital = paper_coordinates(c, site)
    return paper_site(c, x, y + translation, orbital)
end

function _reference_dopant_pair(
    c::InfiniteCylinderConfig,
    cdw_translation::Int,
)
    empty_sites = Int[]
    for x in 0:(c.x_period - 1), orbital in (:A, :B)
        site = paper_site(c, x, 0, orbital)
        paper_cdw_color(c, site) == cdw_translation || push!(empty_sites, site)
    end
    sort!(empty_sites)
    isempty(empty_sites) && error("CDW background has no empty dopant site")
    first_site = first(empty_sites)
    first_x, _, _ = paper_coordinates(c, first_site)
    second_index = findfirst(empty_sites) do site
        site_x, _, _ = paper_coordinates(c, site)
        site_x != first_x
    end
    isnothing(second_index) && error(
        "CDW background has no dopant pair in distinct original cells"
    )
    return first_site, empty_sites[second_index]
end

function _validate_fig2_candidate_config(c::InfiniteCylinderConfig)
    c.geometry == :paper_straight || throw(
        ArgumentError("Fig. 2 candidates require paper_straight geometry")
    )
    c.Ny == 6 || throw(ArgumentError("Fig. 2 candidates require Ny=6"))
    c.x_period == 3 || throw(
        ArgumentError("Fig. 2 candidates require x_period=3")
    )
    (c.filling_num, c.filling_den) == (7, 9) || throw(
        ArgumentError("Fig. 2 candidates require filling 7/9")
    )
    sites_per_cell(c) == 36 || error("Fig. 2 candidate cell must have 36 sites")
    particles_per_cell(c) == 14 || error(
        "Fig. 2 candidate cell must have 14 particles"
    )
    return nothing
end

function fig2_initial_candidates(c::InfiniteCylinderConfig)
    _validate_fig2_candidate_config(c)
    candidates = InitialCandidate[]
    for cdw_translation in 0:2
        background = sort([
            site for site in 1:sites_per_cell(c) if
            paper_cdw_color(c, site) == cdw_translation
        ])
        length(background) == 12 || error(
            "period-three CDW background must contain 12 particles"
        )
        reference_pair = _reference_dopant_pair(c, cdw_translation)
        for transverse_translation in 0:2
            dopants = sort([
                _translate_paper_site_y(c, site, transverse_translation) for
                site in reference_pair
            ])
            isempty(intersect(background, dopants)) || error(
                "dopant overlaps the occupied CDW background"
            )
            occupied = sort([background; dopants])
            length(unique(occupied)) == particles_per_cell(c) || error(
                "Fig. 2 candidate does not have 14 distinct occupied sites"
            )
            push!(
                candidates,
                InitialCandidate(
                    "cdw_t$(cdw_translation)_dopant_y$(transverse_translation)",
                    :cdw_doped,
                    cdw_translation,
                    dopants,
                    occupied,
                ),
            )
        end
    end

    push!(
        candidates,
        InitialCandidate(
            "uniform_control",
            :uniform_control,
            -1,
            Int[],
            default_occupied_sites(c),
        ),
    )
    length(unique(candidate.occupied_sites for candidate in candidates)) ==
        length(candidates) || error("Fig. 2 candidate patterns are not unique")
    return candidates
end
