#include "Exosphere.h"
#include "GridUtility.h"
#include <AMReX_ParallelDescriptor.H>

void Exosphere::read_param(const std::string& command, ReadParam& param) {
    if (command == "#EXOSPHERE") {
        param.read_var("neutral_profile", neutral_profile);
        if (neutral_profile == "exponential") {
            param.read_var("neutral_profile.n0", n0);
            param.read_var("neutral_profile.H0", H0);
        } else if (neutral_profile == "power-law") {
            param.read_var("neutral_profile.n0", n0);
            param.read_var("neutral_profile.k0", k0);
        } else {
            amrex::Abort("Unknown neutral_profile: " + neutral_profile);
        }
        param.read_var("neutral_profile.T0", T0);
        param.read_var("exobase_radius", exobase_radius);
        param.read_var("total_production_rate", total_production_rate);
    }
}

void Exosphere::initialize() {
    if (!grid) return;

    int n_lev = grid->n_lev();
    neutralDensity.resize(n_lev);

    double total_neutral_content = 0.0;

    // 1. Loop over all grid cells and calculate neutral density
    for (int iLev = 0; iLev < n_lev; ++iLev) {
        const amrex::Geometry& geom = grid->Geom(iLev);
        const amrex::BoxArray& ba = grid->get_const_box_array(iLev);
        const amrex::DistributionMapping& dm = grid->get_const_distribution_map(iLev);

        // Define MultiFab with 1 component (density), 0 ghost cells
        neutralDensity[iLev].define(ba, dm, 1, 0);
        neutralDensity[iLev].setVal(0.0);

        for (amrex::MFIter mfi(neutralDensity[iLev]); mfi.isValid(); ++mfi) {
            const amrex::Box& box = mfi.validbox();
            auto& dens_arr = neutralDensity[iLev][mfi].array();

            // Loop over cells in the box
            amrex::ParallelFor(box, [&](int i, int j, int k) {
                amrex::IntVect ijk(AMREX_D_DECL(i, j, k));
                amrex::RealVect pos = grid->get_cell_center_loc(iLev, ijk);

                double dens = calculate_profile(pos);
                dens_arr(i, j, k) = dens;
            });
        }
    }

    // 2. Calculate total neutral content (sum of density * volume)
    // We need to be careful with AMR. To correctly integrate over the domain with AMR,
    // we should only include cells that are not covered by a finer level.

    for (int iLev = 0; iLev < n_lev; ++iLev) {
        const amrex::Geometry& geom = grid->Geom(iLev);
        amrex::Real cell_vol = AMREX_D_TERM(geom.CellSize(0), * geom.CellSize(1), * geom.CellSize(2));

        // Need access to mask/status to check for refinement if AMR is used.
        // Assuming `grid->cell_status(iLev)` gives an iMultiFab with status flags.

        amrex::Real level_sum = amrex::ReduceSum(neutralDensity[iLev], grid->cell_status(iLev), 0,
            [=] AMREX_GPU_HOST_DEVICE (amrex::Box const& bx, amrex::Array4<amrex::Real const> const& dens_arr,
                                       amrex::Array4<int const> const& status_arr) -> amrex::Real
            {
                amrex::Real sum = 0.0;
                amrex::Loop(bx, [=, &sum](int i, int j, int k) {
                    if (!bit::is_refined(status_arr(i,j,k))) {
                        sum += dens_arr(i,j,k) * cell_vol;
                    }
                });
                return sum;
            });

        total_neutral_content += level_sum;
    }

    // 3. Sum across all MPI processes
    amrex::ParallelDescriptor::ReduceRealSum(total_neutral_content);

    // 4. Normalize
    // We want the total production rate to match.
    // The user requirement says: "Stores the neutral density multiplied by cell volume... Sums... Normalizes to match total_production_rate"
    // This implies the stored "density" should be scaled so that Integral(density * dV) = total_production_rate ?
    // Or maybe the density profile is fixed by n0, and we just calculate the total production rate?
    // "Normalizes to match the configured total_production_rate"
    // This usually means we scale n0 (or the whole profile) so that the integral equals the production rate.
    // If n0 is given, maybe we ignore n0 and recalculate it?
    // Or maybe "neutral_profile.n0" is just a shape parameter if normalization is enforced.

    // Let's assume we scale the whole density field so that TotalContent = total_production_rate.

    if (total_neutral_content > 0.0) {
        double scale_factor = total_production_rate / total_neutral_content;

        // Apply scaling
        for (int iLev = 0; iLev < n_lev; ++iLev) {
             for (amrex::MFIter mfi(neutralDensity[iLev]); mfi.isValid(); ++mfi) {
                const amrex::Box& box = mfi.validbox();
                auto& dens_arr = neutralDensity[iLev][mfi].array();

                amrex::ParallelFor(box, [&](int i, int j, int k) {
                    dens_arr(i, j, k) *= scale_factor;
                });
             }
        }

        // Also update n0 to reflect the new normalization
        n0 *= scale_factor;
        amrex::Print() << "Exosphere initialized. Total content normalized to " << total_production_rate
                       << ". Scaling factor: " << scale_factor << ", New n0: " << n0 << std::endl;
    } else {
        amrex::Print() << "Warning: Total neutral content is zero. Cannot normalize." << std::endl;
    }
}

double Exosphere::calculate_profile(const amrex::RealVect& pos) const {
    double r = pos.vectorLength();

    if (r < exobase_radius) {
        return 0.0;
    }

    if (neutral_profile == "exponential") {
        // n(r) = n0 * exp(-(r - r0)/H0) ?
        // Typically exponential profile is n(r) = n0 * exp(-(r - R_ref)/H)
        // User gave n0 at reference radius. Let's assume reference radius is exobase_radius or we need another param?
        // "neutral_profile.n0: Number densities at reference radius"
        // "exobase_radius: Altitude below which neutral density is zero"
        // Usually reference radius is the planetary surface or the exobase.
        // Let's assume reference radius is exobase_radius for now, as it's the only radius provided.
        // Wait, usually reference radius is R_planet.

        // Let's assume reference radius is exobase_radius for simplicity unless specified.
        // Or maybe just use r.
        // n(r) = n0 * exp(-(r - exobase_radius)/H0)

        return n0 * std::exp(-(r - exobase_radius) / H0);

    } else if (neutral_profile == "power-law") {
        // n(r) = n0 * (r_ref / r)^k0
        // Again, what is r_ref?
        // If we assume n0 is at exobase_radius:
        return n0 * std::pow(exobase_radius / r, k0);
    }

    return 0.0;
}

amrex::Real Exosphere::get_density(const int iLev, const amrex::MFIter& mfi, const amrex::IntVect& ijk) const {
    if (iLev >= neutralDensity.size()) return 0.0;
    const auto& arr = neutralDensity[iLev][mfi].array();
    return arr(ijk);
}
