
#pragma once

#include "../../types.hh"
#include <cmath>
#include <cstddef>
#include <spec/utils.hh>
#include <yaml-cpp/yaml.h>

namespace hard::tasks::initial_data {

//
// An initial data for testing the diffusion solver
//
template<std::size_t Dim>
auto
radiation_diffusion(typename mesh<Dim>::template accessor<ro> m,
  field<double>::accessor<wo, na> mass_density_a,
  typename field<vec<Dim>>::template accessor<wo, na> momentum_density_a,
  field<double>::accessor<wo, na> total_energy_density_a,
  field<double>::accessor<wo, na> radiation_energy_density_a,
  single<double>::accessor<ro> gamma_a,
  single<double>::accessor<ro> particle_mass_a) {

  auto mass_density = m.template mdcolex<is::cells>(mass_density_a);
  auto momentum_density = m.template mdcolex<is::cells>(momentum_density_a);
  auto total_energy_density =
    m.template mdcolex<is::cells>(total_energy_density_a);
  auto radiation_energy_density =
    m.template mdcolex<is::cells>(radiation_energy_density_a);

  if constexpr(Dim == 1) {

    using spec::utils::sqr;

    // Parse input parameters
    YAML::Node config = YAML::LoadFile(opt::config.value());
    const double sigma = config["problem_parameters"]["sigma"].as<double>();

    for(auto i : m.template cells<ax::x, dm::quantities>()) {
      // Initialize hydro variables with unit values.
      // It does not matter, since matter-light coupling should be turned off in
      // this test.
      mass_density(i) = 1.0;
      momentum_density(i).x = 0.0;
      total_energy_density(i) = 1.0;

      // const auto x = m.template center<ax::x>(i);
      // radiation_energy_density(i) = 2.0 + sin(2.0 * M_PI * x);

      // D = 1, t_0 = 1.0
      const auto x_shifted = m.template center<ax::x>(i) - 0.5;
      radiation_energy_density(i) =
        exp(-0.5 * spec::utils::sqr(x_shifted) / spec::utils::sqr(sigma));
      radiation_energy_density(i) *= 1 / (std::sqrt(2.0 * M_PI) * sigma);
    }
  }

  else if constexpr(Dim == 2) {

    for(auto j : m.template cells<ax::y, dm::quantities>()) {
      for(auto i : m.template cells<ax::x, dm::quantities>()) {
        // Initialize hydro variables with unit values.
        // It does not matter, since matter-light coupling should be turned off
        // in this test.
        mass_density(i, j) = 1.0;
        momentum_density(i, j).x = 0.0;
        momentum_density(i, j).y = 0.0;
        total_energy_density(i, j) = 1.0;

        const auto x = m.template center<ax::x>(i);
        const auto y = m.template center<ax::y>(j);

        radiation_energy_density(i, j) =
          2.0 + sin(2.0 * M_PI * x) * sin(2.0 * M_PI * y);
      }
    }
  }
}

} // namespace hard::tasks::initial_data
