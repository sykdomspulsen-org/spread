# spread 2020.3.19

typo in variable name

# spread 2020.3.18

Including asymptomatic incidence in asymmetric_mobility_se1e2iiar

# spread 2020.3.14

Initial commit of asymmetric_mobility_se1e2iiar

# spread 2020.3.5

Initial commit of asymmetric_mobility

# spread 2020.1.31

Included branching process functions from GunnarOyvindIsaksson.Ro@fhi.no

- branching_process
- fit_params_bp
- plot_quantiles_bp
- summarize_bp

# spread 2020.1.30

Included option for "aggregate_location" in commuter function. This will automatically aggregate the results over all locations, resulting in a smaller memory footprint. Useful when `simulations` is a large number.

# spread 2020.1.28

Included option for "simulations" (number of simulations to be repeated) in commuter function via a `foreach` loop, allowing for simulations to be run in parallel if desired (backend must be registered externally to package).

# spread 2020.1.27

Included option for "verbose" in commuter function.

Exported functions:

- commuter
- convert_blank_seiiar_with_vax

Exported datasets:

- norway_commuters_2017_b2020
- norway_seiiar_noinfected_2017_b2020
- norway_seiiar_oslo_2017_b2020
- norway_seiiar_measles_noinfected_2017_b2020
- norway_seiiar_measles_oslo_2017_b2020
- single_entity_fake_commuters_2017.rda (new)
- single_entity_seiiar_2017.rda (new)

# spread 2020.1.23

Upgraded to 2020 borders.

Exported functions:

- commuter
- convert_blank_seiiar_with_vax

Exported datasets (these are now tagged with "_b2020" to show they are for 2020 borders):

- norway_commuters_2017_b2020
- norway_seiiar_noinfected_2017_b2020
- norway_seiiar_oslo_2017_b2020
- norway_seiiar_measles_noinfected_2017_b2020
- norway_seiiar_measles_oslo_2017_b2020

# spread 2019.8.5

Submitted to CRAN.

Exported functions:

- commuter
- convert_blank_seiiar_with_vax

Exported datasets:

- norway_commuters_2017
- norway_seiiar_noinfected_2017
- norway_seiiar_oslo_2017
- norway_seiiar_measles_noinfected_2017
- norway_seiiar_measles_oslo_2017
