# spread

An R package that contains different infectious disease spread models.

Currently we have implemented [commuter model](https://folkehelseinstituttet.github.io/spread/articles/commuter_model.html). This model is a stochastic SEIIaR (susceptible, exposed, infectious, infectious asymptomatic, recovered) metapopulation model that including commuting. Each location has a local infection system, while the locations are connected by people who commute each day. The model differentiates between day and night. During the day you can infect/be infected in the location where you work, while during the night you can infect/be infected in the location where you live. It is the same commuters who travel back and forth each day. At the start of a day, all commuters are sent to their work location, where they mix for 12 hours. The commuters are then sent to their respective home locations, where they mix for 12 hours. The model is based upon a published model.

## fhiverse

The `fhiverse` is a set of R packages developed by the Norwegian Institute of Public Health to help solve problems that frequently occur when performing infectious disease surveillance.

If you want to install the dev versions (or access packages that haven't been released on CRAN), run `usethis::edit_r_profile()` to edit your `.Rprofile`. Then write in:

```
options(repos=structure(c(
  FHI="https://folkehelseinstituttet.github.io/drat/",
  CRAN="https://cran.rstudio.com"
)))
```

Save the file and restart R. This will allow you to install `fhiverse` packages from the FHI registry.

Current `fhiverse` packages are:

| Name    	| Info                                                             	|
|---------	|------------------------------------------------------------------	|
| [org](https://folkehelseinstituttet.github.io/org)         	| A system to help you organize projects.  |
| [plnr](https://folkehelseinstituttet.github.io/plnr)    	  | A system to help you plan analyses.  |
| [attrib](https://folkehelseinstituttet.github.io/attrib)  	| Calculating attributable mortalities and incident risk ratios.  |
| [spread](https://folkehelseinstituttet.github.io/spread)  	| Different infectious disease spread models.  |
| [fhidata](https://folkehelseinstituttet.github.io/fhidata) 	| Preformatted structural data for Norway.  |
| [fhimaps](https://folkehelseinstituttet.github.io/fhimaps) 	| Preformatted maps of Norway that generally don't need geolibraries.  |
| [fhiplot](https://folkehelseinstituttet.github.io/fhiplot) 	| Helpful functions for creating outputs in the style used by FHI.  |
