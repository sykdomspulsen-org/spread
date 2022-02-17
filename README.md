# spread <a href="https://docs.sykdomspulsen.no/spread"><img src="man/figures/logo.png" align="right" width="120" /></a>

## Overview 

[spread](https://docs.sykdomspulsen.no/spread) contains different infectious disease spread models.

Currently we have implemented [commuter model](https://docs.sykdomspulsen.no/spread/articles/commuter_model.html). This model is a stochastic SEIIaR (susceptible, exposed, infectious, infectious asymptomatic, recovered) metapopulation model that including commuting. Each location has a local infection system, while the locations are connected by people who commute each day. The model differentiates between day and night. During the day you can infect/be infected in the location where you work, while during the night you can infect/be infected in the location where you live. It is the same commuters who travel back and forth each day. At the start of a day, all commuters are sent to their work location, where they mix for 12 hours. The commuters are then sent to their respective home locations, where they mix for 12 hours. The model is based upon a published model.

Read the introduction vignette [here](https://docs.sykdomspulsen.no/spread/articles/commuter_model.html) or run `help(package="spread")`.

## splverse

<a href="https://docs.sykdomspulsen.no/packages"><img src="https://docs.sykdomspulsen.no/packages/splverse.png" align="right" width="120" /></a>

The [splverse](https://docs.sykdomspulsen.no/packages) is a set of R packages developed to help solve problems that frequently occur when performing infectious disease surveillance.

If you want to install the dev versions (or access packages that haven't been released on CRAN), run `usethis::edit_r_profile()` to edit your `.Rprofile`. 

Then write in:

```
options(
  repos = structure(c(
    SPLVERSE  = "https://docs.sykdomspulsen.no/drat/",
    CRAN      = "https://cran.rstudio.com"
  ))
)
```

Save the file and restart R.

You can now install [splverse](https://docs.sykdomspulsen.no/packages) packages from our [drat registry](https://docs.sykdomspulsen.no/drat).

```
install.packages("spread")
```

