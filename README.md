# GTG-LBSPR_DomeShaped
Original GTG LB-SPR method <https://github.com/AdrianHordyk/GTG_LBSPR> modified to include 
* dome-shaped selectivity
* fixed logistic selection curves 
in the simulation and estimation routine.

We noted in Hommik *et al.* (2020) that 

> The original versions of GTG LB-SPR internally estimates selectivityat-
length. This approach is unresolved for dome-shaped selection, because
it remains difficult to distinguish whether fishing mortality of
larger individuals or the descending limb of the selectivity curve causes
the typical decrease in catch numbers with length.

and for dome-shaped curves

> ...only the ratio of fishing mortality to natural mortality (F/M) is estimated during
maximum likelihood optimisation.

That is, the dome-shaped selection curve must be specified *a priori* based on selectivity data or expert judgment. The code supports normal and lognormal dome-shaped curves. 
We also support multiple mesh sizes; however, the demonstration code GTGLBSPR_DomeShaped_SimulateEstimate assumes a single mesh size for simplicity.

## Code details

GTGLBSPR_Dome.R contains the GTG-LBSPR code including 
 * the per-recruit simulation (operating model) function **GTGDomeLBSPRSim(StockPars, FleetPars, SizeBins)**;
 * the top-level estimation model optimisation function **DoOptDome(StockPars, fixedFleetPars, LenDat, SizeBins, mod)** which requires length composition data LenDat. This function calls optimisation routine nlminb which seeks to minimise the difference between length composition data and expected (per-recruit theory) length composition.
 * the negative log-likelihood calculation function (multinomial likelihood) which calculates the goodness-of-fit of the expected per-recruit length composition (estimated from FM and selectivity-at-length parameters) to the actual length data.


GTGLBSPRDome_SimulateEstimate.Rmd
* demonstrates the dome-shaped selectivity estimation approach using simulated data.

## Resources

* **Publication** Hommik, K., Fitzgerald C. J., Kelly, F. K. and Shephard, S. Dome-shaped selectivity in LB-SPR: Length-Based assessment of data-limited inland fish stocks sampled with gillnets. <https://doi.org/10.1016/j.fishres.2020.105574>
