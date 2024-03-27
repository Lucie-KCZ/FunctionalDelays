# Functional diversity loss and taxonomic delays of European freshwater fish and North American Breeding birds
**Authors** Lucie Kuczynski, Ana Maria Bastidas Urrutia, Helmut Hillebrand

## Abstract
Biodiversity is temporally dynamic, reflecting historical environmental conditions and influencing ecosystem stability. Colonisation and extinction dynamics frequently exhibit asynchronous patterns, resulting in net imbalances and thus to long-lasting richness trends. If these trends are not functionally random, functional Net Imbalances between Colonisations and Extinctions (fNICE) are likely to emerge. Using community time series data of European freshwater fish and North American breeding birds, we investigated how fNICE differs from its taxonomic equivalent (tNICE) to provide a comprehensive picture of biodiversity dynamics. Our findings reveal that taxonomic and functional delays are a prevalent feature, challenging the assumption of an immediate response to environmental changes. Taxonomic delays manifest as extinction debts and colonization credits, while functional delays indicate a shift in the balance between functional gains and losses over time. Moreover, we found that taxonomic and functional imbalances are not always directly correlated, although some specific patterns were found consistently for fish and birds. Early colonisations outpaced functional gains, indicating that although new species arrived earlier than the extinction of other species, the acquisition of new functional traits lagged. Although this may temporarily stabilize communities, as functional redundancy can mitigate loss of function via local extinctions, excessive redundancy can compromise biodiversity's capacity to respond to environmental variations, thereby undermining long-term resilience. In conclusion, understanding the intricate temporal dynamics of biodiversity responses is paramount for effective conservation practices. While short-term observations may suggest an equilibrium between diversity and the environmental conditions, our results underscore the importance of considering long-term dynamics and the interplay between species traits and changing environments. The metrics tNICE and fNICE are valuable tools for quantifying these temporal dynamics and unravelling their consequences for ecosystem stability. Incorporating these insights into conservation strategies can aid in proactively preserving biodiversity and safeguarding the integrity of ecosystems.

## Structure of the script
The `fNICE_script.R` script is composed of different sections:
- Importing the data and running the simulations (L23 to L127)
- Running both NICE metrics (L130 to L185)
- Describing the NICE trends (L188 to L296)

### Data import
The data are formatted as a list for which each element is a list containing i) the community matrix, ii) information about the community (e.g., sampling methods, total species richness, ...), and iii) the occurring species' traits. Mostly, data were processed beforehand so most of the selection criteria, especially for BBS, don't appear here. However, the main aspect is that one can see how to simulate data (for the neutral model) from the observed ones. The data were the same as used in [Kuczynski *et al.* (2023)](https://www.nature.com/articles/s41559-023-02078-w) (code available [here](https://github.com/Lucie-KCZ/NeutralDynamics/tree/main)).


### Notes
L49: `list.files('~/Dropbox/Neutral trends/analysis/functions island package/R/', full.names = T)` is a manual import of the package `island` which can be found here: [CRAN page for the island pkg](https://cran.r-project.org/web/packages/island/index.html).


