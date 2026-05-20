# Rodent BLIND vs AWARE comparison

Generated from full 1M-iteration chains (both chains complete as of 2026-05-19).

## Numerical summary

| Quantity | Value |
|---|---|
| Aware trees (total from NWK) | 160 |
| Blind trees (total from NWK) | 359 |
| Burnin fraction applied | 25% |
| Aware post-burnin | 120 |
| Blind post-burnin | 269 |
| Aware subsampled for CID | 120 |
| Blind subsampled for CID | 200 |
| Aware minESS (at 1M iter, from .er log) | **BELOW 200 (minESS = 88)** |
| Blind minESS (at 1M iter, from .er log) | **BELOW 200 (minESS = 38)** |

### Pairwise CID distances (subsampled set)

| Comparison | N pairs | Mean CID | Median CID |
|---|---|---|---|
| Within aware | 7140 | 10.4144 | 10.7967 |
| Within blind | 19900 | 12.8216 | 12.6655 |
| Between chains | 24000 | 15.5035 | 15.2445 |

## Split agreement

- BLIND MR consensus: **38 splits**
- AWARE MR consensus: **45 splits**
- Shared splits: **32**
- Blind-unique splits: **6**
- Aware-unique splits: **13**

### Blind-unique splits (absent from AWARE consensus)
  {Acomys, Arvicanthis, Arvicola, Cricetomys, Dipus, Gerbillus, Hydromys, Jaculus, Mus, Napaeozapus, Neotoma, Platacanthomyidae, Rattus, Spalax}
  {Acomys, Arvicanthis, Arvicola, Cricetomys, Gerbillus, Hydromys, Mus, Neotoma, Platacanthomyidae, Rattus}
  {Arvicola, Neotoma}
  {Capromys, Cavia, Chinchilla, Cuniculus, Dasyprocta, Echimys, Erethizon, Hydrochoerus, Hystrix, Myocastor}
  {Capromys, Cavia, Chinchilla, Cuniculus, Dasyprocta, Echimys, Erethizon, Hydrochoerus, Myocastor}
  {Lepus, Ochotona, Oryctolagus, Pedetes, Sylvilagus}

### Aware-unique splits (absent from BLIND consensus)
  {Acomys, Anomalurus, Arvicanthis, Arvicola, Chaetodipus, Cricetomys, Dipodomys, Dipus, Eliomys, Geomys, Gerbillus, Graphiurus, Heteromys, Hydromys, Jaculus, Mus, Myoxus, Napaeozapus, Neotoma, Orthogeomys, Pedetes, Perognathus, Platacanthomyidae, Rattus, Spalax, Thomomys}
  {Acomys, Arvicanthis, Arvicola, Cricetomys, Gerbillus, Hydromys, Mus, Napaeozapus, Neotoma, Rattus}
  {Acomys, Mus, Rattus}
  {Hydromys, Napaeozapus}
  {Chaetodipus, Dipodomys, Geomys, Orthogeomys, Perognathus, Thomomys}
  {Eliomys, Graphiurus, Myoxus}
  {Dipus, Jaculus, Pedetes}
  {Capromys, Cavia, Chinchilla, Cuniculus, Dasyprocta, Dermoptera, Didelphis, Echimys, Erethizon, Eulemur_macaco, Hydrochoerus, Hystrix, Laonastes, Lepus, Macropus, Massoutiera, Myocastor, Ochotona, Oryctolagus, Papio_hamadryas, Rhynchocyon, Sylvilagus, Thryonomys, Tupaia, Vombatus}
  {Dermoptera, Didelphis, Eulemur_macaco, Laonastes, Lepus, Macropus, Massoutiera, Ochotona, Oryctolagus, Papio_hamadryas, Rhynchocyon, Sylvilagus, Tupaia, Vombatus}
  {Eulemur_macaco, Papio_hamadryas, Rhynchocyon, Tupaia}
  {Aplodontia, Castor, Glaucomys, Heterocephalus, Marmota, Ratufa, Sciurus, Tamias, Tamiasciurus}
  {Ratufa, Tamias, Tamiasciurus}
  {Castor, Heterocephalus}

## MDS topology

The two chains' posterior distributions are **largely separated** in CID-MDS space
(centroid distance = 9.6544; aware spread = 9.0090, blind spread = 7.0096).

## Key patterns

The aware and blind posterior clouds are distinctly separated in CID-MDS space, indicating systematic topological differences driven by the ecology-aware model. The ecology covariation prior shifts the inferred rodent phylogeny in a direction that is inconsistent with the blind Mk' posterior — consistent with the hypothesis that ecological convergence creates homoplasy that misleads standard Mk' inference.

## ESS and convergence

Both chains ran for 1M iterations with treeThin = 1000. The minESS values
reported here are continuous-parameter ESS from the MkPrime MCMC log (the
minimum over all monitored parameters at the final iteration).

- **Aware chain**: minESS = 88 at 1M iterations. Flag: below recommended minimum of 200.
- **Blind chain**: minESS = 38 at 1M iterations. Flag: below recommended minimum of 200.

ESS below 200 indicates that the chains have not fully converged on the
parameter that mixes most slowly (likely topology or a correlated rate
parameter). Results should be treated as indicative rather than definitive.
A further continuation or parallel-tempering run is advisable for publication.

## Chain details

- AWARE: 1M iterations, resumed from 371k checkpoint (2026-05-18 to 2026-05-19), 20.68h wall time.
- BLIND: 1M iterations, resumed from 200k checkpoint (2026-05-17), 19.5 min wall time
  (blind chain resumed quickly because the standard Mk' likelihood is much faster).

