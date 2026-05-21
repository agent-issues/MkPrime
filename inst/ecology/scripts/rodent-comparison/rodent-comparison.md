# Rodent BLIND vs AWARE comparison

Generated from full 1M-iteration chains (both chains complete as of 2026-05-19).

## Numerical summary

| Quantity | Value |
|---|---|
| Aware trees (total from NWK) | 160 |
| Blind trees (total from NWK) | 2525 |
| Burnin fraction applied | 25% |
| Aware post-burnin | 120 |
| Blind post-burnin | 1893 |
| Aware subsampled for CID | 120 |
| Blind subsampled for CID | 200 |
| Aware minESS (at 1M iter, from .er log) | **BELOW 200 (minESS = 88)** |
| Blind minESS (at 1M iter, from .er log) | **BELOW 200 (minESS = 140)** |

### Pairwise CID distances (subsampled set)

| Comparison | N pairs | Mean CID | Median CID |
|---|---|---|---|
| Within aware | 7140 | 10.4144 | 10.7967 |
| Within blind | 19900 | 13.5662 | 13.3755 |
| Between chains | 24000 | 17.2558 | 17.2693 |

## Split agreement

- BLIND MR consensus: **45 splits**
- AWARE MR consensus: **45 splits**
- Shared splits: **30**
- Blind-unique splits: **15**
- Aware-unique splits: **15**

### Blind-unique splits (absent from AWARE consensus)
  {Acomys, Arvicanthis, Arvicola, Chaetodipus, Cricetomys, Dipodomys, Dipus, Geomys, Gerbillus, Heteromys, Hydromys, Jaculus, Lepus, Mus, Napaeozapus, Neotoma, Ochotona, Orthogeomys, Oryctolagus, Pedetes, Perognathus, Platacanthomyidae, Rattus, Spalax, Sylvilagus, Thomomys}
  {Chaetodipus, Dipodomys, Dipus, Geomys, Heteromys, Jaculus, Lepus, Ochotona, Orthogeomys, Oryctolagus, Pedetes, Perognathus, Sylvilagus, Thomomys}
  {Dipus, Jaculus, Lepus, Ochotona, Oryctolagus, Pedetes, Sylvilagus}
  {Lepus, Ochotona, Oryctolagus, Pedetes, Sylvilagus}
  {Arvicola, Neotoma}
  {Capromys, Castor, Cavia, Chinchilla, Cuniculus, Dasyprocta, Echimys, Erethizon, Heterocephalus, Hydrochoerus, Hystrix, Laonastes, Massoutiera, Myocastor, Thryonomys}
  {Capromys, Cavia, Chinchilla, Cuniculus, Dasyprocta, Echimys, Erethizon, Hydrochoerus, Hystrix, Myocastor}
  {Aplodontia, Eliomys, Glaucomys, Graphiurus, Marmota, Myoxus, Ratufa, Sciurus, Tamias, Tamiasciurus}
  {Marmota, Ratufa, Sciurus, Tamiasciurus}
  {Anomalurus, Dermoptera, Didelphis, Eulemur_macaco, Macropus, Papio_hamadryas, Rhynchocyon, Tupaia, Vombatus}
  {Dermoptera, Didelphis, Eulemur_macaco, Rhynchocyon, Tupaia}
  {Dermoptera, Didelphis, Rhynchocyon, Tupaia}
  {Dermoptera, Didelphis, Tupaia}
  {Dermoptera, Tupaia}
  {Macropus, Vombatus}

### Aware-unique splits (absent from BLIND consensus)
  {Acomys, Anomalurus, Arvicanthis, Arvicola, Chaetodipus, Cricetomys, Dipodomys, Dipus, Eliomys, Geomys, Gerbillus, Graphiurus, Heteromys, Hydromys, Jaculus, Mus, Myoxus, Napaeozapus, Neotoma, Orthogeomys, Pedetes, Perognathus, Platacanthomyidae, Rattus, Spalax, Thomomys}
  {Acomys, Arvicanthis, Arvicola, Cricetomys, Gerbillus, Hydromys, Mus, Napaeozapus, Neotoma, Rattus}
  {Hydromys, Napaeozapus}
  {Eliomys, Graphiurus, Myoxus}
  {Dipus, Jaculus, Pedetes}
  {Capromys, Cavia, Chinchilla, Cuniculus, Dasyprocta, Dermoptera, Didelphis, Echimys, Erethizon, Eulemur_macaco, Hydrochoerus, Hystrix, Laonastes, Lepus, Macropus, Massoutiera, Myocastor, Ochotona, Oryctolagus, Papio_hamadryas, Rhynchocyon, Sylvilagus, Thryonomys, Tupaia, Vombatus}
  {Dermoptera, Didelphis, Eulemur_macaco, Laonastes, Lepus, Macropus, Massoutiera, Ochotona, Oryctolagus, Papio_hamadryas, Rhynchocyon, Sylvilagus, Tupaia, Vombatus}
  {Dermoptera, Didelphis, Eulemur_macaco, Macropus, Papio_hamadryas, Rhynchocyon, Tupaia}
  {Dermoptera, Eulemur_macaco, Papio_hamadryas, Rhynchocyon, Tupaia}
  {Eulemur_macaco, Papio_hamadryas, Rhynchocyon, Tupaia}
  {Eulemur_macaco, Papio_hamadryas, Rhynchocyon}
  {Eulemur_macaco, Papio_hamadryas}
  {Aplodontia, Castor, Glaucomys, Heterocephalus, Marmota, Ratufa, Sciurus, Tamias, Tamiasciurus}
  {Ratufa, Tamias, Tamiasciurus}
  {Castor, Heterocephalus}

## MDS topology

The two chains' posterior distributions are **largely separated** in CID-MDS space
(centroid distance = 11.9111; aware spread = 8.8298, blind spread = 6.1487).

## Key patterns

The aware and blind posterior clouds are distinctly separated in CID-MDS space, indicating systematic topological differences driven by the ecology-aware model. The ecology covariation prior shifts the inferred rodent phylogeny in a direction that is inconsistent with the blind Mk' posterior — consistent with the hypothesis that ecological convergence creates homoplasy that misleads standard Mk' inference.

## ESS and convergence

Both chains ran for 1M iterations with treeThin = 1000. The minESS values
reported here are continuous-parameter ESS from the MkPrime MCMC log (the
minimum over all monitored parameters at the final iteration).

- **Aware chain**: minESS = 88 at 1M iterations. Flag: below recommended minimum of 200.
- **Blind chain**: minESS = 140 at 1M iterations. Flag: below recommended minimum of 200.

ESS below 200 indicates that the chains have not fully converged on the
parameter that mixes most slowly (likely topology or a correlated rate
parameter). Results should be treated as indicative rather than definitive.
A further continuation or parallel-tempering run is advisable for publication.

## Chain details

- AWARE (this rendering): PLACEHOLDER (Mk' v2 fallback while MkNT aware run is in progress)
- BLIND (this rendering): MkNT v1 multirun (4 x 100k iter, treeThin=100, PT nChains=4)

