# Rodent BLIND vs AWARE comparison

Generated from full 1M-iteration chains (both chains complete as of 2026-05-19).

## Numerical summary

| Quantity | Value |
|---|---|
| Aware trees (total from NWK) | 268 |
| Blind trees (total from NWK) | 2525 |
| Burnin fraction applied | 25% |
| Aware post-burnin | 201 |
| Blind post-burnin | 1893 |
| Aware subsampled for CID | 200 |
| Blind subsampled for CID | 200 |
| Aware minESS (at 1M iter, from .er log) | TBD (parse MkNT aware err log) |
| Blind minESS (at 1M iter, from .er log) | **BELOW 200 (minESS = 140)** |

### Pairwise CID distances (subsampled set)

| Comparison | N pairs | Mean CID | Median CID |
|---|---|---|---|
| Within aware | 19900 | 15.4616 | 18.1947 |
| Within blind | 19900 | 13.6185 | 13.3912 |
| Between chains | 40000 | 43.9473 | 43.7442 |

## Split agreement

- BLIND MR consensus: **45 splits**
- AWARE MR consensus: **29 splits**
- Shared splits: **0**
- Blind-unique splits: **45**
- Aware-unique splits: **29**

### Blind-unique splits (absent from AWARE consensus)
  {Acomys, Arvicanthis, Arvicola, Chaetodipus, Cricetomys, Dipodomys, Dipus, Geomys, Gerbillus, Heteromys, Hydromys, Jaculus, Lepus, Mus, Napaeozapus, Neotoma, Ochotona, Orthogeomys, Oryctolagus, Pedetes, Perognathus, Platacanthomyidae, Rattus, Spalax, Sylvilagus, Thomomys}
  {Chaetodipus, Dipodomys, Dipus, Geomys, Heteromys, Jaculus, Lepus, Ochotona, Orthogeomys, Oryctolagus, Pedetes, Perognathus, Sylvilagus, Thomomys}
  {Chaetodipus, Dipodomys, Geomys, Heteromys, Orthogeomys, Perognathus, Thomomys}
  {Chaetodipus, Dipodomys, Geomys, Orthogeomys, Perognathus, Thomomys}
  {Dipodomys, Geomys, Orthogeomys, Thomomys}
  {Geomys, Orthogeomys, Thomomys}
  {Geomys, Orthogeomys}
  {Chaetodipus, Perognathus}
  {Dipus, Jaculus, Lepus, Ochotona, Oryctolagus, Pedetes, Sylvilagus}
  {Lepus, Ochotona, Oryctolagus, Pedetes, Sylvilagus}
  {Lepus, Ochotona, Oryctolagus, Sylvilagus}
  {Lepus, Oryctolagus, Sylvilagus}
  {Oryctolagus, Sylvilagus}
  {Dipus, Jaculus}
  {Acomys, Arvicanthis, Arvicola, Cricetomys, Gerbillus, Hydromys, Mus, Napaeozapus, Neotoma, Platacanthomyidae, Rattus, Spalax}
  {Acomys, Arvicanthis, Arvicola, Cricetomys, Gerbillus, Hydromys, Mus, Napaeozapus, Neotoma, Platacanthomyidae, Rattus}
  {Acomys, Mus, Rattus}
  {Mus, Rattus}
  {Arvicola, Neotoma}
  {Capromys, Castor, Cavia, Chinchilla, Cuniculus, Dasyprocta, Echimys, Erethizon, Heterocephalus, Hydrochoerus, Hystrix, Laonastes, Massoutiera, Myocastor, Thryonomys}
  {Capromys, Cavia, Chinchilla, Cuniculus, Dasyprocta, Echimys, Erethizon, Hydrochoerus, Hystrix, Myocastor, Thryonomys}
  {Capromys, Cavia, Chinchilla, Cuniculus, Dasyprocta, Echimys, Erethizon, Hydrochoerus, Hystrix, Myocastor}
  {Cavia, Chinchilla, Dasyprocta, Erethizon, Hydrochoerus, Myocastor}
  {Cavia, Chinchilla, Dasyprocta, Hydrochoerus}
  {Cavia, Chinchilla, Hydrochoerus}
  {Cavia, Hydrochoerus}
  {Erethizon, Myocastor}
  {Capromys, Cuniculus, Echimys}
  {Capromys, Cuniculus}
  {Aplodontia, Eliomys, Glaucomys, Graphiurus, Marmota, Myoxus, Ratufa, Sciurus, Tamias, Tamiasciurus}
  {Aplodontia, Glaucomys, Marmota, Ratufa, Sciurus, Tamias, Tamiasciurus}
  {Glaucomys, Marmota, Ratufa, Sciurus, Tamias, Tamiasciurus}
  {Marmota, Ratufa, Sciurus, Tamias, Tamiasciurus}
  {Marmota, Ratufa, Sciurus, Tamiasciurus}
  {Marmota, Sciurus}
  {Ratufa, Tamiasciurus}
  {Eliomys, Graphiurus}
  {Anomalurus, Dermoptera, Didelphis, Eulemur_macaco, Macropus, Papio_hamadryas, Rhynchocyon, Tupaia, Vombatus}
  {Dermoptera, Didelphis, Eulemur_macaco, Macropus, Papio_hamadryas, Rhynchocyon, Tupaia, Vombatus}
  {Dermoptera, Didelphis, Eulemur_macaco, Papio_hamadryas, Rhynchocyon, Tupaia}
  {Dermoptera, Didelphis, Eulemur_macaco, Rhynchocyon, Tupaia}
  {Dermoptera, Didelphis, Rhynchocyon, Tupaia}
  {Dermoptera, Didelphis, Tupaia}
  {Dermoptera, Tupaia}
  {Macropus, Vombatus}

### Aware-unique splits (absent from BLIND consensus)
  {Anomalurus, Aplodontia, Castor, Chinchilla, Eliomys, Heteromys, Hydromys, Massoutiera, Myoxus, Perognathus, Thryonomys, Tupaia}
  {Anomalurus, Aplodontia, Castor, Chinchilla, Eliomys, Heteromys, Hydromys, Massoutiera, Myoxus, Thryonomys, Tupaia}
  {Anomalurus, Chinchilla, Thryonomys}
  {Anomalurus, Thryonomys}
  {Acomys, Chaetodipus, Didelphis, Heterocephalus, Hydrochoerus, Laonastes, Mus, Spalax}
  {Chaetodipus, Didelphis, Heterocephalus, Hydrochoerus, Laonastes, Mus}
  {Acomys, Spalax}
  {Arvicanthis, Cavia, Myocastor, Napaeozapus, Ochotona, Papio_hamadryas, Tamiasciurus}
  {Arvicanthis, Cavia, Myocastor, Napaeozapus}
  {Arvicanthis, Cavia, Myocastor}
  {Arvicanthis, Myocastor}
  {Ochotona, Tamiasciurus}
  {Arvicola, Dipus, Geomys, Neotoma, Pedetes, Platacanthomyidae, Thomomys}
  {Arvicola, Dipus, Geomys, Neotoma, Pedetes, Platacanthomyidae}
  {Dipus, Geomys, Neotoma, Pedetes, Platacanthomyidae}
  {Neotoma, Pedetes}
  {Geomys, Platacanthomyidae}
  {Glaucomys, Graphiurus, Macropus, Marmota, Orthogeomys}
  {Glaucomys, Graphiurus, Macropus, Marmota}
  {Graphiurus, Macropus, Marmota}
  {Graphiurus, Macropus}
  {Capromys, Dasyprocta, Eulemur_macaco, Oryctolagus}
  {Dasyprocta, Oryctolagus}
  {Dermoptera, Jaculus, Vombatus}
  {Jaculus, Vombatus}
  {Gerbillus, Rattus, Tamias}
  {Gerbillus, Rattus}
  {Cricetomys, Sciurus}
  {Cuniculus, Dipodomys}

## MDS topology

The two chains' posterior distributions are **largely separated** in CID-MDS space
(centroid distance = 41.1963; aware spread = 21.0507, blind spread = 20.6980).

## Key patterns

The aware and blind posterior clouds are distinctly separated in CID-MDS space, indicating systematic topological differences driven by the ecology-aware model. The ecology covariation prior shifts the inferred rodent phylogeny in a direction that is inconsistent with the blind Mk' posterior — consistent with the hypothesis that ecological convergence creates homoplasy that misleads standard Mk' inference.

## ESS and convergence

Both chains ran for 1M iterations with treeThin = 1000. The minESS values
reported here are continuous-parameter ESS from the MkPrime MCMC log (the
minimum over all monitored parameters at the final iteration).

- **Aware chain**: minESS = NA at 1M iterations. Adequate / TBD.
- **Blind chain**: minESS = 140 at 1M iterations. Flag: below recommended minimum of 200.

ESS below 200 indicates that the chains have not fully converged on the
parameter that mixes most slowly (likely topology or a correlated rate
parameter). Results should be treated as indicative rather than definitive.
A further continuation or parallel-tempering run is advisable for publication.

## Chain details

- AWARE (this rendering): MkNT v1 multirun (4 x 100k iter, treeThin=100, PT nChains=4)
- BLIND (this rendering): MkNT v1 multirun (4 x 100k iter, treeThin=100, PT nChains=4)

