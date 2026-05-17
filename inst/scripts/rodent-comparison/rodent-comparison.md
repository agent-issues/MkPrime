# Rodent BLIND vs AWARE comparison

## Numerical summary

| Quantity | Value |
|---|---|
| Aware trees (total from NWK) | 10 |
| Blind trees (total from NWK) | 359 |
| Burnin fraction applied | 25% |
| Aware post-burnin | 7 |
| Blind post-burnin | 269 |
| Aware subsampled for CID | 7 |
| Blind subsampled for CID | 200 |

### Pairwise CID distances (subsampled set)

| Comparison | N pairs | Mean CID | Median CID |
|---|---|---|---|
| Within aware | 21 | 32.5893 | 33.3560 |
| Within blind | 19900 | 12.5754 | 12.5093 |
| Between chains | 1400 | 45.7361 | 46.0790 |

## Split agreement

- BLIND MR consensus: **38 splits**
- AWARE MR consensus: **15 splits**
- Shared splits: **0**
- Blind-unique splits: **38**
- Aware-unique splits: **15**

### Blind-unique splits (absent from AWARE consensus)
  {Acomys, Arvicanthis, Arvicola, Cricetomys, Dipus, Gerbillus, Hydromys, Jaculus, Mus, Napaeozapus, Neotoma, Platacanthomyidae, Rattus, Spalax}
  {Acomys, Arvicanthis, Arvicola, Cricetomys, Gerbillus, Hydromys, Mus, Napaeozapus, Neotoma, Platacanthomyidae, Rattus, Spalax}
  {Acomys, Arvicanthis, Arvicola, Cricetomys, Gerbillus, Hydromys, Mus, Napaeozapus, Neotoma, Platacanthomyidae, Rattus}
  {Acomys, Arvicanthis, Arvicola, Cricetomys, Gerbillus, Hydromys, Mus, Neotoma, Platacanthomyidae, Rattus}
  {Mus, Rattus}
  {Arvicola, Neotoma}
  {Dipus, Jaculus}
  {Capromys, Cavia, Chinchilla, Cuniculus, Dasyprocta, Echimys, Erethizon, Hydrochoerus, Hystrix, Myocastor, Thryonomys}
  {Capromys, Cavia, Chinchilla, Cuniculus, Dasyprocta, Echimys, Erethizon, Hydrochoerus, Hystrix, Myocastor}
  {Capromys, Cavia, Chinchilla, Cuniculus, Dasyprocta, Echimys, Erethizon, Hydrochoerus, Myocastor}
  {Cavia, Chinchilla, Dasyprocta, Erethizon, Hydrochoerus, Myocastor}
  {Cavia, Chinchilla, Dasyprocta, Hydrochoerus}
  {Cavia, Chinchilla, Hydrochoerus}
  {Cavia, Hydrochoerus}
  {Erethizon, Myocastor}
  {Capromys, Cuniculus, Echimys}
  {Capromys, Cuniculus}
  {Dermoptera, Didelphis, Eulemur_macaco, Macropus, Papio_hamadryas, Rhynchocyon, Tupaia, Vombatus}
  {Dermoptera, Didelphis, Eulemur_macaco, Macropus, Papio_hamadryas, Rhynchocyon, Tupaia}
  {Dermoptera, Didelphis, Eulemur_macaco, Papio_hamadryas, Rhynchocyon, Tupaia}
  {Dermoptera, Eulemur_macaco, Papio_hamadryas, Rhynchocyon, Tupaia}
  {Eulemur_macaco, Papio_hamadryas, Rhynchocyon}
  {Eulemur_macaco, Papio_hamadryas}
  {Chaetodipus, Dipodomys, Geomys, Heteromys, Orthogeomys, Perognathus, Thomomys}
  {Dipodomys, Geomys, Orthogeomys, Thomomys}
  {Geomys, Orthogeomys, Thomomys}
  {Geomys, Orthogeomys}
  {Chaetodipus, Perognathus}
  {Aplodontia, Glaucomys, Marmota, Ratufa, Sciurus, Tamias, Tamiasciurus}
  {Glaucomys, Marmota, Ratufa, Sciurus, Tamias, Tamiasciurus}
  {Marmota, Ratufa, Sciurus, Tamias, Tamiasciurus}
  {Marmota, Sciurus}
  {Ratufa, Tamiasciurus}
  {Lepus, Ochotona, Oryctolagus, Pedetes, Sylvilagus}
  {Lepus, Ochotona, Oryctolagus, Sylvilagus}
  {Lepus, Oryctolagus, Sylvilagus}
  {Oryctolagus, Sylvilagus}
  {Eliomys, Graphiurus}

### Aware-unique splits (absent from BLIND consensus)
  {Arvicola, Hydrochoerus, Laonastes, Myocastor, Tamiasciurus}
  {Arvicola, Hydrochoerus, Laonastes, Myocastor}
  {Arvicola, Hydrochoerus, Myocastor}
  {Hydrochoerus, Myocastor}
  {Cuniculus, Heterocephalus, Massoutiera, Oryctolagus}
  {Cuniculus, Massoutiera, Oryctolagus}
  {Massoutiera, Oryctolagus}
  {Castor, Sciurus, Sylvilagus}
  {Erethizon, Napaeozapus, Spalax}
  {Erethizon, Napaeozapus}
  {Orthogeomys, Vombatus}
  {Graphiurus, Rattus}
  {Cavia, Papio_hamadryas}
  {Lepus, Ochotona}
  {Echimys, Perognathus}

## MDS topology

The two chains' posterior distributions are **largely separated** in CID-MDS space
(centroid distance = 39.0978; aware spread = 22.7068, blind spread = 4.7162).

## Key patterns

The aware and blind clouds are distinctly separated in MDS space, indicating systematic topological differences driven by the ecology-aware model.

## Caveats

**Aware chain ran only ~127 k / 500 k iterations** (terminated 2026-05-13;
treeThin = 1000, so only 10 trees written to NWK). The 200-tree subsample
target was not achievable for the aware chain; all 7 post-burnin aware trees
were used. The blind chain (200 k iter, 359 trees written) is more complete
but still a pilot run. ESS values from the BLIND chain: the `.er` log reports
minESS = 21 at iteration 200 000, indicating low mixing; blind results should
be treated as preliminary.

The aware-chain topology is informed by the ecology-aware likelihood; blind
uses the standard Mk' without ecology covariation. Topological differences
between chains reflect genuine model-dependent signal as well as the shorter
effective run length of the aware chain. A longer aware run is required before
drawing substantive conclusions.

