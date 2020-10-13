# pg-pathcomp (Pan-Genome Path Comparison scripts)

## Overview

Retrieve basic statistics from a pan-genome graph (in GFA format) by comparing paths.  Coded as part of the [SV Biohackathon](https://github.com/collaborativebioinformatics/swagg)

## Build

```
git clone https://github.com/glennhickey/pg-pathcomp.git --recursive
make
```

## Use

### Create a heatmap that shows (reference-free, SV aware) similarity between all samples in the graph
(from the `pg-pathcomp/` directory)
```
. venv/bin/activate
bin/pg-pathcomp graph.gfa > path-comp.tsv
bin/plotHeatmap.py path-comp.tsv path-comp.pdf
```
