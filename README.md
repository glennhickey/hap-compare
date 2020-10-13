# pg-pathcomp (Pan-Genome Path Comparison scripts)

## Overview

Retrieve basic statistics from a pan-genome graph (in GFA format) by comparing paths.  Hacked as part of the [SV Biohackathon](https://github.com/collaborativebioinformatics/swagg)

## Dependencies

```
sudo apt-get install build-essential git protobuf-compiler libprotoc-dev libhts-dev libjansson-dev
```
## Build

```
git clone https://github.com/glennhickey/pg-pathcomp.git --recursive
cd pg-pathcomp
make
```

## Use

### Create a heatmap that shows (reference-free, SV aware) similarity between all samples in the graph
(from the `pg-pathcomp/` directory)
```
. venv/bin/activate
pg-pathcomp -g sars.gfa -m sars-comp.tsv
plotHeatMap.py sars-comp.tsv sars-comp.pdf --log_scale
```

### Create Wiggle SV coverage on given reference path
```
vg view -Fv sars.gfa | vg snarls - > sars.snarls
pg-pathcomp -g sars.gfa -r sars.snarls -p NC_045512 -u sars.NC_045512.wig -M 2
```
