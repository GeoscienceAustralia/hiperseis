#!/usr/bin/env bash

# gather
mpirun --mca mpi_warn_on_fork 0 cluster gather tests/mocks/events/ -w "P S"

# sort and filter
cluster sort outfile_P.csv 5. -s sorted_P.csv
cluster sort outfile_S.csv 10. -s sorted_S.csv

# match `P` with `S` counterparts
cluster match sorted_P.csv sorted_S.csv -p matched_P.csv -s matched_S.csv

# zones
cluster zone -z '0 -50.0 100 190' matched_P.csv -r region_P.csv -g global_P.csv
cluster zone -z '0 -50.0 100 190' matched_S.csv -r region_S.csv -g global_S.csv
