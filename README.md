# locprivacy

Programs used in "Location Privacy under Continual Observation". 

The results used in the paper can be found in `/out/summary_old.csv` and `/out/summary2_old.csv`.

Please note this is not licensed for unauthorised sharing or use due to unlicensed third-party code.

# Usage

The executables are experiments.py and experiments2.py, which can be used as follows:
```
experiments [NITERS]

arguments:
  NITERS    number of iterations for each unique combination of epsilon and length
```

Output is written to `/out`. Summary files are output here, with full records of each run of the experiment written to `/out/rawdists`.
