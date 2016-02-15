GO-SMART-IRE
============

go-smart-ire is a simulation of irreversible electroporation, a minimally invasive cancer treatment (MICT), as part of the Go-Smart EU FP7 project. This Python software uses the FEniCS libraries for finite element solution. There is also an Elmer-based solver against which this compares favourably.

EXECUTION
---------

These instructions create necessary directories and execute the IRE solver in the same container that is used on the Go-Smart server.

```bash
sudo ./setup.sh
sudo ./run.sh
```

Output appears in the output folder, created by `setup.sh`.
