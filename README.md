# DM-HEOM Readme

DM-HEOM comprises a software suite to compute Open Quantum System 
Dynamics with the Hierarchical Equations of Motion approach [1].

DM-HEOM was developed by Matthias Noack based on physical models
and prototypes provided by Tobias Kramer, Yaroslav Zelinskyi, and Mirta 
Rodríguez. Implementation details are described in [2] (computational part)
and in [3] (equations used, physical reference spectra).

Please see the Wiki for the documentation and examples.

## Acknowledgements

DM-HEOM has been developed at the Zuse Institute Berlin (ZIB), Germany,
in the framework of:
- the DFG Grant "Realistische Simulationen photoaktiver Systeme auf 
  Hochleistungsrechnern mit Vielkernprozessoren" (Realistic Simulations of
  Photoactive Systems on HPC-systems with Many-Core Processors),
  (PIs: Dr. Tobias Kramer, Prof. Dr. Alexander Reinefeld), and 
- the Intel Parallel Computing Center (IPCC) "Research Center for Many-core
  High-Performance Computing" at ZIB (PI: Dr. Thomas Steinke).

We thank the North-German Supercomputing Alliance (HLRN) for providing
computational resources.

We further thank L. Deecke, J. Launer, L. Gaedke-Merzhäuser, and D. Reusche
for also contributing to DM-HEOM. A full list of contributors can be found
in the CONTRIBUTORS.md file.

## License and Citing

DM-HEOM is released under the 3-clause BSD License (see LICENSE file).

If you use DM-HEOM for your scientific research, we respectfully ask you
to cite [2,3].

## References

[1] Y. Tanimura, R. Kubo "Time evolution of a quantum system in contact
with a nearly Gaussian-Markoffian noise bath" J. Phys. Soc. Jpn. 1989,
58, 101.

[2] M. Noack, A. Reinefeld, T. Kramer and Th. Steinke "DM-HEOM: A
Portable and Scalable Solver-Framework for the Hierarchical Equations of
Motion", 2018 IEEE International Parallel and Distributed Processing
Symposium Workshops (IPDPSW), Vancouver, BC, 2018, pp. 947-956.
https://doi.org/10.1109/IPDPSW.2018.00149

[3] T. Kramer, M. Noack, A. Reinefeld, M. Rodríguez, Y. Zelinskyi,
"Efficient calculation of open quantum system dynamics and 
time-resolved spectroscopy with distributed memory HEOM (DM-HEOM)",
2018, Journal of Computational Chemistry.
https://doi.org/10.1002/jcc.25354

