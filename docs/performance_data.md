# Rakali — performance data

Headline solver throughput on NVIDIA V100 (OpenACC, NVHPC). Setup + I/O
excluded unless noted.

## Coastal benchmarks

| Configuration | Grid | GPUs | Throughput |
|---|---|--:|--:|
| 2D barotropic | 1652×1396 | 1 | ~707 Mcells/s |
| 3D multilayer (5 σ) | 1652×1396 | 1 | ~130 Mcells/s |
| 3D multilayer bench | 2000×2000 | 1 | ~340 Mcells/s |
