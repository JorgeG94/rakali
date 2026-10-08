"""Python interface to the Rakali shallow water solver.

Usage (from namelist):
    from rakali import Rakali

    sim = Rakali("river_plume.nml")
    sim.advance(3600.0)       # run for 1 hour
    h = sim.h_phys            # numpy view of water depth (physical domain)
    print(f"Max depth: {h.max():.2f} m")

Usage (step-by-step):
    sim = Rakali("dam_break.nml")
    for _ in range(100):
        sim.step()
        print(f"t = {sim.time:.3f}, dt = {sim.dt:.4f}")
"""
from rakali._core import Rakali
from rakali.unstr import UnstructuredMesh, UnstructuredSimulation

__version__ = "0.1.0"
__all__ = [
    "Rakali",
    "UnstructuredMesh",
    "UnstructuredSimulation",
    "__version__",
]
