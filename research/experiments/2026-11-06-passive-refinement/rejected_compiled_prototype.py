"""Optional strict, serial compilation of the passive diagnostic pair calculation."""
import numpy as np
from numba import njit


@njit(fastmath=False,cache=False)
def ordered_pairs(positions,radii,adhesion,radius,repulsion_strength,damping,dt):
    count=len(positions)
    forces=np.zeros((count,3),dtype=np.float64)
    for i in range(count):
        for j in range(i+1,count):
            delta=positions[j]-positions[i]
            # Match numpy.linalg.norm's real 1D dot/sqrt path.
            distance=np.sqrt(np.dot(delta,delta))
            if distance<1e-8:
                continue
            direction=delta/distance
            contact=radii[i]+radii[j]
            if distance>radius:
                continue
            overlap=contact-distance
            if overlap>0.:
                repulsion=repulsion_strength*overlap
                forces[i]-=direction*repulsion
                forces[j]+=direction*repulsion
            attraction=adhesion[i]*adhesion[j]*min(1.,radius/max(distance,1e-8))
            forces[i]+=direction*attraction
            forces[j]-=direction*attraction
    velocities=np.zeros((count,3),dtype=np.float64)
    pressure=np.zeros(count,dtype=np.float64)
    for i in range(count):
        pressure[i]=np.sqrt(np.dot(forces[i],forces[i]))
        velocities[i]=forces[i]/max(damping,1e-8)
        positions[i]+=velocities[i]*dt
    return positions,velocities,pressure


def compiled_pair_step(cells,mechanics,dt):
    if len(cells)<2:
        return
    positions=np.array([c.position for c in cells],dtype=np.float64)
    radii=np.array([c.radius for c in cells],dtype=np.float64)
    adhesion=np.array([c.adhesion_strength for c in cells],dtype=np.float64)
    positions,velocity,pressure=ordered_pairs(positions,radii,adhesion,
        mechanics.interaction_radius,mechanics.repulsion_strength,mechanics.damping,dt)
    for i,cell in enumerate(cells):
        cell.position[:]=positions[i]
        cell.velocity=velocity[i].copy()
        cell.pressure=float(pressure[i])
