import luminescent as lumi
import gdsfactory as gf
import numpy as np
import os


# Before creating the component, clear the component cache
gf.clear_cache()

# function to create a straight waveguide


def straight(length: float = 20, width: float = 1, layer=(1, 0)):
    c = gf.Component()
    c.add_polygon([(0, 0), (length, 0), (length, width),
                  (0, width)], layer=layer)
    c.add_port(
        name="o1", center=(0, width / 2), width=width, orientation=180, layer=layer
    )
    c.add_port(
        name="o2", center=(length, width / 2), width=width, orientation=0, layer=layer
    )
    return c


# creates straight waveguide
c = straight(length=5, width=0.05, layer=(1, 0))

# creates ring resonator
ring = gf.components.ring(radius=0.5, width=0.05,
                          angle_resolution=2.5, layer=(1, 0))

# this adds the ring to the straight waveguide
the_ring = c.add_ref(ring)

# Now that the geometry has been added to "c", we can move everything around:
the_ring.move((5/2, 0.6))

c.show()


keys = ["2,1"]  # same as keys=["o2@0,o1@0"]
nres = 6
gpu = None
materials = lumi.MATERIALS
layer_stack = lumi.SOI
wavelengths = np.linspace(1.5, 1.6, 5)


path = os.path.join("runs", "RingR")

lumi.make_prob(
    path, c, wavelengths=wavelengths, nres=nres, keys=keys, gpu=gpu,
    materials=materials, layer_stack=layer_stack)
# lumi.solve(path)
# sol=lumi.load(path)
