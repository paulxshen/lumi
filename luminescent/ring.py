import luminescent as lumi
from gdsfactory.technology import LogicalLayer, LayerLevel, LayerStack
import gdsfactory as gf
import numpy as np
import os


# Before creating the component, clear the component cache
gf.clear_cache()

layers = {
    "core": LayerLevel(
        layer=LogicalLayer(layer=(1, 0)),
        thickness=0.18,
        zmin=0.0,
        material="Si",
        mesh_order=2,
    ),
}


SOI180 = LayerStack(layers=layers)
SOI180.layers['default'] = {
    'material': 'SiO2'
}

radius = 3.1
width = 0.4
gap = .1
extlength = 1

c = gf.Component()
wg = gf.components.straight(length=2*radius+width, width=width)
ring = gf.components.ring(radius=radius, width=width)
ext = gf.components.straight(length=extlength, width=width)

_wg = c.add_ref(wg)
_ext1 = c.add_ref(ext)
_ext2 = c.add_ref(ext)
_ext1.connect("o2", _wg.ports["o1"])
_ext2.connect("o1", _wg.ports["o2"])

_ring = c.add_ref(ring)
_ring.move((radius+width/2, width+gap+radius))
c.add_port("o1", port=wg.ports["o1"])
c.add_port("o2", port=wg.ports["o2"])

c.show()


keys = ["2,1"]  # same as keys=["o2@0,o1@0"]
nres = 3
gpu = None
materials = lumi.MATERIALS
layer_stack = SOI180
wavelengths = np.linspace(1.5, 1.6, 5)
margins = [[0, 0], [1, 1]]


path = os.path.join("runs", "ring")

lumi.make_prob(
    path, c, wavelengths=wavelengths, keys=keys, gpu=gpu,
    nres=nres, margins=margins,
    materials=materials, layer_stack=layer_stack)
# lumi.solve(path)
# sol=lumi.load(path)
