from gdsfactory.technology import LogicalLayer, LayerLevel, LayerStack
from gdsfactory.generic_tech.layer_map import LAYER
import gdsfactory as gf
import copy
from gdsfactory.generic_tech import LAYER_STACK, LAYER

MATKEYS = {
    "si": "cSi",
    "Si": "cSi",

    "sio2": "SiO2",
    "sin": "SiN",
    "ge": "Ge",

}


def matname(k):
    if k in MATKEYS:
        return MATKEYS[k]
    return k.capitalize()


MATERIALS = {
    "cSi": {
        'epsilon': 3.48**2},
    "SiO2": {'epsilon': 1.44**2},
    "SiN": {'epsilon': 2.0**2},
    "Ge": {'epsilon': 4.0**2},
    "Si": {'epsilon': 3.48**2},
    'PEC': {'epsilon': 10000},
}

ks = copy.deepcopy(list(MATERIALS.keys()))
for k in ks:
    MATERIALS[k.lower()] = MATERIALS[k]


nm = 1e-3
thickness_wg = 220 * nm
thickness_slab_deep_etch = 90 * nm
thickness_slab_shallow_etch = 150 * nm

sidewall_angle_wg = 0
layer_core = LogicalLayer(layer=LAYER.WG)
layer_shallow_etch = LogicalLayer(layer=LAYER.SHALLOW_ETCH)
layer_deep_etch = LogicalLayer(layer=LAYER.DEEP_ETCH)

layers = {
    "core": LayerLevel(
        layer=LogicalLayer(layer=LAYER.WG),
        thickness=thickness_wg,
        zmin=0.0,
        material="si",
        mesh_order=2,
        # sidewall_angle=sidewall_angle_wg,
        # width_to_z=0.5,
    ),
    # "shallow_etch": LayerLevel(
    #     layer=LogicalLayer(layer=LAYER.SHALLOW_ETCH),
    #     thickness=thickness_wg - thickness_slab_shallow_etch,
    #     zmin=0.0,
    #     material="si",
    #     mesh_order=1,
    #     derived_layer=LogicalLayer(layer=LAYER.SLAB150),
    # ),
    # "deep_etch": LayerLevel(
    #     layer=LogicalLayer(layer=LAYER.DEEP_ETCH),
    #     thickness=thickness_wg - thickness_slab_deep_etch,
    #     zmin=0.0,
    #     material="si",
    #     mesh_order=1,
    #     derived_layer=LogicalLayer(layer=LAYER.SLAB90),
    # ),
    # "slab150": LayerLevel(
    #     layer=LogicalLayer(layer=LAYER.SLAB150),
    #     thickness=150e-3,
    #     zmin=0,
    #     material="si",
    #     mesh_order=3,
    # ),
    # "slab90": LayerLevel(
    #     layer=LogicalLayer(layer=LAYER.SLAB90),
    #     thickness=thickness_slab_deep_etch,
    #     zmin=0.0,
    #     material="si",
    #     mesh_order=2,
    # ),
}


SOI = LayerStack(layers=layers)
SOI.layers['default'] = {
    'material': 'SiO2'
}
BBOX_LAYER = (8888, 8888)
