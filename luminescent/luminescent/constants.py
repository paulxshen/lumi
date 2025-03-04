from copy import deepcopy
from gdsfactory.generic_tech import get_generic_pdk
from gdsfactory.technology import (
    LayerLevel,
    LayerStack,
    LayerView,
    LayerViews,
    LayerMap,
)
CANVAS_LAYER = (1000, 1)

pdk = get_generic_pdk()
pdk.activate()
LAYER_VIEWS = pdk.layer_views
LAYER_VIEWS.layer_views["WGCLAD"].visible = True
LAYER_VIEWS.layer_views["canvas"] = deepcopy(LAYER_VIEWS.layer_views["WGCLAD"])
LAYER_VIEWS.layer_views["canvas"].layer = CANVAS_LAYER
# LayerView(layer=CANVAS_LAYER, visible=True)

eps0 = 8.854187817e-12
