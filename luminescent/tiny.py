import itertools
import os
import luminescent as lumi
import gdsfactory as gf

# dir = "runs"
dir = os.path.join("build", "precompile_execution")
# lumi.solve(os.path.join(dir, "tiny_2_float32_None"))

c = gf.Component()
wg = gf.components.straight(length=2,)
ext = gf.components.straight(length=.5)
_ext = c.add_ref(ext)
_wg = c.add_ref(wg)
_ext.connect("o2", _wg.ports["o1"])
c.add_port("o1", port=_wg.ports["o1"])
c.add_port("o2", port=_wg.ports["o2"])

c << gf.components.bbox(component=c, layer=lumi.BBOX, top=1, bottom=1)


c.show()

wavelengths = [1, 1.5]
margins = [[0, 0], [1, 1]]

for N, dtype, gpu in itertools.product([2, 3], ["float32"], [None, "CUDA"]):
    # for N, dtype, gpu in itertools.product([2, 3], ["float32"], [None]):
    path = os.path.join(dir, f"a_{N}_{dtype}_{gpu}")
    approx_2D_mode = "TE" if N == 2 else None
    lumi.make_prob(path, c, wavelengths=wavelengths, keys=[
        "2,1"], nres=4, approx_2D_mode=approx_2D_mode, gpu=gpu, dtype=dtype)


# raise NotImplementedError("This is a stub")
c = gf.Component()
d = lumi.mimo(west=1, east=1, l=.5, w=.5,  wwg=.5)
ext = gf.components.straight(length=.5)
_ext = c.add_ref(ext)
_d = c.add_ref(d)
_ext.connect("o2", _d.ports["o1"])
c.add_port("o1", port=_d.ports["o1"])
c.add_port("o2", port=_d.ports["o2"])
c << gf.components.bbox(component=c, layer=lumi.BBOX, top=1, bottom=1)
c.show()

targets = {"tparams": {
    1.5: {
        "2,1": 1.0
    }}}
for dtype in ['float16']:
    path = os.path.join(dir, f"b_{dtype}")
    lumi.make_design_prob(
        path,  c, targets,
        fill_material="Si", void_material="SiO2",
        lvoid=0.2, iters=2, nres=4,
        approx_2D_mode="TE", dtype=dtype)


# lumi.finetune(os.path.join("runs", "back"), iters=2)
0
