import platform
import subprocess
from ..constants import *
from ..layers import *
from ..utils import *
from ..materials import *
import json
import gdsfactory as gf
from copy import deepcopy
# from time import time
import datetime
from math import cos, pi, sin
import os
import numpy as np

from sortedcontainers import SortedDict, SortedSet
from gdsfactory.generic_tech import LAYER_STACK, LAYER


def setup(path, c, study, nres, center_wavelength,
          zlims=None, core_layer=LAYER.WG,
          port_margin="auto",
          runs=[],  sources=[],
          layer_stack=SOI, materials=dict(),
          default_material='SiO2',
          exclude_layers=[], Courant=None,
          gpu=None, dtype=np.float32,
          plot=False, framerate=None,
          magic="", wd=os.path.join(os.getcwd(), "runs"), name=None,
          source_port_margin=.1,
          Ttrans=None,
          approx_2D_mode=False):
    materials = {**MATERIALS, **materials}
    prob = {
        'nres': nres,
        'center_wavelength': center_wavelength,
        'name': name,
        'path': path,
    }
    if approx_2D_mode:
        N = 2
        prob["approx_2D_mode"] = approx_2D_mode
    else:
        N = 3
        prob["approx_2D_mode"] = None
    prob["class"] = "pic"
    prob["Ttrans"] = Ttrans
    prob["dtype"] = str(dtype)
    prob["timestamp"] = datetime.datetime.now().isoformat(
        timespec="seconds").replace(":", "-")
    prob["magic"] = magic
    prob["framerate"] = framerate

    gpu_backend = gpu
    # if gpu_backend:s
    prob["gpu_backend"] = gpu_backend

    ports = {
        p.name: {
            "center": (np.array(p.center)/1e0).tolist(),
            "normal": [cos(p.orientation/180*pi), sin(p.orientation/180*pi)],
            "tangent": [-sin(p.orientation/180*pi), cos(p.orientation/180*pi)],
            'width': p.width/1e0,
        }
        for p in c.get_ports_list(prefix="o")
    }
    for (k, v) in ports.items():
        z, n = [0, 0, 1], [*v['normal'], 0]
        t = np.cross(z, n).tolist()
        v['frame'] = [t, z, n]
        v["dimensions"] = [v['width']],
    prob["ports"] = ports
    mode_solutions = []

    d = get_layers(layer_stack, core_layer)[0]
    hcore = d.thickness
    zcore = d.zmin

    zmargin = 4*hcore
    zcenter = zcore+hcore/2

    h = hcore+2*zmargin
    zmin = zcore-zmargin
    zmax = zmin+h
#
    port_width = max([p.width/1e0 for p in c.ports])
    ps = portsides(c)
    xmargin = ymargin = 2*port_width

    # source_port_margin = .9*maxmargin
    # source_port_margin = 2 * port_width if N == 2 else 6*port_width
    #
    modexmargin = .6*xmargin
    modezmargin = .6*zmargin
    wmode = port_width+2*modexmargin
    hmode = hcore+2*modezmargin
    zmode = zcore-modezmargin

    prob["hcore"] = hcore
    prob["zcenter"] = zcenter
    prob["zmin"] = zmin
    prob["zmax"] = zmax
    prob["zcore"] = zcore
    # prob["L"] = [l, w, h]

    layers = set(c.layers)-set(exclude_layers)

    TEMP = os.path.join(path, "temp")
    os.makedirs(TEMP, exist_ok=True)
    BODIES = os.path.join(path, 'bodies')
    os.makedirs(BODIES, exist_ok=True)

    layer_stack_info = material_voxelate(
        c,  zmin, zmax, layers, layer_stack, BODIES)
    dir = os.path.dirname(os.path.realpath(__file__))
    print(dir)
    fn = os.path.join(dir, "solvemodes.py")
    assert os.path.exists(fn)
    assert os.path.exists(TEMP)
    if platform.system() == "Windows":
        os.system(f"copy /Y {fn} {TEMP}")
    else:
        subprocess.run(["cp", fn, TEMP])
    prob["layer_stack"] = layer_stack_info
    prob["study"] = study
    prob["materials"] = materials

    prob["N"] = N

    prob["hmode"] = hmode
    prob["wmode"] = wmode
    prob["zmode"] = zmode
    wmode = port_width+2*modexmargin
    wavelengths = []
    # _c = add_bbox(c, layer=BBOX, nonport_margin=margin)
    for run in runs:
        for k, v in list(run["sources"].items())+list(run["monitors"].items()):
            p = ports[k]
            # v['dimensions'] = p['dimensions'][0]+[hmode]
            v['dimensions'] = [wmode, hmode]
            v['frame'] = p['frame']
            for center_wavelength in v["wavelength_mode_numbers"]:
                wavelengths.append(center_wavelength)
        for k, v in run['sources'].items():
            p = ports[k]
            ct = np.array(p['center'])
            n = np.array(p['normal'])
            v["center"] = (ct+n*source_port_margin).tolist()
            v['center'] += [zcenter]
        for k, v in run['monitors'].items():
            p = ports[k]
            v['center'] = copy.deepcopy(p['center'])
            v['center'] += [zcenter]

    wavelengths = sorted(set(wavelengths))

    prob["mode_solutions"] = mode_solutions
    prob["runs"] = runs
    bbox = c.bbox_np().tolist()
    bbox[0].append(zmin)
    bbox[1].append(zmax)

    prob['bbox'] = bbox
    # prob['epdefault'] = materials[layer_stack['default']['material']]['epsilon']
    prob['epdefault'] = materials[default_material]['epsilon']

    prob["Courant"] = Courant
    if not os.path.exists(path):
        os.makedirs(path)
    # prob["wavelengths"] = wavelengths
    c.write_gds(os.path.join(path, "component.gds"))
    prob["mode_solutions"] = mode_solutions
    # prob["path_length"] = 2*(l+r+lwg)
    return prob


def port_name(port):
    s = str(port).split("@")[0]
    if s[0] == "o":
        return s
    return f"o{s}"


def port_number(port):
    s = str(port).split("@")[0]
    if s[0] == "o":
        s = s[1:]
    return int(s)


def mode_number(port):
    l = str(port).split("@")
    return 0 if len(l) == 1 else int(l[1])


def unpack_sparam_key(k):
    o, i = k.split(",")
    po, pi = port_name(o), port_name(i)
    mo, mi = mode_number(o), mode_number(i)
    return po, mo, pi, mi


def long_sparam_key(k):
    po, mo, pi, mi = unpack_sparam_key(k)
    return f"{po}@{mo},{pi}@{mi}"
