import matplotlib.pyplot as plt
import ateam.data.make_upright as mu 
import allensdk.core.swc as swc
from os.path import splitext

def plot_cell_lims(specimen_id, scale_factor=None, scalebar=True):
    """Plot morphology from LIMS by specimen_id
    
    Parameters
    ----------
    specimen_id : int
    scale_factor : int, optional
        microns/inch, or autoscales by default
    scalebar : bool, optional
        add 100 um scale bar, by default True
    """
    nrn = mu.make_upright_morphology(specimen_id)
    plot_morph(nrn, scale_factor=scale_factor, scalebar=scalebar)

def swc_to_svg(swc_path, out_path=None, scale_factor=None, scalebar=True, transparent=True):
    out_path = out_path or splitext(swc_path)[0] + ".svg"
    nrn = swc.read_swc(swc_path)
    plot_morph(nrn, scale_factor=scale_factor, scalebar=scalebar)
    plt.savefig(out_path, transparent=transparent)

def plot_morph(nrn, scale_factor=None, scalebar=True):
    """Plot morphology from AllenSDK SWC object
    
    Parameters
    ----------
    nrn : AllenSDK SWC instance
    scale_factor : int, optional
        microns/inch, or autoscales by default
    scalebar : bool, optional
        add 100 um scale bar, by default True
    """
    fig, ax = plt.subplots()
    MORPH_COLORS = {3: "firebrick", 4: "salmon", 2: "steelblue"}
    for compartment, color in MORPH_COLORS.items():
        lines_x = []
        lines_y = []
        for c in nrn.compartment_list_by_type(compartment):
            if c["parent"] == -1:
                continue
            p = nrn.compartment_index[c["parent"]]
            lines_x += [p["x"], c["x"], None]
            lines_y += [p["y"], c["y"], None]
        plt.plot(lines_x, lines_y, c=color, linewidth=1)
    ax.set_aspect("equal")
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    if scalebar:
        bar_length = 100
        plt.plot([x1-bar_length, x1], [y0, y0], 'k', linewidth=4)
        plt.text(x1-bar_length, y0+5, '100 $\mu$m')
    if scale_factor:
        w = x1-x0
        h = y1-y0
        set_size(w/scale_factor, h/scale_factor, ax=ax)
    ax.axis('off')

def set_size(w, h, ax=None):
    """ w, h: width, height in inches """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)