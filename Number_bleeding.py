"""
Number_bleeding.py
==================
Bleeding COUNT metric.

Where `Volume_bleeding.py` counts how many vessel *voxels* a probe destroys,
this module counts how many *distinct vessel objects* it touches. The binary
volume is first split into 3-D connected components ("labelled"), and an
insertion is then scored by the number of different labels that fall inside
the probe body.

    raw binary tif
        -> vessel_seg()                  label components, drop tiny ones
        -> process_cone_positions_num()  count distinct labels per insertion

Geometry, axis order and units come from `Volume_bleeding.py`; see that
module's docstring for the (depth, height, width) = (z, x, y) convention and
the 1 um isotropic sampling. The geometry helpers are imported from there
rather than redefined, so there is exactly one implementation of each.

Importing this module does not pull in PyVista. The only function that needs
it is the debug helper `visualize_labeled_3d`, which imports it on call.
"""

import numpy as np
from scipy.ndimage import label, sum as ndi_sum, binary_erosion, generate_binary_structure
from scipy.spatial import cKDTree
from collections import defaultdict

# Single source of truth for the probe geometry.
from Volume_bleeding import create_cone_mask, calculate_cone_radius


# ========================================
# Label Vessels & Filter
# ========================================
def label_vessel(img, min_size, connectivity):
    """Split the binary volume into 3-D connected components and drop small ones.

    Parameters
    ----------
    img : ndarray, shape (depth, height, width)
        Binary vasculature volume (non-zero = vessel).
    min_size : int
        Components with fewer than this many voxels are removed. The pipeline
        uses 2, i.e. only isolated single voxels are discarded (about 2,600
        per volume, roughly 0.001 % of the foreground).
    connectivity : int
        2 selects a 3x3x3 structuring element (26-connectivity, diagonals
        count as connected); anything else falls back to the SciPy default of
        6-connectivity (faces only).

    Returns
    -------
    filtered : ndarray of int32, same shape as `img`
        Label image: 0 = background, n = component n. Labels are not
        renumbered after filtering, so the surviving ids are sparse.
    n_kept : int
        Number of components that passed the size filter.
    """
    print("    Labelling the vessels...")
    structure = np.ones((3, 3, 3)) if connectivity == 2 else None
    labeled, num_features = label(img, structure=structure)
    # For a binary image, summing the image over each label equals its voxel count.
    sizes = ndi_sum(img, labeled, index=np.arange(1, num_features + 1))
    keep = np.where(sizes >= min_size)[0] + 1          # +1: label ids are 1-based
    mask = np.isin(labeled, keep)
    filtered = np.where(mask, labeled, 0)

    print(f"    Total number of vessel: {len(keep)}")
    return filtered, len(keep)


def vessel_seg(img_data, min_size=100, connectivity=2, distance=2):
    """Turn a binary volume into a label image ready for the count metric.

    Thin wrapper around `label_vessel`. The `merge_labels` step below is left
    disabled, so `distance` is currently ignored.

    Parameters
    ----------
    img_data : ndarray
        Binary vasculature volume.
    min_size : int
        Minimum component size in voxels; the notebooks pass 2.
    connectivity : int
        2 for 26-connectivity.
    distance : float
        Would-be merge distance; unused.

    Returns
    -------
    ndarray of int32
        Label image.
    """
    filtered, count = label_vessel(img_data, min_size, connectivity)
    # merged = merge_labels(filtered, distance)   # see merge_labels: does not scale
    return filtered


# ========================================
# Label merging - written but not used
# ========================================
def get_all_surface_voxels_vectorized(labeled_data):
    """Collect the surface voxels of every label.

    A voxel is on the surface if it is foreground but disappears after a
    single binary erosion. Restricting the later distance test to surface
    voxels keeps the KD-tree small.

    Returns
    -------
    dict[int, ndarray]
        Label id -> (n, 3) array of surface voxel coordinates.
    """
    structure = generate_binary_structure(3, 1)
    mask = labeled_data > 0
    eroded = binary_erosion(mask, structure=structure, border_value=0)
    surface_mask = mask & (~eroded)
    surface_coords = np.argwhere(surface_mask)
    surface_labels = labeled_data[tuple(surface_coords.T)]

    label_surface_pts = defaultdict(list)
    for lbl in np.unique(surface_labels):
        label_surface_pts[lbl] = surface_coords[surface_labels == lbl]
    return label_surface_pts


def merge_labels(labeled_data, distance=2):
    """Merge labels whose surfaces come within `distance` voxels of each other.

    Intended to repair vessels broken into pieces by imperfect segmentation.
    Bounding boxes are compared first as a cheap rejection test, then a
    KD-tree checks the actual surface voxels, and a union-find structure
    records which labels should collapse together.

    NOT USED by the pipeline. The pairwise loop is O(n^2) over labels, and
    these volumes contain roughly 41,000 of them, so it does not finish in
    reasonable time.

    Parameters
    ----------
    labeled_data : ndarray
        Label image.
    distance : float
        Maximum surface-to-surface gap, in voxels, that still merges.

    Returns
    -------
    ndarray
        Label image with merged components sharing one id.
    """
    labels = np.unique(labeled_data)
    labels = labels[labels != 0]

    # Bounding box per label - used to skip pairs that cannot possibly be close.
    bounding_boxes = {}
    for lbl in labels:
        coords = np.argwhere(labeled_data == lbl)
        bounding_boxes[lbl] = (coords.min(axis=0), coords.max(axis=0))

    label_surface_pts = get_all_surface_voxels_vectorized(labeled_data)

    # Union-find: parent[x] points at a representative of x's group.
    parent = {lbl: lbl for lbl in labels}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]      # path halving
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[ry] = rx

    label_list = list(labels)
    for i in range(len(label_list)):
        for j in range(i + 1, len(label_list)):
            lbl_i, lbl_j = label_list[i], label_list[j]
            min_i, max_i = bounding_boxes[lbl_i]
            min_j, max_j = bounding_boxes[lbl_j]
            # Cheap test: do the bounding boxes overlap once padded by `distance`?
            if np.all(max_i + distance >= min_j) and np.all(max_j + distance >= min_i):
                pts_i = label_surface_pts[lbl_i]
                pts_j = label_surface_pts[lbl_j]
                if len(pts_i) > 0 and len(pts_j) > 0:
                    if any(cKDTree(pts_i).query_ball_point(pts_j, r=distance)):
                        union(lbl_i, lbl_j)

    label_map = {lbl: find(lbl) for lbl in labels}
    out = labeled_data.copy()
    for lbl in labels:
        out[labeled_data == lbl] = label_map[lbl]

    print(f"    Vessels after merging: {len(np.unique(out)) - 1}")
    return out


# ========================================
# Debug helpers
# ========================================
def get_top_n_labels(labeled_data, n=30):
    """Return the ids of the `n` largest components, for quick inspection."""
    labels, counts = np.unique(labeled_data[labeled_data > 0], return_counts=True)
    return labels[np.argsort(counts)[-n:]]


def visualize_labeled_3d(labeled_data):
    """Render every label as its own randomly coloured mesh (debug helper).

    Opens a blocking PyVista window, one mesh per label, so pass a volume that
    has already been reduced with `get_top_n_labels`. PyVista is imported here
    rather than at module level so the metrics do not depend on it.
    """
    import pyvista as pv

    pl = pv.Plotter()
    for lbl in np.unique(labeled_data):
        if lbl == 0:
            continue
        mesh = pv.wrap((labeled_data == lbl).astype(np.uint8)).threshold(0.5)
        pl.add_mesh(mesh, color=np.random.rand(3), opacity=1.0, show_scalar_bar=False)
    pl.show()


# ========================================
# Simulation
# ========================================
def find_first_black_pixel_slice(img_data, x_center, y_center):
    """Cortical surface depth under (x_center, y_center); 0 if never found.

    Same idea as `Volume_bleeding.find_first_black_pixel_slice`, but returns 0
    rather than None for a column that is foreground all the way down, so the
    caller always gets a usable index. That is the only difference between the
    two, and it is why this one is not simply imported.
    """
    for z in range(img_data.shape[0]):
        if img_data[z, x_center, y_center] == 0:
            return z
    return 0


def simulate_cone_insertion_num(img_data, x_center, y_center, shank_length, shank_base_diameter,
                                shank_top_diameter, tip_length, tip_base_diameter, tip_top_diameter,
                                start_slice, depth_limit):
    """Number of distinct labels touched by one insertion.

    Same geometry walk as `Volume_bleeding.simulate_cone_insertion`, but
    instead of adding up voxels it collects label ids into a set, so a vessel
    crossed at many depths still counts once.

    Parameters
    ----------
    img_data : ndarray, shape (depth, height, width)
        LABEL image from `vessel_seg` - not the raw binary volume.
    x_center, y_center : int
        Insertion site, in voxels.
    shank_length ... tip_top_diameter : float
        Probe geometry, in voxels.
    start_slice : int
        Cortical surface depth for this column.
    depth_limit : int
        Maximum insertion depth, in voxels.

    Returns
    -------
    int
        Number of distinct non-zero labels intersected.
    """
    depth, height, width = img_data.shape
    contacted_labels = set()

    for z in range(start_slice, start_slice + depth_limit):
        if z >= depth:
            break                                   # bounds guard

        relative_depth = z - start_slice
        if relative_depth < shank_length:
            radius = calculate_cone_radius(relative_depth, shank_base_diameter, shank_top_diameter, shank_length)
        elif relative_depth < shank_length + tip_length:
            radius = calculate_cone_radius(relative_depth - shank_length, tip_base_diameter, tip_top_diameter, tip_length)
        else:
            break

        slice_masked = img_data[z][create_cone_mask(height, width, x_center, y_center, radius)]

        # Set semantics: re-touching a vessel at another depth adds nothing.
        unique_labels = np.unique(slice_masked)
        contacted_labels.update(unique_labels[unique_labels != 0])

    return len(contacted_labels)


def process_cone_positions_num(img_data, y_center, shank_length, shank_base_diameter,
                               shank_top_diameter, tip_length, tip_base_diameter, tip_top_diameter,
                               depth_limit, pos_num):
    """Count metric repeated at `pos_num` evenly spaced sites along x.

    Mirror of `Volume_bleeding.process_cone_positions`; the only differences
    are that `img_data` is a label image and the per-site value is a count of
    distinct labels rather than a voxel sum.

    Returns
    -------
    ndarray of int, shape (pos_num,)
    """
    print("    Counting the number of intersected vessels...")
    height = img_data.shape[1]
    space = int(height / (pos_num + 1))
    counts = np.zeros(pos_num, dtype=int)

    for pos in range(pos_num):
        x_center = space * (pos + 1)
        start_slice = find_first_black_pixel_slice(img_data, x_center, y_center)
        counts[pos] = simulate_cone_insertion_num(
            img_data, x_center, y_center, shank_length, shank_base_diameter,
            shank_top_diameter, tip_length, tip_base_diameter, tip_top_diameter,
            start_slice, depth_limit
        )

    return counts


# ========================================
# Example
# ========================================
if __name__ == "__main__":
    # Standalone smoke test: one animal, Blackrock UEA geometry, 100 sites.
    # Reads only the first 2000 depth slices to keep it quick.
    import os
    import tifffile as tiff

    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mouse_data', 'Reslice of 9.tif')
    img = tiff.imread(path).astype(np.uint16)
    img = np.transpose(img, axes=(0, 2, 1))            # (z, y, x) -> (z, x, y)
    seg = vessel_seg(img[0:2000], min_size=2, connectivity=2)

    counts = process_cone_positions_num(
        seg, y_center=50,
        shank_length=950, shank_base_diameter=90, shank_top_diameter=28,
        tip_length=50, tip_base_diameter=28, tip_top_diameter=3,
        depth_limit=1000, pos_num=100)

    print(counts)
