"""
Number_bleeding.py
==================
Bleeding COUNT metric.

Where `Volume_bleeding.py` counts how many vessel *voxels* a probe destroys,
this module counts how many *distinct vessel objects* it touches. The binary
volume is first split into 3-D connected components ("labelled"), and an
insertion is then scored by the number of different labels that fall inside
the probe body.

    binary volume
        -> load_segmentation()           label components (cached on disk)
        -> process_cone_positions_num()  count distinct labels per insertion

Geometry, axis order and units come from `Volume_bleeding.py`.
"""

import os

import numpy as np
from scipy.ndimage import label, sum as ndi_sum

from Volume_bleeding import create_cone_mask, calculate_cone_radius, load_volume


def vessel_seg(img, min_size=2, connectivity=2):
    """Split the binary volume into 3-D connected components and drop small ones.

    Parameters
    ----------
    img : ndarray, shape (depth, height, width)
        Binary vasculature volume (non-zero = vessel).
    min_size : int
        Components with fewer than this many voxels are removed. The default
        of 2 discards only isolated single voxels (about 2,600 per volume).
    connectivity : int
        2 selects a 3x3x3 structuring element (26-connectivity, diagonals
        count as connected); anything else falls back to the SciPy default of
        6-connectivity (faces only).

    Returns
    -------
    ndarray of int32, same shape as `img`
        Label image: 0 = background, n = component n. Labels are not
        renumbered after filtering, so the surviving ids are sparse.
    """
    print("    Labelling the vessels...")
    structure = np.ones((3, 3, 3)) if connectivity == 2 else None
    labeled, num_features = label(img, structure=structure)
    # For a binary image, summing the image over each label equals its voxel count.
    sizes = ndi_sum(img, labeled, index=np.arange(1, num_features + 1))
    keep = np.where(sizes >= min_size)[0] + 1          # +1: label ids are 1-based
    filtered = np.where(np.isin(labeled, keep), labeled, 0)

    print(f"    Total number of vessel: {len(keep)}")
    return filtered


def load_segmentation(tiff_path, cache_dir, img=None):
    """Label image for one animal: read it from the cache, or build and cache it.

    Labelling a full volume is the slowest step in the pipeline and yields a
    ~6 GB array, so it is done once per animal and stored as a compressed
    .npz (~15-20 MB). Delete the cache directory to force a rebuild.

    Parameters
    ----------
    tiff_path : str
        Path to the animal's binary tif stack.
    cache_dir : str
        Directory holding `<name>_seg.npz` files. Created if missing.
    img : ndarray, optional
        The already-loaded volume, to avoid reading the tif twice on a cache
        miss. Loaded from `tiff_path` if not given.

    Returns
    -------
    ndarray of int32
        Label image, shape (z, x, y).
    """
    name = os.path.splitext(os.path.basename(tiff_path))[0]
    cache = os.path.join(cache_dir, f"{name}_seg.npz")
    if os.path.exists(cache):
        return np.load(cache)["seg"]

    print(f"    Segmenting {name} (first run only, then cached) ...")
    seg = vessel_seg(load_volume(tiff_path) if img is None else img)
    os.makedirs(cache_dir, exist_ok=True)
    np.savez_compressed(cache, seg=seg)
    return seg


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
        LABEL image from `load_segmentation` - not the raw binary volume.
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
