"""
Volume_bleeding.py
==================
Bleeding VOLUME metric.

Estimates how much vascular tissue a probe destroys during a straight,
perpendicular insertion into a binary vasculature volume, by counting the
number of vessel voxels that fall inside the probe body.

The probe is modelled as two stacked conical frusta (truncated cones):

        z = start_slice        <-- cortical surface, probe enters here
        |  shank    (length = shank_length, tapers shank_base -> shank_top)
        |  tip      (length = tip_length,   tapers tip_base  -> tip_top)
        v  z increasing = deeper into the brain

At every depth slice the cross-section of the probe is a circle; the number
of vessel voxels inside that circle is summed over all slices.

Input volume convention
-----------------------
Load volumes with `load_volume`, which returns shape (depth, height, width)
= (z, x, y):

    z (depth)  : 0 .. 2999   depth into the tissue, 1 um per slice
    x (height) : 0 .. 4999   medial-lateral axis, 1 um per voxel
    y (width)  : 0 .. 99     thickness of the slab, 1 um per voxel

Voxel values are binary: 1 = vessel (or space above the cortical surface),
0 = non-vascular tissue.  All lengths and diameters are in micrometres,
which equals voxels because the sampling is 1 um isotropic.

This module is the single source of truth for the probe geometry and for the
volume layout. `Number_bleeding.py` and `model_3D_visualization.py` import
from here rather than defining their own.
"""

import numpy as np


def load_volume(path):
    """Read one binary vasculature stack in the (z, x, y) layout used everywhere.

    The tif is stored as (z, y, x); every function in this project expects
    (z, x, y), so this is the one place the transpose happens.

    Returns
    -------
    ndarray of uint16, shape (depth, height, width) = (z, x, y)
    """
    import tifffile

    return np.transpose(tifffile.imread(path), axes=(0, 2, 1)).astype(np.uint16)


def create_cone_mask(height, width, x_center, y_center, radius):
    """Return a boolean (height, width) mask of the probe cross-section.

    The cross-section is a filled disc of `radius` voxels centred on
    (x_center, y_center).  `np.ogrid` broadcasts two 1-D ranges instead of
    building two full 2-D coordinate grids, which keeps this cheap even
    though it is called once per depth slice.

    Parameters
    ----------
    height, width : int
        Shape of one depth slice, i.e. img_data.shape[1:3] = (x, y).
    x_center, y_center : int
        Centre of the insertion, in voxels.
    radius : float
        Probe radius at this depth, in voxels (= micrometres).

    Returns
    -------
    ndarray of bool, shape (height, width)
    """
    i, j = np.ogrid[:height, :width]
    return (i - x_center) ** 2 + (j - y_center) ** 2 <= radius ** 2


def calculate_cone_radius(depth, base_diameter, top_diameter, total_length):
    """Radius of a linearly tapering segment at `depth` voxels from its base.

    The segment goes from `base_diameter` at depth 0 to `top_diameter` at
    depth `total_length`, so the returned radius runs from base_diameter / 2
    to top_diameter / 2.  A constant-diameter (cylindrical) segment is
    obtained by passing base_diameter == top_diameter.

    Parameters
    ----------
    depth : float
        Distance from the base of this segment, in voxels.
    base_diameter, top_diameter : float
        Diameters at the two ends of the segment, in voxels.
    total_length : float
        Length of the segment, in voxels.  Must be non-zero.

    Returns
    -------
    float
        Radius in voxels.
    """
    slope = (top_diameter - base_diameter) / total_length
    radius = base_diameter / 2 + slope * depth / 2

    return radius


def simulate_cone_insertion(img_data, x_center, y_center, shank_length, shank_base_diameter,
                            shank_top_diameter, tip_length, tip_base_diameter, tip_top_diameter,
                            start_slice, depth_limit):
    """Total number of vessel voxels swept by one insertion at (x_center, y_center).

    Walks down from `start_slice`, computes the probe radius at each depth,
    and adds up the vessel voxels inside the corresponding disc.  The walk
    stops as soon as the accumulated depth exceeds shank_length + tip_length,
    so `depth_limit` only acts as an upper bound.

    Parameters
    ----------
    img_data : ndarray, shape (depth, height, width)
        Binary vasculature volume (1 = vessel).
    x_center, y_center : int
        Insertion site, in voxels.
    shank_length, shank_base_diameter, shank_top_diameter : float
        Geometry of the upper (shank) segment, in voxels.
    tip_length, tip_base_diameter, tip_top_diameter : float
        Geometry of the lower (tip) segment, in voxels.
    start_slice : int
        Depth index of the cortical surface for this insertion; normally
        obtained from `find_first_black_pixel_slice`.
    depth_limit : int
        Maximum number of slices to walk.  Also sizes the internal buffer,
        so it must be >= shank_length + tip_length.

    Returns
    -------
    float
        Summed vessel-voxel count over the whole insertion path.

    Precondition: start_slice + shank_length + tip_length <= depth.
    """
    depth, height, width = img_data.shape

    # One accumulator per depth slice; summed at the end.
    overlap_areas = np.zeros(depth_limit)
    for z in range(start_slice, start_slice + depth_limit):
        # Work out which segment of the probe sits at this depth.
        if z - start_slice < shank_length:
            # Shank: the part closest to the cortical surface.
            radius = calculate_cone_radius(z - start_slice, shank_base_diameter, shank_top_diameter, shank_length)
        elif z - start_slice < tip_length + shank_length:
            # Tip: the deepest part; depth is measured from the tip base.
            radius = calculate_cone_radius(z - start_slice - shank_length, tip_base_diameter, tip_top_diameter, tip_length)

        else:
            # Past the end of the probe - nothing left to intersect.
            break

        # Disc covered by the probe at this depth.
        circle_mask = create_cone_mask(height, width, x_center, y_center, radius)

        # Volume is binary, so summing the masked voxels counts vessel voxels.
        overlap_area = np.sum(img_data[z][circle_mask])
        overlap_areas[z - start_slice] = overlap_area

    # Total vessel volume intersected by this single insertion.
    return np.sum(overlap_areas)


def find_first_black_pixel_slice(img_data, x_center, y_center):
    """Depth index of the cortical surface directly under (x_center, y_center).

    In these volumes everything above the cortical surface is filled with 1s,
    so the first voxel equal to 0 along the vertical column marks the point
    where the probe enters real tissue.  The surface is not flat: measured
    across the 100 insertion sites it varies from roughly z = 80 to z = 1700
    depending on the animal.

    Parameters
    ----------
    img_data : ndarray, shape (depth, height, width)
        Binary vasculature volume.
    x_center, y_center : int
        Column to probe, in voxels.

    Returns
    -------
    int or None
        Depth index of the surface, or None if the column never reaches 0.

    Only the central column is inspected; the surface depth is not resolved
    per column across the probe cross-section.
    """
    depth, height, width = img_data.shape
    for z in range(depth):
        if img_data[z, x_center, y_center] == 0:
            return z
    return None  # If no slice with a black pixel is found


def process_cone_positions(img_data, y_center, shank_length, shank_base_diameter,
                           shank_top_diameter, tip_length, tip_base_diameter, tip_top_diameter,
                           depth_limit, pos_num):
    """Repeat the insertion at `pos_num` evenly spaced sites along the x axis.

    Sites are placed at x = space, 2 * space, ... where space = height /
    (pos_num + 1), so no insertion lands on the edge of the volume.  Each
    site gets its own surface depth before the probe is dropped in.

    Parameters
    ----------
    img_data : ndarray, shape (depth, height, width)
        Binary vasculature volume.
    y_center : int
        Fixed position along the thin axis; 50 puts the probe in the middle
        of the 100-voxel-thick slab.
    shank_length ... tip_top_diameter : float
        Probe geometry, in voxels (see `simulate_cone_insertion`).
    depth_limit : int
        Maximum insertion depth in voxels.
    pos_num : int
        Number of insertion sites.

    Returns
    -------
    ndarray of float, shape (pos_num,)
        Intersected vessel volume for each insertion site.
    """

    print("    Counting the volume of intersected vessels...")
    height, width = img_data.shape[1:3]
    space = int(height / (pos_num + 1))

    sum_overlap_areas_by_position = np.zeros(pos_num)

    for pos in range(pos_num):
        x_center = space * (pos + 1)
        start_slice = find_first_black_pixel_slice(img_data, x_center, y_center)
        # print(f"Position {pos+1}: x_center = {x_center}, start_slice = {start_slice}")  # Debug statement

        # Simulate the cone insertion
        sum_overlap_areas_by_position[pos] = simulate_cone_insertion(
            img_data, x_center, y_center, shank_length, shank_base_diameter,
            shank_top_diameter, tip_length, tip_base_diameter, tip_top_diameter,
            start_slice, depth_limit
        )

    return sum_overlap_areas_by_position
