"""
model_3D_visualization.py
=========================
3-D rendering of one probe insertion, using PyVista.

Produces the figures saved into `output_images/`: the vasculature around the
insertion site, with the vessels that fall inside the probe highlighted.

    gold   vessels outside the probe
    red    vessels inside the probe (the ones counted as damage)
    grey   the probe itself

Rendering the full 3000 x 5000 x 100 volume is not practical, so the caller
first cuts a slab around the insertion with `cropping_img` and passes only
that.  Axis order and units follow `Volume_bleeding.py`.
"""

import pyvista as pv
import numpy as np
import warnings
warnings.filterwarnings("ignore", message="Failed to use notebook backend")


# Function to crop______________________________________________________________________________________________________________________________________________________________
def cropping_img(img_data, x_center, y_center, crop_margin_x, crop_margin_y, start_slice):
    """Cut a slab around the insertion site so the renderer stays responsive.

    Takes `crop_margin_x` voxels either side of the insertion along x and
    `crop_margin_y` either side along y, starting at the cortical surface and
    running 1000 slices deep.

    Parameters
    ----------
    img_data : ndarray, shape (depth, height, width)
        Binary vasculature volume.
    x_center, y_center : int
        Insertion site in the full volume, in voxels.
    crop_margin_x, crop_margin_y : int
        Half-width of the crop along x and y, in voxels.
    start_slice : int
        Cortical surface depth, from `find_first_black_pixel_slice`.

    Returns
    -------
    cropped_img_data : ndarray, shape (1000, 2 * crop_margin_x, 2 * crop_margin_y)
        The slab. With the notebook defaults this is (1000, 160, 100).
    adjusted_x_center, adjusted_y_center : int
        Insertion site expressed in the cropped coordinate frame.

    The crop depth is hard-coded to 1000 slices. Precondition:
    y_center +/- crop_margin_y stays inside the volume.
    """
    x_min = max(x_center - crop_margin_x , 0)
    x_max = min(x_center + crop_margin_x , img_data.shape[1])
    y_min = y_center - crop_margin_y
    y_max = y_center + crop_margin_y 
    cropped_img_data = img_data[start_slice:start_slice+1000, x_min:x_max, y_min:y_max]  # shape: [depth, H, W]

    # Adjust the cone's center to match the cropped region
    adjusted_x_center = crop_margin_x
    adjusted_y_center = crop_margin_y

    return cropped_img_data, adjusted_x_center, adjusted_y_center


# Function to calculate radius for a cone at a given depth_________________________________________________________________________________________________________________________
def calculate_cone_radius(depth, base_diameter, top_diameter, total_length):
    """Radius of a linearly tapering segment; see Volume_bleeding for details."""
    slope = (top_diameter - base_diameter) / total_length
    radius = base_diameter / 2 + slope * depth / 2
    return radius


# Function to find the first slice where the center pixel is 0______________________________________________________________________________________________________________________
def find_first_black_pixel_slice(img_data, x_center, y_center):
    """Cortical surface depth under (x_center, y_center); see Volume_bleeding."""
    depth, height, width = img_data.shape
    for z in range(depth):
        if img_data[z, x_center, y_center] == 0:
            return z
    return None  # If no slice with a black pixel is found


# Function to visualize_______________________________________________________________________________________________________________________________________________________________
def visualize_cone_pyvista(img_data, x_center, y_center, shank_length, shank_base_diameter, 
                                          shank_top_diameter, tip_length, tip_base_diameter, tip_top_diameter, 
                                          start_slice, depth_limit, ui):
    """Render the cropped slab with the probe and the vessels it intersects.

    Builds a scalar field where 1 = vessel outside the probe and 2 = vessel
    inside it, thresholds each into its own mesh, and adds a wireframe cone
    for the probe itself.  Two key bindings are attached to the window:
    "s" saves a screenshot, "r" resets the camera.

    Parameters
    ----------
    img_data : ndarray, shape (depth, height, width)
        CROPPED binary volume from `cropping_img`, not the full volume.
    x_center, y_center : int
        Insertion site in cropped coordinates.
    shank_length ... tip_top_diameter : float
        Probe geometry, in voxels.
    start_slice : int
        Accepted for signature compatibility only; overwritten with 0 below
        because `cropping_img` has already trimmed everything above the
        surface.
    depth_limit : int
        Maximum depth to draw, in voxels.
    ui : int
        1 opens an interactive window (blocking); 0 builds the scene and
        returns it so the caller can call `.screenshot(...)` headlessly.

    Returns
    -------
    pyvista.Plotter
        The populated plotter, shown already if ui == 1.
    """
    import numpy as np
    import pyvista as pv

    depth, height, width = img_data.shape
    img_data = np.transpose(img_data, (1, 2, 0))  # Reorder axes for PyVista (Y, X, Z)

    # --- Build cone mask ---
    inside_cone_mask = np.zeros_like(img_data, dtype=bool)
    
    # The crop already starts at the cortical surface, so depth is relative.
    start_slice = 0
    for z in range(start_slice, min(start_slice + depth_limit, depth)):
        relative_z = z - start_slice
        if relative_z < shank_length:
            radius = calculate_cone_radius(relative_z, shank_base_diameter, shank_top_diameter, shank_length)
        elif relative_z < shank_length + tip_length:
            radius = calculate_cone_radius(relative_z - shank_length, tip_base_diameter, tip_top_diameter, tip_length)
        else:
            break

        z_idx = relative_z
        for x in range(max(0, x_center - int(radius)), min(width, x_center + int(radius) + 1)):
            for y in range(max(0, y_center - int(radius)), min(height, y_center + int(radius) + 1)):
                if (x - x_center)**2 + (y - y_center)**2 <= radius**2:
                    inside_cone_mask[x, y, z_idx] = True

    # --- Build scalar field ---
    # 0 = tissue, 1 = vessel outside the probe, 2 = vessel inside the probe.
    scalar_field = np.zeros_like(img_data, dtype=np.uint8)
    scalar_field[img_data == 1] = 1
    scalar_field[(img_data == 1) & (inside_cone_mask)] = 2  # Inside cone

    # --- Create PyVista volume ---
    plotter = pv.Plotter()
    grid = pv.ImageData()
    grid.dimensions = scalar_field.shape
    grid.spacing = (1, 1, 1)                 # 1 um per voxel
    grid.origin = (0, 0, 0)
    grid.point_data["values"] = scalar_field.flatten(order="F")

    # Realistic: use isosurface to show vessels
    vessel_threshold = 0.5
    vessel_surface = grid.threshold(value=1, scalars="values")  # Extract vessels
    vessel_inside_surface = grid.threshold(value=1.5, scalars="values")  # Extract inside cone

    plotter.add_mesh(vessel_surface, color="#FFD700", opacity=0.5, name="vessels_outside")  # Outside cone: gold
    plotter.add_mesh(vessel_inside_surface, color="red", opacity=0.9, name="vessels_inside")  # Inside cone: red

    # --- Cone mesh ---
    # Stack of 36-gon rings, one per depth, stitched into quad faces.
    cone_points = []
    cone_faces = []
    last_radius = None

    for z in range(start_slice, min(start_slice + depth_limit, depth)):
        relative_z = z - start_slice
        if relative_z < shank_length:
            radius = calculate_cone_radius(relative_z, shank_base_diameter, shank_top_diameter, shank_length)
        elif relative_z < tip_length + shank_length:
            radius = calculate_cone_radius(relative_z - shank_length, tip_base_diameter, tip_top_diameter, tip_length)
        else:
            break

        slice_points = []
        for angle in np.linspace(0, 2 * np.pi, num=36, endpoint=False):
            x = x_center + radius * np.cos(angle)
            y = y_center + radius * np.sin(angle)
            slice_points.append([x, y, z])

        if last_radius is not None:
            # Connect this ring to the previous one with quads.
            n = len(slice_points)
            for i in range(n):
                next_i = (i + 1) % n
                cone_faces.append([4, len(cone_points) - n + i, len(cone_points) - n + next_i,
                                   len(cone_points) + next_i, len(cone_points) + i])
        cone_points.extend(slice_points)
        last_radius = radius

    cone_faces = np.array(cone_faces).flatten()
    cone_mesh = pv.PolyData(np.array(cone_points), faces=cone_faces)
    plotter.add_mesh(cone_mesh, color="gray", opacity=0.3)

    # --- Add bounding box and scalebar ---
    outline = grid.outline()
    plotter.add_mesh(outline, color="black", line_width=0.5)

    scale_length = 100                        # 100 um scale bar
    bar_start = np.array([width + 50, height - 10, 885])
    bar_end = bar_start + np.array([0, 0, scale_length])
    scale_line = pv.Line(bar_start, bar_end)
    plotter.add_mesh(scale_line, color="black", line_width=5)

    # --- Camera and lighting ---
    plotter.camera_position = 'xz'
    plotter.camera.roll += 180
    plotter.camera.azimuth = 143 + 180
    plotter.enable_lightkit()
    plotter.add_light(pv.Light(position=(10, 10, 10), intensity=1))
    plotter.add_legend([
        ['Vessels outside electrode', '#FFD700'],
        ['Vessels inside electrode', 'red'],
        ['Electrode', 'gray']
    ])

    

    # Screenshoot
    def save_screenshot():
        """Key "s": write the current view to screenshot.png."""
        plotter.screenshot("screenshot.png")
        print("Screenshot saved as screenshot.png")

    # Reset Camera
    def reset_camera():
        """Key "r": restore the default viewing angle."""
        plotter.camera_position = 'xz'
        plotter.camera.roll += 180
        plotter.camera.azimuth = 143 + 180
        plotter.reset_camera()
        print("Camera reset.")


    plotter.add_key_event("s", save_screenshot)
    plotter.add_key_event("r", reset_camera)
    

    if ui == 1:
        # auto_close=False keeps the render window usable after show() returns.
        plotter.show(auto_close=False)
    return plotter
