import bpy
import os
import numpy as np
import scipy.interpolate as interpol
import matplotlib as mpl
import matplotlib.cm as cm
import sys
import argparse

# function that is called when blender advances to the next frame
def frame_change_handler(scene):
    global sphere, nstep, framecount
    # we rotate the torus around the z axis
    # note that the initial position for the rotation is always the same,
    # so the rotation angle needs to increase with frame number
    # we do not use the frame count stored in 'scene' as that is not necessarily
    # an integer
    sphere.rotation_euler = (0.0, 0.0, 2.0 * np.pi / nstep * framecount)
    # advance our own framecount
    framecount += 1


# main driver routine
if __name__ == "__main__":

    # parse command line arguments
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--file", "-f", required=True)
    argparser.add_argument("--npix", "-p", type=int, default=100)
    argparser.add_argument("--nstep", "-s", type=int, default=100)
    argparser.add_argument("--nvert", "-v", type=int, default=10)
    argparser.add_argument("--resolution", "-r", type=int, default=500)
    argparser.add_argument("--output-folder", "-o", required=True)
    argparser.add_argument("--output-prefix", "-x", default="render")
    # make sure only the arguments for the Python script are parsed by argparse
    args = argparser.parse_args(sys.argv[sys.argv.index("--") + 1 :])

    npix = args.npix
    nstep = args.nstep
    nvert = args.nvert
    resolution = args.resolution
    filename = args.file
    output_folder = args.output_folder
    output_prefix = args.output_prefix

    # initialise the frame count
    framecount = 0

    # Remove all default elements in the scene
    bpy.ops.object.select_by_layer()
    bpy.ops.object.delete(use_global=False)

    # Create target for camera and lighting
    target = bpy.data.objects.new("Target", None)
    bpy.context.scene.objects.link(target)
    target.location = (0.0, 0.0, 0.0)

    # Create camera
    camera = bpy.data.cameras.new("Camera")
    camera.lens = 35
    camera.clip_start = 0.1
    camera.clip_end = 200.0
    camera.type = "PERSP"

    # Link camera to scene by creating a camera object
    camera_obj = bpy.data.objects.new("CameraObj", camera)
    camera_obj.location = (-10.0, -10.0, 10.0)
    bpy.context.scene.objects.link(camera_obj)
    bpy.context.scene.camera = camera_obj
    # make sure the camera tracks the target
    constraint = camera_obj.constraints.new("TRACK_TO")
    constraint.target = target
    constraint.track_axis = "TRACK_NEGATIVE_Z"
    constraint.up_axis = "UP_Y"

    # Set cursor to (0, 0, 0)
    bpy.context.scene.cursor_location = (0, 0, 0)

    # Create lighting
    bpy.ops.object.add(type="LAMP", location=(-10.0, -10.0, 10.0))
    lamp_obj = bpy.context.object
    lamp_obj.data.type = "SUN"
    lamp_obj.data.energy = 1
    lamp_obj.data.color = (1, 1, 1)
    # make sure the lamp tracks the target
    constraint = lamp_obj.constraints.new("TRACK_TO")
    constraint.target = target
    constraint.track_axis = "TRACK_NEGATIVE_Z"
    constraint.up_axis = "UP_Y"

    # Create the sphere
    verts = list()
    faces = list()

    # generate the vertices and faces of the sphere from a uv parametrisation
    duv = 1.0 / nvert
    u = np.linspace(0.0, 1.0 - duv, nvert)
    v = np.linspace(0.0, 1.0 - duv, nvert)
    ug, vg = np.meshgrid(u, v)
    tau = 2.0 * np.pi
    cosv = np.cos(tau * vg)
    sinv = np.sin(tau * vg)
    cosu = np.cos(tau * ug)
    sinu = np.sin(tau * ug)
    points = (6.0 * sinv * cosu, 6.0 * sinv * sinu, 6.0 * cosv)
    # we need to loop over the vertices to generate the correct faces
    for col in range(nvert):
        for row in range(nvert):
            # Surface parameterization
            point = (
                points[0][row, col],
                points[1][row, col],
                points[2][row, col],
            )
            verts.append(point)

            # Connect first and last vertices on the u and v axis
            rowNext = (row + 1) % nvert
            colNext = (col + 1) % nvert
            # Indices for each qued
            faces.append(
                (
                    (col * nvert) + rowNext,
                    (colNext * nvert) + rowNext,
                    (colNext * nvert) + row,
                    (col * nvert) + row,
                )
            )

    # Create sphere mesh and object
    mesh = bpy.data.meshes.new("SurfaceMesh")
    sphere = bpy.data.objects.new("Surface", mesh)
    sphere.location = (0.0, 0.0, 0.0)
    # Link sphere to scene
    bpy.context.scene.objects.link(sphere)
    # Create mesh from given verts and faces
    mesh.from_pydata(verts, [], faces)
    # Update mesh with new data
    mesh.update(calc_edges=True)

    # Add subsurf modifier to increase the smoothness of the sphere
    modifier = sphere.modifiers.new("Subsurf", "SUBSURF")
    modifier.levels = 2
    modifier.render_levels = 2

    # Set the sphere rendering to smooth
    mesh = sphere.data
    for p in mesh.polygons:
        p.use_smooth = True

    # read texture data
    data = np.loadtxt(filename)

    # parse the data columns
    theta = data[:, 0]
    phi = data[:, 1]
    Zsingle = data[:, 2]
    # convert theta to sin(theta - pi/2) to account for curvature
    theta = np.sin(theta - 0.5 * np.pi)

    # map the unequally spaced Z samples to a regular pixel grid
    gx, gy = np.meshgrid(
        np.linspace(phi.min(), phi.max(), npix),
        np.linspace(theta.min(), theta.max(), npix),
    )
    grid = interpol.griddata(
        (phi, theta), np.log10(Zsingle), (gx, gy), method="cubic"
    )
    # convert the mapped Z values to RGBA colours using a matplotlib colormap
    norm = mpl.colors.Normalize(
        vmin=min(np.log10(Zsingle)), vmax=max(np.log10(Zsingle))
    )
    mapper = cm.ScalarMappable(norm=norm)
    # create a texture image in blender and overwrite its pixels with the
    # colour map we just created
    textimg = bpy.data.images.new("TextureImage", width=npix, height=npix)
    pixels = mapper.to_rgba(grid).flatten()
    textimg.pixels[:] = pixels
    textimg.update()

    # Create a texture from the image
    spheretex = bpy.data.textures.new("SphereTexture", type="IMAGE")
    spheretex.image = textimg

    # Create a material for the sphere surface
    mat = bpy.data.materials.new("BumpMapMaterial")
    mat.specular_intensity = 0.0

    # Add the texture to the material
    slot = mat.texture_slots.add()
    slot.texture = spheretex
    # make sure the texture uses the coordinates of the spherical surface, so
    # that it rotates along with the sphere
    slot.texture_coords = "ORCO"
    # make sure the texture coordinates are computed using spherical projection
    slot.mapping = "SPHERE"

    # Append the material to the surface of the sphere
    sphere.data.materials.append(mat)

    # add our custom frame change handler to update the frames during the
    # animation
    bpy.app.handlers.frame_change_pre.append(frame_change_handler)

    # set the resolution of the rendered image and make sure the background is
    # transparent
    scn = bpy.context.scene
    scn.render.resolution_x = resolution
    scn.render.resolution_y = resolution
    scn.render.resolution_percentage = 100.0
    scn.render.alpha_mode = "TRANSPARENT"
    # set the number of frames in the animation
    scn.frame_end = nstep

    # Specify folder to save rendering and check if it exists
    render_folder = os.path.join(os.getcwd(), output_folder)
    if not os.path.exists(render_folder):
        os.mkdir(render_folder)

    # Render animation
    scn.render.filepath = os.path.join(render_folder, output_prefix)
    bpy.ops.render.render(animation=True)
