import vispy

def sheet_view(sheet, coords=None, interactive=True, **draw_specs_kw):
    """Uses VisPy to display an epithelium
    """
    draw_specs = sheet_spec()
    spec_updater(draw_specs, draw_specs_kw)

    if coords is None:
        coords = ["x", "y", "z"]
    canvas = scene.SceneCanvas(keys="interactive", show=True, size=(1240, 720))
    view = canvas.central_widget.add_view()
    view.camera = "turntable"
    view.camera.aspect = 1
    view.bgcolor = vp.color.Color("#222222")
    if draw_specs["edge"]["visible"]:
        wire = edge_visual(sheet, coords, **draw_specs["edge"])
        view.add(wire)
    if draw_specs["face"]["visible"]:
        mesh = face_visual(sheet, coords, **draw_specs["face"])
        view.add(mesh)

    canvas.show()
    view.camera.set_range()
    if interactive:
        app.run()
    return canvas, view 