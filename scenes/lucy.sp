version: 1

scene_parameters {
    output_file_name: "image.pfm"
    width: 1350
    height: 2000
}

perspective_camera {
    origin: 690.756 500.0 -2000.0
    look_at: 690.756 200.0 192.627
    fov: 45
}

material_glossy {
    name: "material_glossy_base"
    diffuse: 0.7 0.7 0.7
    ior: 1.3
    roughness: 0.75
}

material_glossy {
    name: "material_glossy_plane"
    diffuse: 0.4 0.1 0.1
    ior: 1.8
    roughness: 0.01
}

material_clearcoat {
    name: "material_glossy_clearcoat"
    base: "material_glossy_base"
    ior: 1.5
    color: 1.0 1.0 1.0
}

mesh {
    file: "ply_files/lucy.ply"
    rotate: 1.0 0.0 0.0 -90.0
    ##material: "material_glossy"
    material: "material_glossy_clearcoat"
}

#sphere {
    #material: "material_glossy_base"
    #translate: 100.756 -607.893 800.531
    #scale: 100 100 100
#}

plane {
    material: "material_glossy_plane"
    translate: 0.0 -605.893 0.0
}

environment_light {
    rotate: 0.0 1.0 0.0 45.0
    radiance: 1.0 1.0 1.3
}
