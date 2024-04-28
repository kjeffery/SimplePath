version: 1

scene_parameters {
    output_file_name: "image.pfm"
    width: 1350
    height: 2000
}

perspective_camera {
    origin: -1.79536 -0.0338669 130.0
    look_at: -1.79536, -0.0338669, 13.8378
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
    file: "stl_files/elf/nude-body.stl"
    #rotate: 1.0 0.0 0.0 -90.0
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
    translate: 0.0 -42.7188 0.0
}

environment_light {
    rotate: 0.0 1.0 0.0 45.0
    radiance: 0.75 0.75 0.75
}
