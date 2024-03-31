version: 1

scene_parameters {
    output_file_name: "image.pfm"
    width: 1000
    height: 600
}

perspective_camera {
    origin: 0.0 2.0 5.0
    look_at: -0.25 1.0 0.0
    fov: 45
}

material_lambertian {
    name: "material_lambertian"
    diffuse: 0.1 0.8 0.8
}

material_lambertian {
    name: "material_lambertian_base"
    diffuse: 0.1 0.2 0.8
}

material_glossy {
    name: "material_glossy_base"
    diffuse: 0.8 0.2 0.8
    ior: 1.8
    roughness: 0.25
}

material_glossy {
    name: "material_glossy"
    diffuse: 0.8 0.2 0.2
    ior: 1.8
    roughness: 0.75
}

material_glossy {
    name: "material_glossy_plane"
    diffuse: 0.6 0.6 0.6
    ior: 1.8
    roughness: 0.01
}

material_clearcoat {
    name: "material_lambertian_clearcoat"
    base: "material_lambertian_base"
    ior: 1.5
    color: 1.0 0.8 0.8
}

material_clearcoat {
    name: "material_glossy_clearcoat"
    base: "material_glossy_base"
    ior: 1.3
    color: 1.0 1.0 1.0
}

mesh {
    file: "ply_files/bunny/reconstruction/bun_zipper.ply"
    translate: 2.25 0.0 0.0
    scale: 10.0 10.0 10.0
    material: "material_glossy_clearcoat"
}

mesh {
    file: "ply_files/bunny/reconstruction/bun_zipper.ply"
    translate: 0.75 0.0 0.0
    scale: 10.0 10.0 10.0
    material: "material_lambertian_clearcoat"
}

mesh {
    file: "ply_files/bunny/reconstruction/bun_zipper.ply"
    translate: -0.75 0.0 0.0
    scale: 10.0 10.0 10.0
    material: "material_lambertian"
}

mesh {
    file: "ply_files/bunny/reconstruction/bun_zipper.ply"
    translate: -2.25 0.0 0.0
    scale: 10.0 10.0 10.0
    material: "material_glossy"
}

plane {
    material: "material_glossy_plane"
    translate: 0.0 0.329874 0.0
}

#environment_light {
    #rotate: 0.0 1.0 0.0 45.0
    #radiance: 1.0 1.0 1.0
#}

sphere_light {
    translate: 0.0 3.0 0.0
    scale: 0.5 0.5 0.5
    radiance: 10.0 10.0 10.0
}
