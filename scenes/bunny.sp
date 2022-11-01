version: 1

scene_parameters {
    output_file_name: "image.pfm"
    width: 800
    height: 600
}

perspective_camera {
    origin: 0.0 2.0 5.0
    look_at: 0.0 1.0 0.0
    fov: 45
}

material_lambertian {
    name: "material_lambertian"
    diffuse: 0.1 0.2 0.8
}

material_glossy {
    name: "material_glossy"
    diffuse: 0.8 0.2 0.8
    ior: 1.8
    roughness: 0.01
}

material_clearcoat {
    name: "material1"
    base: "material_lambertian"
    ior: 1.5
    color: 1.0 0.8 0.8
}

material_clearcoat {
    name: "material2"
    base: "material_glossy"
    ior: 1.3
    color: 1.0 1.0 1.0
}

mesh {
    file: "ply_files/bunny/reconstruction/bun_zipper.ply"
    translate: 1.0 0.0 0.0
    scale: 10.0 10.0 10.0
    material: "material1"
}

mesh {
    file: "ply_files/bunny/reconstruction/bun_zipper.ply"
    translate: -1.0 0.0 0.0
    scale: 10.0 10.0 10.0
    material: "material2"
}

plane {
    material: "material_glossy"
    translate: 0.0 0.329874 0.0
}

environment_light {
    rotate: 0.0 1.0 0.0 45.0
    radiance: 1.0 1.0 1.3
}
