version: 1

perspective_camera {
    origin: 0.0 10.0 -10.0
    look_at: 0.0 0.0 0.0
    fov: 45
    focal_distance: 10
}

material_transmissive_dielectric {
    name: "clearcoat"
    ior: 1.3
}

material_lambertian {
    name: "base"
    diffuse: 0.1 0.2 0.8
}

# Layered materials are specified with topmost layer first
material_layered {
    name: "material0"
    layer: "clearcoat"
    layer: "base"
}

mesh {
    name: "bunny mesh"
    file: "bunny.ply"
}

mesh {
    name: "lucy mesh"
    file: "lucy.ply"
}

sphere {
    name: "my sphere"
}

plane {
    name: "my plane"
}

sphere_light {
    translate: 10.0 15.0 0.0
    radiance: 10.0 10.0 15.0
}

environment_light {
    rotate: 0.0 1.0 0.0 45.0
    radiance: 1.0 1.0 1.3
}

primitive {
    transform4x4: 1.0 0.0 0.0 0.0 # We need 16 values
                  0.0 1.0 0.0 0.0
                  0.0 0.0 1.0 0.0
                  0.0 0.0 0.0 1.0
    geometry: "bunny mesh"
    material: "material0"
}

primitive {
    # Transformations are applied in order listed in the file
    translate: 10.0 0.5 0.0
    rotate: 0.0 1.0 0.0 45.0
    scale: 10.0 10.0 10.0
    geometry: "my sphere"
    material: "base"
}

primitive {
    # Transformations are applied in order listed in the file
    geometry: "my plane"
    material: "base"
}

instance {
    # Transformations are applied in order listed in the file
    translate: 10.0 -10.5 0.0
    rotate: 0.0 1.0 0.0 45.0
    scale: 10.0 10.0 10.0
    geometry: "lucy mesh"
    material: "base"
}

instance {
    # Transformations are applied in order listed in the file
    translate: 10.0 10.5 0.0
    rotate: 0.0 1.0 0.0 -45.0
    scale: 10.0 10.0 10.0
    geometry: "lucy mesh"
    material: "base"
}
