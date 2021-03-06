import * as THREE from 'three';

function Sphere(center, radius) {
    this.center = center || new THREE.Vector3();
    this.radius = radius || 1.0;
}

Sphere.prototype.constructor = Sphere;

Sphere.prototype.setCenter = function (center) {
    this.center.copy(center);
};

Sphere.prototype.setRadius = function (radius) {
    this.radius = radius;
};

var vector = new THREE.Vector3();

//
Sphere.prototype.intersectWithRayNoMiss = function (ray) {
    let pc = ray.closestPointToPoint(this.center);
    let a = pc.length(),
        d,
        b;

    // TODO: recompute mirror ray
    // If the ray miss sphere, we recompute the new ray with point symetric to tangent sphere
    if (a > this.radius) {
        // mirror point is symetric of pc
        // The mirror ray must pass through the point mirrorPoint
        const mirrorPoint = pc.clone().setLength(this.radius * 2 - a);

        // Compute the new direction
        d = ray.direction.subVectors(mirrorPoint, ray.origin).normalize();

        // Classic intersection with the new ray
        pc = ray.closestPointToPoint(this.center);
        a = pc.length();

        b = Math.sqrt(this.radius * this.radius - a * a);
        d.setLength(b);

        return vector.addVectors(pc, d);
    }

    // TODO: check all intersections : if (ray.origin.length() > this.radius)
    d = ray.direction.clone();
    b = Math.sqrt(this.radius * this.radius - a * a);
    d.setLength(b);

    return vector.subVectors(pc, d);
};

export default Sphere;
