// This set of controls performs orbiting, dollying (zooming), and panning. It maintains
// the "up" direction as +Y, unlike the TrackballControls. Touch on tablet and phones is
// supported.
//
//    Orbit - left mouse / touch: one finger move
//    Zoom - middle mouse, or mousewheel / touch: two finger spread or squish
//    Pan - right mouse, or arrow keys / touch: three finter swipe

import * as THREE from 'three';
import Sphere from 'Core/Math/Sphere';
import CustomEvent from 'custom-event';
import AnimationPlayer, { Animation, AnimatedExpression } from 'Scene/AnimationPlayer';

var selectClick = new CustomEvent('selectClick');

// TODO:
// Recast touch for globe
// Fix target problem with pan and panoramic (when target isn't on globe)
// Fix problem with space
// Add real collision

// FIXME:
// when move globe in damping orbit, there isn't move!!

var CONTROL_STATE = {
    NONE: -1,
    ORBIT: 0,
    DOLLY: 1,
    PAN: 2,
    TOUCH_ROTATE: 3,
    TOUCH_DOLLY: 4,
    TOUCH_PAN: 5,
    MOVE_GLOBE: 6,
    PANORAMIC: 7,
};

// The control's keys
var CONTROL_KEYS = {
    LEFT: 37,
    UP: 38,
    RIGHT: 39,
    BOTTOM: 40,
    SPACE: 32,
    SHIFT: 16,
    CTRL: 17,
    S: 83,
};

var presiceSlerp = function (qb, t) {
    if (t === 0) {
        return this;
    }

    if (t === 1) {
        return this.copy(qb);
    }

    const x = this._x;
    const y = this._y;
    const z = this._z;
    const w = this._w;

    // http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/

    var cosHalfTheta = w * qb._w + x * qb._x + y * qb._y + z * qb._z;

    if (cosHalfTheta < 0) {
        this._w = -qb._w;
        this._x = -qb._x;
        this._y = -qb._y;
        this._z = -qb._z;

        cosHalfTheta = -cosHalfTheta;
    } else {
        this.copy(qb);
    }

    if (cosHalfTheta >= 1.0) {
        this._w = w;
        this._x = x;
        this._y = y;
        this._z = z;

        return this;
    }

    const sinHalfTheta = Math.sqrt(1.0 - cosHalfTheta * cosHalfTheta);
    const halfTheta = Math.atan2(sinHalfTheta, cosHalfTheta);
    const ratioA = Math.sin((1 - t) * halfTheta) / sinHalfTheta;
    const ratioB = Math.sin(t * halfTheta) / sinHalfTheta;

    this._w = (w * ratioA + this._w * ratioB);
    this._x = (x * ratioA + this._x * ratioB);
    this._y = (y * ratioA + this._y * ratioB);
    this._z = (z * ratioA + this._z * ratioB);

    this.onChangeCallback();

    return this;
};

// private members
var space = false;
var EPS = 0.000001;

// Orbit
const rotateStart = new THREE.Vector2();
const rotateEnd = new THREE.Vector2();
const rotateDelta = new THREE.Vector2();
const spherical = new THREE.Spherical(1.0, 0.01, 0);
const sphericalDelta = new THREE.Spherical(1.0, 0, 0);

// Pan
const panStart = new THREE.Vector2();
const panEnd = new THREE.Vector2();
const panDelta = new THREE.Vector2();
const panOffset = new THREE.Vector3();

const offset = new THREE.Vector3();

// Dolly
const dollyStart = new THREE.Vector2();
const dollyEnd = new THREE.Vector2();
const dollyDelta = new THREE.Vector2();
// var scale = 1.0;
const sphericalTo = new THREE.Spherical();

// Globe move
const quatGlobe = new THREE.Quaternion();
const globeTarget = new THREE.Object3D();
const movingGlobeTarget = new THREE.Vector3();
var animatedScale = 0.0;

const ctrl = {
    progress: 0,
    quatGlobe,
    qDelta: new THREE.Quaternion(),
    dampingFactor: 0.25,
};

ctrl.qDelta.presiceSlerp = presiceSlerp;
quatGlobe.presiceSlerp = presiceSlerp;

// AnimationPlayer
const player = new AnimationPlayer();
var lastRotation = [];

var dampingMoveAnimatedExpression = (function () {
    const damp = new THREE.Quaternion(0, 0, 0, 1);
    return function (root) {
        root.qDelta.presiceSlerp(damp, root.dampingFactor * 0.2);
        root.quatGlobe.multiply(root.qDelta);
    };
}());

var zoomCenterAnimatedExpression = function (root, progress) {
    root.quatGlobe.set(0, 0, 0, 1);
    root.progress = 1 - Math.pow((1 - (Math.sin((progress - 0.5) * Math.PI) * 0.5 + 0.5)), 2);
    root.quatGlobe.presiceSlerp(root.qDelta, root.progress);
};

var animationOrbitExpression = function (root, progress) {
    root.scale = 1.0 - (1.0 - root.sphericalTo.radius / root.spherical.radius) * progress;
    root.sphericalDelta.theta = root.sphericalTo.theta;
    root.sphericalDelta.phi = root.sphericalTo.phi;
};

const orbit = {
    spherical,
    sphericalDelta,
    sphericalTo,
    scale: 1,
};

var animationDampingMove = new AnimatedExpression(120, ctrl, dampingMoveAnimatedExpression, 'Damping Move');
var animationZoomCenter = new AnimatedExpression(45, ctrl, zoomCenterAnimatedExpression, 'Zoom Center');
var animationOrbit = new AnimatedExpression(30, orbit, animationOrbitExpression, 'set Orbit');
var OrbitDamping = new Animation(60, {}, 'orbit damping');

// Replace matrix float by matrix double
globeTarget.matrixWorld.elements = new Float64Array(16);
globeTarget.matrixWorldInverse = new THREE.Matrix4();
globeTarget.matrixWorldInverse.elements = new Float64Array(16);

// Pan Move
const panVector = new THREE.Vector3();

// Save last transformation
const lastPosition = new THREE.Vector3();
const lastQuaternion = new THREE.Quaternion();

// State control
var state = CONTROL_STATE.NONE;

// Initial transformation
var initialTarget;
var initialPosition;
var initialZoom;

// picking
const ptScreenClick = new THREE.Vector2();
const sizeRendering = new THREE.Vector2();

// Tangent sphere to ellispoid
const tSphere = new Sphere();
tSphere.picking = { position: new THREE.Vector3(), normal: new THREE.Vector3() };

// Special key
var keyCtrl = false;
var keyShift = false;
var keyS = false;

// Set to true to enable target helper
var enableTargetHelper = false;

// Handle function
var _handlerMouseMove;
var _handlerMouseUp;

// Pseudo collision
const radiusCollision = 50;

// SnapCamera saves transformation's camera
// It's use to globe move
function SnapCamera(camera) {
    camera.updateMatrixWorld();

    this.matrixWorld = new THREE.Matrix4();
    this.projectionMatrix = new THREE.Matrix4();
    this.invProjectionMatrix = new THREE.Matrix4();
    this.position = new THREE.Vector3();

    this.matrixWorld.elements = new Float64Array(16);
    this.projectionMatrix.elements = new Float64Array(16);
    this.invProjectionMatrix.elements = new Float64Array(16);

    this.init = function (camera)
    {
        this.matrixWorld.elements.set(camera.matrixWorld.elements);
        this.projectionMatrix.elements.set(camera.projectionMatrix.elements);
        this.position.copy(camera.position);
        this.invProjectionMatrix.getInverse(this.projectionMatrix);
    };

    this.init(camera);

    this.shot = function (objectToSnap)
    {
        objectToSnap.updateMatrixWorld();
        this.matrixWorld.elements.set(objectToSnap.matrixWorld.elements);
        this.position.copy(objectToSnap.position);
    };

    const matrix = new THREE.Matrix4();
    matrix.elements = new Float64Array(16);

    this.updateRay = function (ray, mouse)
    {
        ray.origin.copy(this.position);
        ray.direction.set(mouse.x, mouse.y, 0.5);
        matrix.multiplyMatrices(this.matrixWorld, this.invProjectionMatrix);
        ray.direction.applyProjection(matrix);
        ray.direction.sub(ray.origin).normalize();
    };
}

var snapShotCamera;

// ///////////////////////

/* globals document,window */

function GlobeControls(camera, domElement, engine) {
    player.dom = domElement;
    var scene = engine.scene;
    this.camera = camera;
    snapShotCamera = new SnapCamera(camera);

    this.domElement = (domElement !== undefined) ? domElement : document;

    // Set to false to disable this control
    this.enabled = true;

    // This option actually enables dollying in and out; left as "zoom" for
    // backwards compatibility
    this.enableZoom = true;
    this.zoomSpeed = 1.0;

    // Limits to how far you can dolly in and out ( PerspectiveCamera only )
    this.minDistance = radiusCollision;
    this.maxDistance = Infinity;

    // Limits to how far you can zoom in and out ( OrthographicCamera only )
    this.minZoom = 0;
    this.maxZoom = Infinity;

    // Set to true to disable this control
    this.enableRotate = true;
    this.rotateSpeed = 1.0;

    // Set to true to disable this control
    this.enablePan = true;
    this.keyPanSpeed = 7.0; // pixels moved per arrow key push

    // Set to true to automatically rotate around the target
    this.autoRotate = false;
    this.autoRotateSpeed = 2.0; // 30 seconds per round when fps is 60

    // How far you can orbit vertically, upper and lower limits.
    // Range is 0 to Math.PI radians.
    // TODO Warning minPolarAngle = 0.01 -> it isn't possible to be perpendicular on Globe
    this.minPolarAngle = 0.01; // radians
    this.maxPolarAngle = Math.PI * 0.5; // radians

    // How far you can orbit horizontally, upper and lower limits.
    // If set, must be a sub-interval of the interval [ - Math.PI, Math.PI ].
    this.minAzimuthAngle = -Infinity; // radians
    this.maxAzimuthAngle = Infinity; // radians

    // Set to true to disable use of the keys
    this.enableKeys = true;

    // Enable Damping
    this.enableDamping = true;

    if (enableTargetHelper) {
        this.pickingHelper = new THREE.AxisHelper(500000);
    }

    // Mouse buttons
    this.mouseButtons = {
        PANORAMIC: THREE.MOUSE.LEFT,
        ZOOM: THREE.MOUSE.MIDDLE,
        PAN: THREE.MOUSE.RIGHT,
    };

    // Radius tangent sphere
    tSphere.setRadius(engine.size);
    spherical.radius = tSphere.radius;

    sizeRendering.set(engine.width, engine.height);

    // Note A
    // TODO: test before remove test code
    // so camera.up is the orbit axis
    // var quat = new THREE.Quaternion().setFromUnitVectors(camera.up, new THREE.Vector3(0, 1, 0));
    // var quatInverse = quat.clone().inverse();

    // events
    this.changeEvent = {
        type: 'change',
    };
    this.startEvent = {
        type: 'start',
    };
    this.endEvent = {
        type: 'end',
    };

    //
    this.updateCamera = function (camera) {
        snapShotCamera.init(camera.camera3D);
        sizeRendering.width = camera.width;
        sizeRendering.height = camera.height;
    };

    this.getAutoRotationAngle = function () {
        return 2 * Math.PI / 60 / 60 * this.autoRotateSpeed;
    };

    this.getZoomScale = function () {
        return Math.pow(0.95, this.zoomSpeed);
    };

    this.rotateLeft = function (angle) {
        if (angle === undefined) {
            angle = this.getAutoRotationAngle();
        }
        sphericalDelta.theta -= angle;
    };

    this.rotateUp = function (angle) {
        if (angle === undefined) {
            angle = this.getAutoRotationAngle();
        }

        sphericalDelta.phi -= angle;
    };

    // pass in distance in world space to move left
    this.panLeft = function (distance) {
        var te = this.camera.matrix.elements;

        // get X column of matrix
        panOffset.set(te[0], te[1], te[2]);
        panOffset.multiplyScalar(-distance);

        panVector.add(panOffset);
    };

    // pass in distance in world space to move up
    this.panUp = function (distance) {
        var te = this.camera.matrix.elements;

        // get Y column of matrix
        panOffset.set(te[4], te[5], te[6]);
        panOffset.multiplyScalar(distance);

        panVector.add(panOffset);
    };

    // pass in x,y of change desired in pixel space,
    // right and down are positive
    this.mouseToPan = function (deltaX, deltaY) {
        var element = this.domElement === document ? this.domElement.body : this.domElement;

        if (this.camera instanceof THREE.PerspectiveCamera) {
            // perspective
            var position = this.camera.position;

            // var offset = position.clone().sub(this.target);
            var offset = position.clone().sub(this.getTargetCameraPosition());

            var targetDistance = offset.length();

            // half of the fov is center to top of screen
            targetDistance *= Math.tan((this.camera.fov / 2) * Math.PI / 180.0);

            // we actually don't use screenWidth, since perspective camera is fixed to screen height
            this.panLeft(2 * deltaX * targetDistance / element.clientHeight);
            this.panUp(2 * deltaY * targetDistance / element.clientHeight);
        } else if (this.camera instanceof THREE.OrthographicCamera) {
            // orthographic
            this.panLeft(deltaX * (this.camera.right - this.camera.left) / element.clientWidth);
            this.panUp(deltaY * (this.camera.top - this.camera.bottom) / element.clientHeight);
        } else {

            // camera neither orthographic or perspective
            // console.warn('WARNING: GlobeControls.js encountered an unknown camera type - this.mouseToPan disabled.');

        }
    };

    this.dollyIn = function (dollyScale) {
        if (dollyScale === undefined) {
            dollyScale = this.getZoomScale();
        }

        if (this.camera instanceof THREE.PerspectiveCamera) {
            orbit.scale /= dollyScale;
        } else if (this.camera instanceof THREE.OrthographicCamera) {
            this.camera.zoom = Math.max(this.minZoom, Math.min(this.maxZoom, this.camera.zoom * dollyScale));
            this.camera.updateProjectionMatrix();
            this.dispatchEvent(this.changeEvent);
        } else {

            // console.warn('WARNING: GlobeControls.js encountered an unknown camera type - dolly/zoom disabled.');

        }
    };

    this.dollyOut = function (dollyScale) {
        if (dollyScale === undefined) {
            dollyScale = this.getZoomScale();
        }

        if (this.camera instanceof THREE.PerspectiveCamera) {
            orbit.scale *= dollyScale;
        } else if (this.camera instanceof THREE.OrthographicCamera) {
            this.camera.zoom = Math.max(this.minZoom, Math.min(this.maxZoom, this.camera.zoom / dollyScale));
            this.camera.updateProjectionMatrix();
            this.dispatchEvent(this.changeEvent);
        } else {

            // console.warn('WARNING: GlobeControls.js encountered an unknown camera type - dolly/zoom disabled.');

        }
    };

    var getPickingPosition = (function () {
        var engineGfx = engine;
        var position;

        return function (coords)
        {
            position = engineGfx.getPickingPositionFromDepth(coords);
            engineGfx.renderScene();

            return position;
        };
    }());

    // introduction collision
    // Not use for the moment
    // eslint-disable-next-line
    var collision = function(position)
    {
        if (scene.getMap())
        {
            var coord = scene.getMap().projection.cartesianToGeo(position);
            var bbox = scene.getMap().getTile(coord).bbox;
            var delta = coord.altitude() - (bbox.top() + radiusCollision);

            if (delta < 0) {
                position.setLength(position.length() - delta);
            }
        }

        return position;
    };

    const quaterPano = new THREE.Quaternion();
    const quaterAxis = new THREE.Quaternion();
    const axisX = new THREE.Vector3(1, 0, 0);

    var update = function () {
        // MOVE_GLOBE
        // Rotate globe with mouse
        if (state === CONTROL_STATE.MOVE_GLOBE) {
            movingGlobeTarget.copy(this.getTargetCameraPosition()).applyQuaternion(quatGlobe);
            this.camera.position.copy(snapShotCamera.position).applyQuaternion(quatGlobe);
            if (ctrl.progress > 0) {
                this.camera.position.lerp(movingGlobeTarget, ctrl.progress * animatedScale);
            }
            this.camera.up.copy(movingGlobeTarget.clone().normalize());
        // PAN
        // Move camera in projection plan
        } else if (state === CONTROL_STATE.PAN) {
            this.camera.position.add(panVector);
            movingGlobeTarget.add(panVector);
        // PANORAMIC
        // Move target camera
        } else if (state === CONTROL_STATE.PANORAMIC) {
            // TODO: this part must be reworked
            this.camera.worldToLocal(movingGlobeTarget);
            var normal = this.camera.position.clone().normalize().applyQuaternion(this.camera.quaternion.clone().inverse());
            quaterPano.setFromAxisAngle(normal, sphericalDelta.theta).multiply(quaterAxis.setFromAxisAngle(axisX, sphericalDelta.phi));
            movingGlobeTarget.applyQuaternion(quaterPano);
            this.camera.localToWorld(movingGlobeTarget);
            this.camera.up.copy(movingGlobeTarget.clone().normalize());
        } else {
            // ZOOM/ORBIT
            // Move Camera around the target camera

            // TODO: test before remove test code see (Note A)
            // offset.applyQuaternion( quat );

            // get camera position in local space of target
            offset.copy(this.camera.position).applyMatrix4(globeTarget.matrixWorldInverse);

            // angle from z-axis around y-axis
            if (sphericalDelta.theta || sphericalDelta.phi) {
                spherical.setFromVector3(offset);
            }

            if (this.autoRotate && state === CONTROL_STATE.NONE) {
                this.rotateLeft(this.getAutoRotationAngle());
            }

            spherical.theta += sphericalDelta.theta;
            spherical.phi += sphericalDelta.phi;

            // restrict spherical.theta to be between desired limits
            spherical.theta = Math.max(this.minAzimuthAngle, Math.min(this.maxAzimuthAngle, spherical.theta));

            // restrict spherical.phi to be between desired limits
            spherical.phi = Math.max(this.minPolarAngle, Math.min(this.maxPolarAngle, spherical.phi));

            spherical.radius = offset.length() * orbit.scale;

            // restrict spherical.phi to be betwee EPS and PI-EPS
            spherical.makeSafe();

            // restrict radius to be between desired limits
            spherical.radius = Math.max(this.minDistance, Math.min(this.maxDistance, spherical.radius));

            offset.setFromSpherical(spherical);

            // rotate point back to "camera-up-vector-is-up" space
            // offset.applyQuaternion( quatInverse );

            this.camera.position.copy(globeTarget.localToWorld(offset));
        }

        this.camera.lookAt(movingGlobeTarget);

        if (!this.enableDamping) {
            sphericalDelta.theta = 0;
            sphericalDelta.phi = 0;
        } else {
            sphericalDelta.theta *= (1 - ctrl.dampingFactor);
            sphericalDelta.phi *= (1 - ctrl.dampingFactor);
        }

        orbit.scale = 1;
        panVector.set(0, 0, 0);

        // update condition is:
        // min(camera displacement, camera rotation in radians)^2 > EPS
        // using small-angle approximation cos(x/2) = 1 - x^2 / 8

        if (lastPosition.distanceToSquared(this.camera.position) > EPS || 8 * (1 - lastQuaternion.dot(this.camera.quaternion)) > EPS) {
            this.dispatchEvent(this.changeEvent);

            lastPosition.copy(this.camera.position);
            lastQuaternion.copy(this.camera.quaternion);
        }
        if (this.enableDamping && state === CONTROL_STATE.ORBIT && player.isStopped() && (sphericalDelta.theta > EPS || sphericalDelta.phi > EPS)) {
            window.clearInterval(this.timer);
            this.timer = window.setTimeout(() => { player.startAnimation(OrbitDamping); }, 1000 / 30);
        }
    }.bind(this);

    this.getSpace = function () {
        return space;
    };


    this.getSphericalDelta = function () {
        return sphericalDelta;
    };

    // Position object on globe
    var positionObject = (function ()
    {
        const quaterionX = new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1, 0, 0), Math.PI / 2);
        return function (newPosition, object)
        {
            object.position.copy(newPosition);
            object.lookAt(newPosition.clone().multiplyScalar(1.1));
            object.quaternion.multiply(quaterionX);
            object.updateMatrixWorld();
        };
    }());

    // set new globe target
    var setGlobleTarget = function (newPosition) {
    // Compute the new target center position

        positionObject(newPosition, globeTarget);

        globeTarget.matrixWorldInverse.getInverse(globeTarget.matrixWorld);
    };

    const cT = new THREE.Vector3();
    // update globe target
    var updateGlobeTarget = function () {
        // Get distance camera DME
        var pickingPosition = getPickingPosition();

        if (!pickingPosition) {
            return;
        }
        var distanceTarget = Math.floor(pickingPosition.distanceTo(this.camera.position));

        // Position movingGlobeTarget on DME
        cT.subVectors(movingGlobeTarget, this.camera.position);
        cT.setLength(distanceTarget);
        movingGlobeTarget.addVectors(this.camera.position, cT);

        // set new globe target
        setGlobleTarget(movingGlobeTarget);

        // update spherical from target
        offset.copy(this.camera.position);
        offset.applyMatrix4(globeTarget.matrixWorldInverse);
        spherical.setFromVector3(offset);
        state = CONTROL_STATE.NONE;
        lastRotation = [];
    };

    // Update helper
    var updateHelper = enableTargetHelper ? function (position, helper) {
        positionObject(position, helper);
        this.dispatchEvent(this.changeEvent);
    } : function () {};

    this.getPickingPositionOnSphere = function () {
        return tSphere.picking.position;
    };

    // Update radius's sphere : the sphere must cross the point
    // Return intersection with mouse and sphere
    var updateSpherePicking = (function () {
        var mouse = new THREE.Vector2();
        var ray = new THREE.Ray();

        return function (point, screenCoord) {
            tSphere.setRadius(point.length());

            mouse.x = (screenCoord.x / sizeRendering.width) * 2 - 1;
            mouse.y = -(screenCoord.y / sizeRendering.height) * 2 + 1;

            snapShotCamera.updateRay(ray, mouse);
            // pick position on tSphere
            tSphere.picking.position.copy(tSphere.intersectWithRayNoMiss(ray));
            tSphere.picking.normal = tSphere.picking.position.clone().normalize();

            // lastRotation.splice(0);
            lastRotation.push(tSphere.picking.normal);

            updateHelper.bind(this)(tSphere.picking.position, this.pickingHelper);
        };
    }());

    var onMouseMove = (function () {
        const ray = new THREE.Ray();
        const mouse = new THREE.Vector2();

        return function (event) {
            if (player.isPlayed()) {
                window.clearInterval(this.timer);
                player.stop();
            }
            if (this.enabled === false) return;

            event.preventDefault();

            if (state === CONTROL_STATE.ORBIT || state === CONTROL_STATE.PANORAMIC) {
                if (this.enableRotate === false) return;

                rotateEnd.set(event.clientX - event.target.offsetLeft, event.clientY - event.target.offsetTop);
                rotateDelta.subVectors(rotateEnd, rotateStart);

                // rotating across whole screen goes 360 degrees around
                if (!space) {
                    this.rotateLeft(2 * Math.PI * rotateDelta.x / sizeRendering.width * this.rotateSpeed);

                    // rotating up and down along whole screen attempts to go 360, but limited to 180
                    this.rotateUp(2 * Math.PI * rotateDelta.y / sizeRendering.height * this.rotateSpeed);
                } else {
                    this.rotateLeft(rotateDelta.x);

                    // rotating up and down along whole screen attempts to go 360, but limited to 180
                    this.rotateUp(rotateDelta.y);
                }

                rotateStart.copy(rotateEnd);
            } else if (state === CONTROL_STATE.DOLLY) {
                if (this.enableZoom === false) return;

                dollyEnd.set(event.clientX - event.target.offsetLeft, event.clientY - event.target.offsetTop);
                dollyDelta.subVectors(dollyEnd, dollyStart);

                if (dollyDelta.y > 0) {
                    this.dollyIn();
                } else if (dollyDelta.y < 0) {
                    this.dollyOut();
                }

                dollyStart.copy(dollyEnd);
            } else if (state === CONTROL_STATE.PAN) {
                if (this.enablePan === false) return;

                panEnd.set(event.clientX - event.target.offsetLeft, event.clientY - event.target.offsetTop);
                panDelta.subVectors(panEnd, panStart);

                this.mouseToPan(panDelta.x, panDelta.y);

                panStart.copy(panEnd);
            } else if (state === CONTROL_STATE.MOVE_GLOBE) {
                mouse.x = ((event.clientX - event.target.offsetLeft) / sizeRendering.width) * 2 - 1;
                mouse.y = -((event.clientY - event.target.offsetTop) / sizeRendering.height) * 2 + 1;

                snapShotCamera.updateRay(ray, mouse);

                var intersection = tSphere.intersectWithRay(ray);

                if (intersection) {
                    var normalizedIntersection = intersection.normalize();
                    quatGlobe.setFromUnitVectors(normalizedIntersection, tSphere.picking.normal);
                    lastRotation.push(normalizedIntersection.clone());
                    if (lastRotation.length > 2) {
                        lastRotation.splice(0, 1);
                    }
                } else {
                    onMouseUp.bind(this)();
                }
            }

            if (state !== CONTROL_STATE.NONE) {
                update();
            }
        };
    }());

    var onMouseDown = function (event) {
        player.stop().then(() => {
            if (this.enabled === false) return;

            event.preventDefault();

            snapShotCamera.shot(this.camera);

            if (event.button === this.mouseButtons.PANORAMIC) {
                if (this.enableRotate === false) return;

                if (keyCtrl) {
                    state = CONTROL_STATE.ORBIT;
                } else if (keyShift) {
                    state = CONTROL_STATE.PANORAMIC;
                } else if (keyS) {
                    // If the key 'S' is down, the engine selects node under mouse
                    selectClick.mouse = new THREE.Vector2(event.clientX - event.target.offsetLeft, event.clientY - event.target.offsetTop);

                    domElement.dispatchEvent(selectClick);
                } else {
                    snapShotCamera.shot(this.camera);

                    ptScreenClick.x = event.clientX - event.target.offsetLeft;
                    ptScreenClick.y = event.clientY - event.target.offsetTop;

                    var point = getPickingPosition(ptScreenClick);
                    lastRotation = [];
                    // update tangent sphere which passes through the point
                    if (point) {
                        ctrl.range = this.getRange();
                        updateSpherePicking.bind(this)(point, ptScreenClick);
                        state = CONTROL_STATE.MOVE_GLOBE;
                    }
                }

                rotateStart.set(event.clientX - event.target.offsetLeft, event.clientY - event.target.offsetTop);
            } else if (event.button === this.mouseButtons.ZOOM) {
                if (this.enableZoom === false) return;

                state = CONTROL_STATE.DOLLY;

                dollyStart.set(event.clientX - event.target.offsetLeft, event.clientY - event.target.offsetTop);
            } else if (event.button === this.mouseButtons.PAN) {
                if (this.enablePan === false) return;

                state = CONTROL_STATE.PAN;

                panStart.set(event.clientX - event.target.offsetLeft, event.clientY - event.target.offsetTop);
            }

            if (state !== CONTROL_STATE.NONE) {
                this.domElement.addEventListener('mousemove', _handlerMouseMove, false);
                this.domElement.addEventListener('mouseup', _handlerMouseUp, false);
                this.dispatchEvent(this.startEvent);
            }
        });
    };

    var ondblclick = function (event) {
        if (!keyCtrl && !keyShift) {
            ptScreenClick.x = event.clientX - event.target.offsetLeft;
            ptScreenClick.y = event.clientY - event.target.offsetTop;

            const point = getPickingPosition(ptScreenClick);

            if (point) {
                animatedScale = 0.6;
                this.setCenter(point, true);
            }
        }
    };

    var onMouseUp = function (/* event */) {
        if (this.enabled === false) return;

        this.domElement.removeEventListener('mousemove', _handlerMouseMove, false);
        this.domElement.removeEventListener('mouseup', _handlerMouseUp, false);
        this.dispatchEvent(this.endEvent);

        if (this.enableDamping) {
            if (state === CONTROL_STATE.ORBIT && (sphericalDelta.theta > EPS || sphericalDelta.phi > EPS)) {
                player.startAnimation(OrbitDamping).then(() => { this.stopAnimation(); });
            } else if (state === CONTROL_STATE.MOVE_GLOBE && lastRotation.length === 2 && !(lastRotation[1].equals(lastRotation[0]))) {
                ctrl.qDelta.setFromUnitVectors(lastRotation[1], lastRotation[0]);
                player.startAnimation(animationDampingMove).then(() => { this.stopAnimation(); });
            } else {
                updateGlobeTarget.bind(this)();
            }
        } else {
            updateGlobeTarget.bind(this)();
        }
    };

    var onMouseWheel = function (event) {
        player.stop();

        if (this.enabled === false || this.enableZoom === false || state !== CONTROL_STATE.NONE) return;

        event.preventDefault();
        event.stopPropagation();

        var delta = 0;

        // WebKit / Opera / Explorer 9
        if (event.wheelDelta !== undefined) {
            delta = event.wheelDelta;
        // Firefox
        } else if (event.detail !== undefined) {
            delta = -event.detail;
        }

        if (delta > 0) {
            this.dollyOut();
        } else if (delta < 0) {
            this.dollyIn();
        }

        update();

        this.dispatchEvent(this.startEvent);
        this.dispatchEvent(this.endEvent);
    };

    var onKeyUp = function (/* event*/) {
        if (this.enabled === false || this.enableKeys === false || this.enablePan === false) return;


        if (state === CONTROL_STATE.PAN)
        {
            updateGlobeTarget.bind(this)();
        }

        keyCtrl = false;
        keyShift = false;
        keyS = false;
    };

    var onKeyDown = function (event) {
        player.stop().then(() => {
            if (this.enabled === false || this.enableKeys === false || this.enablePan === false) return;
            keyCtrl = false;
            keyShift = false;

            switch (event.keyCode) {
                case CONTROL_KEYS.UP:
                    this.mouseToPan(0, this.keyPanSpeed);
                    state = CONTROL_STATE.PAN;
                    update();
                    break;
                case CONTROL_KEYS.BOTTOM:
                    this.mouseToPan(0, -this.keyPanSpeed);
                    state = CONTROL_STATE.PAN;
                    update();
                    break;
                case CONTROL_KEYS.LEFT:
                    this.mouseToPan(this.keyPanSpeed, 0);
                    state = CONTROL_STATE.PAN;
                    update();
                    break;
                case CONTROL_KEYS.RIGHT:
                    this.mouseToPan(-this.keyPanSpeed, 0);
                    state = CONTROL_STATE.PAN;
                    update();
                    break;
                // TODO Why space key, looking for movement
                case CONTROL_KEYS.SPACE:
                    space = !space;
                    // this.updateTarget();
                    update();
                    break;
                case CONTROL_KEYS.CTRL:
                    // computeVectorUp();
                    keyCtrl = true;
                    break;
                case CONTROL_KEYS.SHIFT:
                    // computeVectorUp();
                    keyShift = true;
                    break;
                case CONTROL_KEYS.S:
                    // WARNING loop !!!
                    keyS = true;
                    break;

            }
        });
    };

    var onTouchStart = function (event) {
        if (this.enabled === false) return;

        switch (event.touches.length) {

            case 1: // one-fingered touch: rotate

                if (this.enableRotate === false) return;

                state = CONTROL_STATE.TOUCH_ROTATE;

                rotateStart.set(event.touches[0].pageX, event.touches[0].pageY);
                break;

            case 2: // two-fingered touch: dolly

                if (this.enableZoom === false) return;

                state = CONTROL_STATE.TOUCH_DOLLY;

                var dx = event.touches[0].pageX - event.touches[1].pageX;
                var dy = event.touches[0].pageY - event.touches[1].pageY;
                var distance = Math.sqrt(dx * dx + dy * dy);
                dollyStart.set(0, distance);
                break;

            case 3: // three-fingered touch: this.mouseToPan

                if (this.enablePan === false) return;

                state = CONTROL_STATE.TOUCH_PAN;

                panStart.set(event.touches[0].pageX, event.touches[0].pageY);
                break;

            default:

                state = CONTROL_STATE.NONE;

        }

        if (state !== CONTROL_STATE.NONE) this.dispatchEvent(this.startEvent);
    };

    var onTouchMove = function (event) {
        if (this.enabled === false) return;

        event.preventDefault();
        event.stopPropagation();

        var element = this.domElement === document ? this.domElement.body : this.domElement;

        switch (event.touches.length) {

            case 1: // one-fingered touch: rotate

                if (this.enableRotate === false) return;
                if (state !== CONTROL_STATE.TOUCH_ROTATE) return;

                rotateEnd.set(event.touches[0].pageX, event.touches[0].pageY);
                rotateDelta.subVectors(rotateEnd, rotateStart);

                // rotating across whole screen goes 360 degrees around
                this.rotateLeft(2 * Math.PI * rotateDelta.x / element.clientWidth * this.rotateSpeed);
                // rotating up and down along whole screen attempts to go 360, but limited to 180
                this.rotateUp(2 * Math.PI * rotateDelta.y / element.clientHeight * this.rotateSpeed);

                rotateStart.copy(rotateEnd);

                update();
                break;

            case 2: // two-fingered touch: dolly

                if (this.enableZoom === false) return;
                if (state !== CONTROL_STATE.TOUCH_DOLLY) return;

                var dx = event.touches[0].pageX - event.touches[1].pageX;
                var dy = event.touches[0].pageY - event.touches[1].pageY;
                var distance = Math.sqrt(dx * dx + dy * dy);

                dollyEnd.set(0, distance);
                dollyDelta.subVectors(dollyEnd, dollyStart);

                if (dollyDelta.y > 0) {
                    this.dollyOut();
                } else if (dollyDelta.y < 0) {
                    this.dollyIn();
                }

                dollyStart.copy(dollyEnd);

                update();
                break;

            case 3: // three-fingered touch: this.mouseToPan

                if (this.enablePan === false) return;
                if (state !== CONTROL_STATE.TOUCH_PAN) return;

                panEnd.set(event.touches[0].pageX, event.touches[0].pageY);
                panDelta.subVectors(panEnd, panStart);

                this.mouseToPan(panDelta.x, panDelta.y);

                panStart.copy(panEnd);

                update();
                break;

            default:

                state = CONTROL_STATE.NONE;

        }
    };

    var onTouchEnd = function (/* event */) {
        if (this.enabled === false) return;

        this.dispatchEvent(this.endEvent);
        state = CONTROL_STATE.NONE;
        keyCtrl = false;
        keyShift = false;
        keyS = false;
    };

    this.stopAnimation = function () {
        lastRotation.splice(0);
        ctrl.progress = 0;
        updateGlobeTarget.bind(this)();
    };

    // update object camera position
    this.updateCameraTransformation = function (controlState)
    {
        const d = this.enableDamping;
        this.enableDamping = false;
        state = controlState || CONTROL_STATE.ORBIT;
        update();
        updateGlobeTarget.bind(this)();
        this.enableDamping = d;
    };

    this.dispose = function () {
        // this.domElement.removeEventListener( 'contextmenu', onContextMenu, false );
        this.domElement.removeEventListener('mousedown', onMouseDown, false);
        this.domElement.removeEventListener('mousewheel', onMouseWheel, false);
        this.domElement.removeEventListener('DOMMouseScroll', onMouseWheel, false); // firefox

        this.domElement.removeEventListener('touchstart', onTouchStart, false);
        this.domElement.removeEventListener('touchend', onTouchEnd, false);
        this.domElement.removeEventListener('touchmove', onTouchMove, false);

        this.domElement.removeEventListener('mousemove', onMouseMove, false);
        this.domElement.removeEventListener('mouseup', onMouseUp, false);

        window.removeEventListener('keydown', onKeyDown, false);

        // this.dispatchEvent( { type: 'dispose' } ); // should this be added here?
    };

    // Instance all
    this.domElement.addEventListener('contextmenu', (event) => {
        event.preventDefault();
    }, false);
    this.domElement.addEventListener('mousedown', onMouseDown.bind(this), false);
    this.domElement.addEventListener('mousewheel', onMouseWheel.bind(this), false);
    this.domElement.addEventListener('dblclick', ondblclick.bind(this), false);
    this.domElement.addEventListener('DOMMouseScroll', onMouseWheel.bind(this), false); // firefox

    this.domElement.addEventListener('touchstart', onTouchStart.bind(this), false);
    this.domElement.addEventListener('touchend', onTouchEnd.bind(this), false);
    this.domElement.addEventListener('touchmove', onTouchMove.bind(this), false);

    // player
    this.domElement.addEventListener('frameAnimation', update.bind(this), false);

    // TODO: Why windows
    window.addEventListener('keydown', onKeyDown.bind(this), false);
    window.addEventListener('keyup', onKeyUp.bind(this), false);

    // Initialisation Globe Target and movingGlobeTarget
    var positionTarget = new THREE.Vector3().copy(camera.position).setLength(tSphere.radius);
    setGlobleTarget(positionTarget);
    movingGlobeTarget.copy(positionTarget);
    this.camera.up.copy(positionTarget.normalize());

    update();

    engine.scene3D.add(globeTarget);

    if (enableTargetHelper) {
        globeTarget.add(new THREE.AxisHelper(500000));
        engine.scene3D.add(this.pickingHelper);
    }

    // Start position
    initialTarget = globeTarget.clone();
    initialPosition = this.camera.position.clone();
    initialZoom = this.camera.zoom;

    _handlerMouseMove = onMouseMove.bind(this);
    _handlerMouseUp = onMouseUp.bind(this);
}

GlobeControls.prototype = Object.create(THREE.EventDispatcher.prototype);
GlobeControls.prototype.constructor = GlobeControls;

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// API Function

GlobeControls.prototype.setTilt = function (tilt, isAnimated) {
    const deltaPhi = (tilt * Math.PI / 180 - this.getTiltRad());
    return this.moveOrbitalPosition(0, 0, deltaPhi, isAnimated);
};

GlobeControls.prototype.setHeading = function (heading, isAnimated) {
    const deltaTheta = (heading * Math.PI / 180 - this.getHeadingRad());
    return this.moveOrbitalPosition(0, deltaTheta, 0, isAnimated);
};

GlobeControls.prototype.setRange = function (pRange, isAnimated) {
    const deltaRange = pRange - this.getRange();
    return this.moveOrbitalPosition(deltaRange, 0, 0, isAnimated);
};

const destSpherical = new THREE.Spherical();

GlobeControls.prototype.setOrbitalPosition = function (range, heading, tilt, isAnimated) {
    const deltaPhi = ((tilt || this.getTilt()) * Math.PI / 180 - this.getTiltRad());
    const deltaTheta = ((heading || this.getHeading()) * Math.PI / 180 - this.getHeadingRad());
    const deltaRange = (range || this.getRange()) - this.getRange();
    return this.moveOrbitalPosition(deltaRange, deltaTheta, deltaPhi, isAnimated);
};

GlobeControls.prototype.moveOrbitalPosition = function (deltaRange, deltaTheta, deltaPhi, isAnimated) {
    const range = deltaRange + this.getRange();
    if (isAnimated) {
        destSpherical.theta = deltaTheta + spherical.theta;
        destSpherical.phi = deltaPhi + spherical.phi;
        sphericalTo.radius = range;
        sphericalTo.theta = deltaTheta / (animationOrbit.duration - 1);
        sphericalTo.phi = deltaPhi / (animationOrbit.duration - 1);
        state = CONTROL_STATE.ORBIT;
        return player.startAnimation(animationOrbit).then(() => {
            if (player.isEnded()) {
                this.moveOrbitalPosition(0, destSpherical.theta - spherical.theta, destSpherical.phi - spherical.phi);
            }
            this.stopAnimation();
        });
    }
    else {
        sphericalDelta.theta = deltaTheta;
        sphericalDelta.phi = deltaPhi;
        orbit.scale = range / this.getRange();
        this.updateCameraTransformation();
        return new Promise((r) => { r(); });
    }
};

/**
 * Gets camera's target position
 *
 * @return     {Vecto3}
 */
GlobeControls.prototype.getTargetCameraPosition = function () {
    return globeTarget.position;
};

GlobeControls.prototype.setCenter = function (position, isAnimated) {
    const center = this.getTargetCameraPosition();

    snapShotCamera.shot(this.camera);

    ptScreenClick.x = this.domElement.width / 2;
    ptScreenClick.y = this.domElement.height / 2;

    const vFrom = center.clone().normalize();
    const vTo = position.normalize();

    if (isAnimated) {
        ctrl.qDelta.setFromUnitVectors(vFrom, vTo);
        if (position.range) {
            animatedScale = 1.0 - position.range / this.getRange();
        }
        state = CONTROL_STATE.MOVE_GLOBE;
        return player.startAnimation(animationZoomCenter).then(() => {
            animatedScale = 0.0;
            this.stopAnimation();
        });
    }
    else {
        quatGlobe.setFromUnitVectors(vFrom, vTo);
        this.updateCameraTransformation(CONTROL_STATE.MOVE_GLOBE);
        return new Promise((r) => { r(); });
    }
};

GlobeControls.prototype.getRange = function () {
    return this.getTargetCameraPosition().distanceTo(this.camera.position);
};

// TODO idea : remove this? used in API
GlobeControls.prototype.getRay = function () {
    var direction = new THREE.Vector3(0, 0, 1);
    this.camera.localToWorld(direction);
    direction.sub(this.camera.position).negate().normalize();

    return {
        origin: this.camera.position,
        direction,
    };
};

GlobeControls.prototype.getTilt = function () {
    return spherical.phi * 180 / Math.PI;
};

GlobeControls.prototype.getHeading = function () {
    return spherical.theta * 180 / Math.PI;
};

// Same functions

GlobeControls.prototype.getTiltRad = function () {
    return spherical.phi;
};

GlobeControls.prototype.getHeadingRad = function () {
    return spherical.theta;
};

GlobeControls.prototype.getPolarAngle = function () {
    return spherical.phi;
};

GlobeControls.prototype.getAzimuthalAngle = function () {
    return spherical.theta;
};

GlobeControls.prototype.moveTarget = function () {
    return movingGlobeTarget;
};

GlobeControls.prototype.pan = function (deltaX, deltaY) {
    this.mouseToPan(deltaX, deltaY);
    this.updateCameraTransformation(CONTROL_STATE.PAN);
};

// End API functions
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

GlobeControls.prototype.reset = function () {
    // TODO not reset target globe

    state = CONTROL_STATE.NONE;

    this.target.copy(initialTarget);
    this.camera.position.copy(initialPosition);
    this.camera.zoom = initialZoom;

    this.camera.updateProjectionMatrix();
    this.dispatchEvent(this.changeEvent);

    this.updateCameraTransformation();
};

export default GlobeControls;
