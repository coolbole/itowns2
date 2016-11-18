/**
 * Generated On: 2016-09-28
 * Class: FeatureToolBox
 * Description:
 */

import * as THREE from 'THREE';
import CVML from 'Core/Math/CVML';
import BoundingBox from 'Scene/BoundingBox';
import Ellipsoid from 'Core/Math/Ellipsoid';
import GeoCoordinate, {UNIT} from 'Core/Geographic/GeoCoordinate';

function FeatureToolBox() {
    this.size       = {x:6378137,y: 6356752.3142451793,z:6378137};
    this.ellipsoid  = new Ellipsoid(this.size);
    this.arrPolygons = null;
    this.arrLines = null;
}

FeatureToolBox.prototype.GeoJSON2Polygon = function(features) {
    var polyGroup = new THREE.Object3D();
    for (var r = 0; r < features.length; r++) {
        var positions = [];
        //var hauteur = (features[r].properties.hauteur) || 0;
        var polygon = features[r].geometry.coordinates[0][0];
        var altitude = features[r].properties.z_min;
        if (polygon.length > 2 && altitude != 9999) {
            for (var j = 0; j < polygon.length; ++j) {
                var pt2DTab = polygon[j]; //.split(' ');
                //long et puis lat
                //var pt = new THREE.Vector3(parseFloat(pt2DTab[1]), hauteur, parseFloat(pt2DTab[0]));
                var geoCoord = new GeoCoordinate(parseFloat(pt2DTab[1]), parseFloat(pt2DTab[0]), altitude, UNIT.DEGREE)
                var spt = this.tool.ellipsoid.cartographicToCartesian(geoCoord);
                positions.push( spt.x, spt.y, spt.z);
            }
            var geometry = new THREE.BufferGeometry();
            var material = new THREE.LineBasicMaterial({ color: 0xff0000, transparent : true, opacity: 0.9}); //side:THREE.DoubleSide, , linewidth: 5,
            geometry.addAttribute( 'position', new THREE.BufferAttribute( new Float32Array( positions ), 3 ) );
            geometry.computeBoundingSphere();
            var poly = new THREE.Line( geometry, material );
            poly.frustumCulled = false;
            polyGroup.add(poly);
        }
    }

    return polyGroup;
};


FeatureToolBox.prototype.GeoJSON2Box = function(features) {
    var bboxGroup = new THREE.Object3D();
    var wallGeometry = new THREE.Geometry(); // for the walls
    var roofGeometry = new THREE.Geometry(); // for the roof
    var suppHeight = 10; // So we don't cut the roof
    //var texture = new THREE.TextureLoader().load( 'data/strokes/wall-texture.jpg');

    for (var r = 0; r < features.length; r++) {
        var hauteur = (features[r].properties.hauteur + suppHeight) || 0;
        var altitude = features[r].properties.z_min;
        var polygon = features[r].geometry.coordinates[0][0];
        var goodAltitude;

        if (polygon.length > 2) {

            if(altitude != 9999) goodAltitude = altitude;

            var arrPoint2D = [];
            // VERTICES
            for (var j = 0; j < polygon.length - 1; ++j) {
                var pt2DTab = polygon[j]; //.split(' ');

                var geoCoord1 = new GeoCoordinate(parseFloat(pt2DTab[1]), parseFloat(pt2DTab[0]), goodAltitude, UNIT.DEGREE);
                var geoCoord2 = new GeoCoordinate(parseFloat(pt2DTab[1]), parseFloat(pt2DTab[0]), goodAltitude + hauteur, UNIT.DEGREE);
                var pgeo1 = this.ellipsoid.cartographicToCartesian(geoCoord1);
                var pgeo2 = this.ellipsoid.cartographicToCartesian(geoCoord2);

                var vector3_1 = new THREE.Vector3(pgeo1.x, pgeo1.y, pgeo1.z);
                var vector3_2 = new THREE.Vector3(pgeo2.x, pgeo2.y, pgeo2.z);

                arrPoint2D.push(CVML.newPoint(parseFloat(pt2DTab[1]), parseFloat(pt2DTab[0])));
                //arrPoint2D.push(CVML.newPoint(pgeo2.x, pgeo2.y, pgeo2.z));

                wallGeometry.vertices.push(vector3_1, vector3_2);

            }

            // FACES
            // indice of the first point of the polygon 3D
            for (var k = wallGeometry.vertices.length - ((polygon.length - 1) * 2); k < wallGeometry.vertices.length; k = k + 2) {

                var l = k; // % (pts2DTab.length);
                if (l > wallGeometry.vertices.length - 4) {
                    l = wallGeometry.vertices.length - ((polygon.length - 1) * 2);
                }
                wallGeometry.faces.push(new THREE.Face3(l, l + 1, l + 3));
                wallGeometry.faces.push(new THREE.Face3(l, l + 3, l + 2));
            }

            var ll = wallGeometry.vertices.length - ((polygon.length - 1) * 2);
            wallGeometry.faces.push(new THREE.Face3(ll, ll + 1, wallGeometry.vertices.length - 1));
            wallGeometry.faces.push(new THREE.Face3(ll, wallGeometry.vertices.length - 1, wallGeometry.vertices.length - 2));
        }

        wallGeometry.computeFaceNormals(); // WARNING : VERY IMPORTANT WHILE WORKING WITH RAY CASTING ON CUSTOM MESH

        //**************** ROOF ****************************

        var triangles = CVML.TriangulatePoly(arrPoint2D);
        triangles.forEach(function(t) {

            var pt1 = t.getPoint(0),
                pt2 = t.getPoint(1),
                pt3 = t.getPoint(2);

            var geoCoord1 = new GeoCoordinate(pt1.x, pt1.y, goodAltitude + hauteur, UNIT.DEGREE);
            var geoCoord2 = new GeoCoordinate(pt2.x, pt2.y, goodAltitude + hauteur, UNIT.DEGREE);
            var geoCoord3 = new GeoCoordinate(pt3.x, pt3.y, goodAltitude + hauteur, UNIT.DEGREE);

            var pgeo1 = this.ellipsoid.cartographicToCartesian(geoCoord1); //{longitude:p1.z, latitude:p1.x, altitude: 0});
            var pgeo2 = this.ellipsoid.cartographicToCartesian(geoCoord2);
            var pgeo3 = this.ellipsoid.cartographicToCartesian(geoCoord3);

            //var geometry = new THREE.Geometry();
            roofGeometry.vertices.push(new THREE.Vector3(pgeo1.x, pgeo1.y, pgeo1.z));
            roofGeometry.vertices.push(new THREE.Vector3(pgeo2.x, pgeo2.y, pgeo2.z));
            roofGeometry.vertices.push(new THREE.Vector3(pgeo3.x, pgeo3.y, pgeo3.z));

            var face = new THREE.Face3(
                roofGeometry.vertices.length - 3,
                roofGeometry.vertices.length - 2,
                roofGeometry.vertices.length - 1
            );
            roofGeometry.faces.push(face);

        }.bind(this));

    }

    roofGeometry.computeFaceNormals();

    var wallMat = new THREE.MeshBasicMaterial({color: 0xcccccc, transparent: true, opacity: 0.8, side : THREE.DoubleSide});  // map : texture,
    var roofMat = new THREE.MeshBasicMaterial({color: 0x660000, transparent: true, opacity: 0.8, side : THREE.DoubleSide});

    var wall  = new THREE.Mesh(wallGeometry, wallMat);
        wall.frustumCulled = false;
    var roof  = new THREE.Mesh(roofGeometry, roofMat);
        roof.frustumCulled = false;

    bboxGroup.add(wall);
    bboxGroup.add(roof);

    return bboxGroup;
};

/**
 * Permit to cut a line at the end of the tile.
 * @param coords: the coords which are tested
 * @param slope: the slope of the current portion of line
 * @param rest:
 * @param bbox; the tile bounding box
 */
FeatureToolBox.prototype.cutLine = function(coords, slope, rest, bbox) {
    var minLong = bbox.west(), maxLong = bbox.east(), minLat  = bbox.south(), maxLat  = bbox.north();
    if(coords[0] < minLong){
        coords[0] = minLong;
        if(coords[1] >= minLat && coords[1] <= maxLat)
            coords[1] = slope * coords[0] + rest;
    }
    else if (coords[0] > maxLong){
        coords[0] = maxLong;
        if(coords[1] >= minLat && coords[1] <= maxLat)
            coords[1] = slope * coords[0] + rest;
    }
    if(coords[1] < minLat){
        coords[1] = minLat;
        if(coords[0] >= minLong && coords[0] <= maxLong)
            coords[0] = (coords[1] - rest) / slope;
    }
    else if (coords[1] > maxLat){
        coords[1] = maxLat;
        if(coords[0] >= minLong && coords[0] <= maxLong)
            coords[0] = (coords[1] - rest) / slope;
    }
};

/**
 * From a single point, the direction of the line and the orientation of the tile
 * compute the two points which will be on the border of the line.
 * @param pt1: one point of the line
 * @param pt2: another point of the line
 * @param isFirstPt: permit to choose to which point we will compute the border points
 * @param offsetValue: Half value of the line size
 */
FeatureToolBox.prototype.computeLineBorderPoints = function(pt1, pt2, isFirstPt, offsetValue) {
    var geoCoord1 = new GeoCoordinate(pt1.x, pt1.y, pt1.z, UNIT.DEGREE);
    var geoCoord2 = new GeoCoordinate(pt2.x, pt2.y, pt2.z, UNIT.DEGREE);

    var cart1 = this.ellipsoid.cartographicToCartesian(geoCoord1);
    var cart2 = this.ellipsoid.cartographicToCartesian(geoCoord2);

    var dx      = cart2.x - cart1.x;
    var dy      = cart2.y - cart1.y;
    var dz      = cart2.z - cart1.z;

    var direct  = new THREE.Vector3(dx, dy, dz);
    direct.normalize();
    var normalGlobe = this.ellipsoid.geodeticSurfaceNormalCartographic(geoCoord1);
    normalGlobe.normalize();

    normalGlobe.cross(direct);
    normalGlobe.normalize();

    //Compute offset to find the left and right point with the given offset value
    var offsetX = normalGlobe.x * offsetValue;
    var offsetY = normalGlobe.y * offsetValue;
    var offsetZ = normalGlobe.z * offsetValue;

    //The first point left and point right of the line
    var left, right;
    if(isFirstPt){
        left    = new THREE.Vector3(cart1.x - offsetX, cart1.y - offsetY, cart1.z - offsetZ);
        right   = new THREE.Vector3(cart1.x + offsetX, cart1.y + offsetY, cart1.z + offsetZ);
    } else {
        left    = new THREE.Vector3(cart2.x - offsetX, cart2.y - offsetY, cart2.z - offsetZ);
        right   = new THREE.Vector3(cart2.x + offsetX, cart2.y + offsetY, cart2.z + offsetZ);
    }
    return {left: left, right: right};
};

/**
 * Process the data received from a WFS request with a tile of feature type 'Line'.
 * Can be used whe the type of the feature tile is a Grid and not a quadTree
 * @param features: the data received as JSON inside a tab
 * @param box: 		the tile bounding box (rad)
 * @param layer: 	the current layer with specific parameters
 */
FeatureToolBox.prototype.GeoJSON2Line = function(features, box, layer) {
    var bbox = new BoundingBox( box.west()  * 180.0 / Math.PI,
								box.east()  * 180.0 / Math.PI,
								box.south() * 180.0 / Math.PI,
								box.north() * 180.0 / Math.PI,
								box.bottom(), box.top());
    var minLong  = bbox.west(), maxLong = bbox.east(), minLat  = bbox.south(),  maxLat  = bbox.north();
    var geometry = new THREE.Geometry();

    for (var i = 0; i < features.length; i++) {
        var feature = features[i];
        var coords  = feature.geometry.coordinates;

        var j = 0;
        var inTile = false;

        //Cut the line according to the tiles limits
        //May not be usefull if we can cut the line before inside the request to the WFS provider
        do{
            var c_1 = coords[j - 1], c = coords[j], cX = c[0], cY = c[1], c1  = coords[j + 1];
            if(c_1 != undefined) var c_1X = c_1[0], c_1Y = c_1[1];
            if(c1 != undefined)  var c1X  = c1[0],  c1Y  = c1[1];

            if (cX < minLong || cX > maxLong || cY < minLat || cY > maxLat) {
                var coeffSlope, rest;
                if(inTile) {
                    coeffSlope = (cY - c_1Y) / (cX - c_1X);
                    rest = cY - coeffSlope * cX;

                    this.cutLine(c, coeffSlope, rest, bbox);
                    j++;
                } else if (c1 != undefined && c1X > minLong && c1X < maxLong &&  c1Y > minLat  && c1Y < maxLat) {
                    coeffSlope = (c1Y - cY) / (c1X - cX);
                    rest = c1Y - coeffSlope * c1X;

                    this.cutLine(c, coeffSlope, rest, bbox);
                    j++;
                } else
                    coords.splice(j, 1);
                inTile = false;
            } else {
                inTile = true;
                j++;
            }
        }while (j < coords.length);

        if(coords.length > 1){
            var resp = this.computeLineBorderPoints(new THREE.Vector3(coords[0][0], coords[0][1], 180),
                                                    new THREE.Vector3(coords[1][0], coords[1][1], 180),
                                                    true, layer.params.length || 10);

            for (j = 0; j < coords.length - 1; j++) {
                var currentGeometry = new THREE.Geometry();
                currentGeometry.vertices.push(resp.left, resp.right);

                resp = this.computeLineBorderPoints(new THREE.Vector3(coords[j][0],     coords[j][1], 180),
                                                    new THREE.Vector3(coords[j + 1][0], coords[j + 1][1], 180),
                                                    false, layer.params.length || 10);

                currentGeometry.vertices.push(resp.left, resp.right);

                currentGeometry.faces.push( new THREE.Face3(0, 2, 1),
                                            new THREE.Face3(2, 3, 1));

                geometry.computeFaceNormals();
                geometry.computeVertexNormals();

                for (var k = 0; k < currentGeometry.faces.length; k++)
                    this.manageColor(feature.properties, currentGeometry.faces[k].color, layer);

                geometry.merge(currentGeometry);
            }
        }
    }
    return geometry;
};

/**
 * Create the entire geometry of the object passed in. Is use to create feature geometry
 * like points or boxes.
 * @param features: the data received as JSON inside a tab
 * @param bbox:  	the tile bounding box (rad)
 * @param type:  	the type of mesh for this layer
 * @param layer: 	the current layer with specific parameters
 * @param node: 	the current node
 */
FeatureToolBox.prototype.GeoJSON2Point = function(features, bbox, type, layer, node) {
    var geometry = new THREE.Geometry();
    for (var i = 0; i < features.length; i++) {
        var feature = features[i];
        var coords = feature.geometry.coordinates;

        var geoCoord = new GeoCoordinate(coords[0], coords[1], ((bbox.bottom() + bbox.top()) / 2) + 3, UNIT.DEGREE);
        var normalGlobe = this.ellipsoid.geodeticSurfaceNormalCartographic(geoCoord);
        var centerPoint = this.ellipsoid.cartographicToCartesian(geoCoord);

        //Change the type of height and radius computation.
        if(layer.params && layer.params.retail) {
            type = layer.params.retail(centerPoint, type, new THREE.Vector3());
            node.retailType = layer.params.getRetailType();
        }

        var currentGeometry;
        var params = layer.params;
        if(type == 'box')
            currentGeometry = new THREE.BoxGeometry(params.boxWidth || 40, params.boxWidth || 40, params.boxHeight || 80);
        else if(type == 'point')
            currentGeometry = new THREE.CircleGeometry(params.radius || 10, params.nbSegment || 3, params.thetaStart || 0, params.thetaLength || 2 * Math.PI);
        else
            continue;

        currentGeometry.lookAt(normalGlobe);
        currentGeometry.translate(centerPoint.x, centerPoint.y, centerPoint.z);

        for (var j = 0; j < currentGeometry.faces.length; j++)
            this.manageColor(feature.properties, currentGeometry.faces[j].color, layer);

        geometry.merge(currentGeometry);
    }
    return geometry;
};

/**
 * Manage to put the colors inside the color manager for a feature type 'Point'.
 * @param properties: properties of the feature
 * @param color: 	  manager of the color of a face
 * @params layer: 	  the current layer with specific parameters
 */
FeatureToolBox.prototype.manageColor = function(properties, color, layer) {
    var colorParams = layer.params.color || undefined;

    if(colorParams !== undefined)
        for (var i = 0; i < colorParams.testTab.length; i++) {
            if(properties[colorParams.property] === colorParams.testTab[i]){
                color.setHex(colorParams.colorTab[i]);
                return;
            }
        }
    color.setHex(new THREE.Color(0xFFFFFF));
};

////addFeature
FeatureToolBox.prototype.createGeometryArray = function(json) {
    var geometry_array = [];

    if (json.type == 'Feature') {
        geometry_array.push(json.geometry);
    } else if (json.type == 'FeatureCollection') {
        for (var feature_num = 0; feature_num < json.features.length; feature_num++) {
            geometry_array.push(json.features[feature_num].geometry);
        }
    } else if (json.type == 'GeometryCollection') {
        for (var geom_num = 0; geom_num < json.geometries.length; geom_num++) {
            geometry_array.push(json.geometries[geom_num]);
        }
    } else {
        throw new Error('The geoJSON is not valid.');
    }
    //alert(geometry_array.length);
    return geometry_array;
};

FeatureToolBox.prototype.convertLonLatToWGS84 = function(coordinates_array) {
    var lon = coordinates_array[0];
    var lat = coordinates_array[1];
    var geoCoord = new GeoCoordinate(lon, lat, 4180, UNIT.DEGREE);
    return this.ellipsoid.cartographicToCartesian(geoCoord);
};

FeatureToolBox.prototype.getMidpoint = function(point1, point2) {
    var midpoint_lon = (point1[0] + point2[0]) / 2;
    var midpoint_lat = (point1[1] + point2[1]) / 2;
    var midpoint = [midpoint_lon, midpoint_lat];

    return midpoint;
}

FeatureToolBox.prototype.needsInterpolation = function(point2, point1) {
    //If the distance between two latitude and longitude values is
    //greater than five degrees, return true.
    var lon1 = point1[0];
    var lat1 = point1[1];
    var lon2 = point2[0];
    var lat2 = point2[1];
    var lon_distance = Math.abs(lon1 - lon2);
    var lat_distance = Math.abs(lat1 - lat2);

    if (lon_distance > 5 || lat_distance > 5) {
        return true;
    } else {
        return false;
    }
};


FeatureToolBox.prototype.interpolatePoints = function(interpolation_array) {
    //This function is recursive. It will continue to add midpoints to the
    //interpolation array until needsInterpolation() returns false.
    var temp_array = [];
    var point1, point2;

    for (var point_num = 0; point_num < interpolation_array.length-1; point_num++) {
        point1 = interpolation_array[point_num];
        point2 = interpolation_array[point_num + 1];

        if (this.needsInterpolation(point2, point1)) {
            temp_array.push(point1);
            temp_array.push(this.getMidpoint(point1, point2));
        } else {
            temp_array.push(point1);
        }
    }

    temp_array.push(interpolation_array[interpolation_array.length-1]);

    if (temp_array.length > interpolation_array.length) {
        temp_array = this.interpolatePoints(temp_array);
    } else {
        return temp_array;
    }
    return temp_array;
};

FeatureToolBox.prototype.createCoordinateArray = function(feature) {
    //Loop through the coordinates and figure out if the points need interpolation.
    var temp_array = [];
    var interpolation_array = [];

        for (var point_num = 0; point_num < feature.length; point_num++) {
            var point1 = feature[point_num];
            var point2 = feature[point_num - 1];

            if (point_num > 0) {
                if (this.needsInterpolation(point2, point1)) {
                    interpolation_array = [point2, point1];
                    interpolation_array = this.interpolatePoints(interpolation_array);

                    for (var inter_point_num = 0; inter_point_num < interpolation_array.length; inter_point_num++) {
                        temp_array.push(interpolation_array[inter_point_num]);
                    }
                } else {
                    temp_array.push(point1);
                }
            } else {
                temp_array.push(point1);
            }
        }
    return temp_array;
};

FeatureToolBox.prototype.intersectsegment = function( a, b, i, p){

   var d = new THREE.Vector2();
   var e = new THREE.Vector2();
   d.x = b.x - a.x;
   d.y = b.y - a.y;
   e.x = p.x - i.x;
   e.y = p.y - i.y;
   var denom = d.x*e.y - d.y*e.x;
   if (denom == 0.)
       return -1;   // erreur, cas limite
   var t = - (a.x*e.y-i.x*e.y-e.x*a.y+e.x*i.y) / denom;
   if (t < 0. || t >= 1.)
      return 0;
   var u = - (-d.x*a.y+d.x*i.y+d.y*a.x-d.y*i.x) / denom;
   if (u < 0. || u >= 1.)
      return 0;
   return 1;
};

/*
 * Function that tests if a point p is inside a polygon using the classic
 * Ray Casting algorithm
 * @param {type} posGeo
 * @returns {undefined}
 */
FeatureToolBox.prototype.inPolygon = function(p, arrPoints){

    var k = new THREE.Vector2(6.,75.); //This point should be out of the polygon
    var nbintersections = 0;
    for(var i=0; i < arrPoints.length -1; i++){

       var a = arrPoints[i];
       var b = arrPoints[i+1];//.xy;
       var iseg = this.intersectsegment(a,b,k,p);
       nbintersections += iseg;
    }
    return ((nbintersections % 2) === 1);
};

/**
 *
 * display feature attribute information at position p (lonlat)
 */
FeatureToolBox.prototype.showFeatureAttributesAtPos = function(p){

    var intersect = false;
    var desc = "";
    var i=0;
    while(!intersect && i< this.arrPolygons.length){
        intersect = this.inPolygon(p, this.arrPolygons[i].polygon);
        i++;
    }

    i--;
    if(intersect && this.arrPolygons[i].properties.description !== undefined){
        desc = this.arrPolygons[i].properties.description;
    }

    if(!intersect) desc = "noIntersect";

    return desc;
}


FeatureToolBox.prototype.drawLine = function(coordOrigin, tileWH, p1, p2, thickness, ctx, prop) {

    var tilePx = 256;
    ctx.strokeStyle = prop.stroke; //"rgba(255, 0, 255, 0.5)";//"rgba(1,1,0,1)";
    ctx.lineWidth = prop["stroke-width"];
    ctx.globalAlpha = prop["stroke-opacity"];
    ctx.beginPath();
    var a = p1.sub(coordOrigin);
    var b = p2.sub(coordOrigin);
    a.divideScalar(tileWH.x).multiplyScalar(tilePx);
    b.divideScalar(tileWH.x).multiplyScalar(tilePx);
   //      console.log("aa",coordOrigin,tileWH,p1,p2,a,b);
      //       if( (a.x>=0 && a.x<= tilePx && a.y >=0 && a.y <=tilePx)
      //          && (b.x>=0 && b.x<= tilePx && b.y >=0 && b.y <=tilePx)){
                ctx.moveTo(a.x, tilePx - a.y);
                ctx.lineTo(b.x, tilePx - b.y);
                ctx.stroke();
                ctx.globalAlpha = 1.;
      //      }
   // }

};

/*
 * Draw 2D polygon in tile (rasterized)
 * @param {type} polygon
 * @param {type} coordOrigin
 * @param {type} tileWH
 * @param {type} ctx
 * @param {type} prop
 * @returns {undefined}
 */
FeatureToolBox.prototype.drawPolygon = function(polygon, coordOrigin, tileWH, ctx, prop) {

    var tilePx = 256;
    ctx.strokeStyle = prop.stroke;
    ctx.lineWidth = prop["stroke-width"];
    ctx.fillStyle = prop.fill; //"rgba(255, 0, 255, 0.5)";//"rgba(1,1,0,1)";
    ctx.globalAlpha = prop["fill-opacity"];
    // ctx["stroke-opacity"] = 0.8;
    // console.log(ctx.strokeStyle);
    ctx.beginPath();

    for (var i = 0; i< polygon.length -1; ++i){

        var p1 = polygon[i];
        var p2 = polygon[i+1];
        var a = p1.clone().sub(coordOrigin);
        var b = p2.clone().sub(coordOrigin);

        a.divideScalar(tileWH.x).multiplyScalar(tilePx);
        b.divideScalar(tileWH.x).multiplyScalar(tilePx);

        if(i === 0) {
            ctx.moveTo(a.x, tilePx - a.y);
            ctx.lineTo(b.x, tilePx - b.y);
        }else{
            ctx.lineTo(b.x, tilePx - b.y);
        }
     }
     ctx.closePath();
     ctx.globalAlpha = prop["fill-opacity"];
     ctx.fill();
     ctx.globalAlpha = prop["stroke-opacity"];
     ctx.stroke();
     ctx.globalAlpha = 1.;
};


// parameters in deg, vec2
FeatureToolBox.prototype.createRasterImage = function(coordOrigin, tileWH, lines, polygons){

    var c = document.createElement('canvas');
    c.width = 256;
    c.height = 256;
    var ctx = c.getContext("2d");
    // Lines
    for (var j= 0; j< lines.length; ++j){
        var line = lines[j].line;
        var properties = lines[j].properties;
        for (var i= 0; i< line.length -1; ++i){
           this.drawLine(coordOrigin, tileWH, line[i].clone(), line[i+1].clone(), 4, ctx, properties);
        }
    }
    // Polygon
    for(i= 0; i< polygons.length; ++i){
        this.drawPolygon(polygons[i].polygon, coordOrigin, tileWH, ctx, polygons[i].properties);
    }

    var texture = new THREE.Texture(c);
    texture.flipY = true;  // FALSE by default on THREE.DataTexture but True by default for THREE.Texture!
    texture.needsUpdate = true;
    texture.name = "featureRaster";
    return texture;
};




 // Extract polygons and lines for raster rendering in GPU
FeatureToolBox.prototype.extractFeatures = function(json) {

    var arrPolygons = [];
    var arrLines = [];

    for (var nFeature = 0; nFeature < json.features.length; nFeature++) {
        var feat = json.features[nFeature];
        if (feat.geometry.type === 'LineString') {
             var arrLine = [];
             for (var point_num = 0; point_num < feat.geometry.coordinates.length; point_num++) {
                 var v = feat.geometry.coordinates[point_num];
                 arrLine.push(new THREE.Vector2(v[0], v[1]));
             }
             arrLines.push({line:arrLine, properties: feat.properties});
         }

         if (feat.geometry.type === 'Polygon') {
            var arrPolygon = [];
            for (point_num = 0; point_num < feat.geometry.coordinates.length; point_num++) {
                for( var p = 0; p < feat.geometry.coordinates[point_num].length; ++p ){
                   v = feat.geometry.coordinates[point_num][p];
                   arrPolygon.push(new THREE.Vector2(v[0], v[1]));
                }
             }
             arrPolygons.push({polygon:arrPolygon, properties: feat.properties});
         }
    }
    this.arrPolygons = arrPolygons;//console.log(this.arrPolygons);
    this.arrLines = arrLines;
    return {lines: arrLines, polygons: arrPolygons};
};


FeatureToolBox.prototype.createFeaturesPoints = function(json){

    var globalObject = new THREE.Object3D();

    for (var nFeature = 0; nFeature < json.features.length; nFeature++) {
           var feat = json.features[nFeature];
           if (feat.geometry.type === 'Point') {
               //console.log(feat);
               var vertex = this.convertLonLatToWGS84(feat.geometry.coordinates);
               if (feat.properties.name != null){
                   //console.log(feat.properties.name);
                   globalObject.add(this.createText(vertex, feat.properties));
               }else{
                   globalObject.add(this.createIcon(vertex));
               }
            }
       }
     //  globalObject.add(new THREE.LineSegments(geometry, basicLineMaterial));
       return globalObject;
};



FeatureToolBox.prototype.processingGeoJSON = function(json) {

    //console.log(json);
    var jsonFeatures = this.createGeometryArray(json);
    //console.log(jsonFeatures);
    var coordinate_array = [];
    var globalObject = new THREE.Object3D();
    var geometry = new THREE.Geometry();
    var basicLineMaterial = new THREE.LineBasicMaterial({color : 0xff0000});
    for (var nFeature = 0; nFeature < jsonFeatures.length; nFeature++) {
	var point_num, segment_num, vertex;
        if (jsonFeatures[nFeature].type == 'Point') {
            vertex = this.convertLonLatToWGS84(jsonFeatures[nFeature].coordinates);
            globalObject.add(this.createIcon(vertex));
           // geometry.vertices.push(vertex);
           // geometry.vertices.push(vertex);
        }
        else if (jsonFeatures[nFeature].type == 'MultiPoint') {
            for (point_num = 0; point_num < jsonFeatures[nFeature].coordinates.length; point_num++) {
                vertex = this.convertLonLatToWGS84(jsonFeatures[nFeature].coordinates[point_num]);
                geometry.vertices.push(vertex);
                if(point_num>0 && point_num<coordinate_array.length-1) geometry.vertices.push(vertex);
            }

        }
        else if (jsonFeatures[nFeature].type == 'LineString') {
            coordinate_array = this.createCoordinateArray(jsonFeatures[nFeature].coordinates);

            for (point_num = 0; point_num < coordinate_array.length; point_num++) {
                vertex = this.convertLonLatToWGS84(coordinate_array[point_num]);
                geometry.vertices.push(vertex);
                if(point_num>0 && point_num<coordinate_array.length-1) geometry.vertices.push(vertex);
            }

        }

        else if (jsonFeatures[nFeature].type == 'Polygon') {
            for (segment_num = 0; segment_num < jsonFeatures[nFeature].coordinates.length; segment_num++) {
                coordinate_array = this.createCoordinateArray(jsonFeatures[nFeature].coordinates[segment_num]);

                for (point_num = 0; point_num < coordinate_array.length; point_num++) {
                    vertex = this.convertLonLatToWGS84(coordinate_array[point_num]);
                    geometry.vertices.push(vertex);
                    if(point_num>0 && point_num<coordinate_array.length-1) geometry.vertices.push(vertex);
                }
            }
        }
        else if (jsonFeatures[nFeature].type == 'MultiLineString') {
            for (segment_num = 0; segment_num < jsonFeatures[nFeature].coordinates.length; segment_num++) {
                coordinate_array = this.createCoordinateArray(jsonFeatures[nFeature].coordinates[segment_num]);

                for (point_num = 0; point_num < coordinate_array.length; point_num++) {
                    vertex = this.convertLonLatToWGS84(coordinate_array[point_num]);
                    geometry.vertices.push(vertex);
                    if(point_num>0 && point_num<coordinate_array.length-1) geometry.vertices.push(vertex);
                }
            }
        }
        else if (jsonFeatures[nFeature].type == 'MultiPolygon') {
            for (var polygon_num = 0; polygon_num < jsonFeatures[nFeature].coordinates.length; polygon_num++) {
                for (segment_num = 0; segment_num < jsonFeatures[nFeature].coordinates[polygon_num].length; segment_num++) {
                    coordinate_array = this.createCoordinateArray(jsonFeatures[nFeature].coordinates[polygon_num][segment_num]);

                    for (point_num = 0; point_num < coordinate_array.length; point_num++) {
                        vertex = this.convertLonLatToWGS84(coordinate_array[point_num]);
                        geometry.vertices.push(vertex);
                        if(point_num>0 && point_num<coordinate_array.length-1) geometry.vertices.push(vertex);
                    }
                }
            }
        } else {
            throw new Error('The geoJSON is not valid.');
        }
    }


   // if(!bpoint){
        globalObject.add(new THREE.LineSegments(geometry, basicLineMaterial));
        return globalObject;//new THREE.LineSegments(geometry, material);
        //return new THREE.Mesh(new THREE.SphereGeometry(7500000,10,10),new THREE.MeshBasicMaterial());
 /*   }else{
        material = new THREE.PointsMaterial({color : 0xff0000, size : 100});
        return new THREE.Points(geometry, material);
    }
*/
};

/*
 * Create icon at position for point geometry
 * Icon size should be constant. -> need a specific shader
 * @returns {FeatureToolBox.prototype.createIcon.texture|THREE.Texture}
 */
FeatureToolBox.prototype.createIcon = function(pos){

    var image = document.createElement( 'img' );
    var texture = new THREE.Texture( image );
    image.onload = function()  {
        texture.needsUpdate = true;
    };
    image.src = 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QAAAAAAAD5Q7t/AAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4AsVDRIojx5s5AAABPVJREFUWMOtlmtsk1UYx3+nb7uVdu0ulg2oso2xjXF1kzmMsE3GRQSBRRMxJvhBY0xITKYoKHzx8sEoihIjGoyJBiKoYMQwQwK4DVHwE1nGpbKxdVvoZrbR0W20b9v3+KGdtHTbu1X+ycmbN+fy/5/ncp5HSCkZDxK4NXSbyy3NtLW1UbyowtHe0b6ipDC3qtBpW2IW/izQDCEx7Za7T716saX199y8gt86W5s7rBYzS8sfITPDhkGMS4HQEzDoG8HlcqUrqWkvz7YPbc0ePjOf7nPQdxmp9oPUECY7ZBaCs4KhzFW97f7ZP3oH+vcWF81tc2Rl/D8B3mGtVAkN7LO7v1pOy37o74xMKoAhulADwtEN1gyYt4XbxXXukWmFdRnTxE+KIXkBNdpg6yGl4YUc3E1gihJPBA0IApn5aNVfhLh/zesG+CQZAQvxdZzi10059DSDmclDAAHAnA7rfgjiXL0VODwVAVbC/nrqN1TSfjqRXAKh6FdEv8YYl4xCBewPQG1DH7Y5jwEtdxMleCcYDCFhC66vK+kYgzwIKCmQvxat9BUtXPpqiOKnwXJfhDD2PinAYBf8tdsB1PkDqtC1wOCw6rTTd1wcX1nGgCvi99gbOSvwLXz7kq1w7YIEu53eAlePRKwhYmJCpMLGkz6PqaxipsN2ZUILDPvVGtlz9kEGXJGDRhECchZC7XnGJAeoOQwlz0WsFMugBqDjmE1VA4/ruiDNFKwSNxoNEeUxPscA1R/f0A3AlQchLTty89ig9FzAYQk/qivAnhouFjdd8TMhIGcJTF89azJJIHPXRfaMQgF8HViN/jm6AggOW2WgL3FmetGkszCYXuKL2y+AoBfUWxZ9AYoxjFDio3mqkJq4e79EAaFougK6+oLXySiKPK2xN/jnyqT5Td5LaXECwiDS8+jyal26Ai619TSEndUScZcPe5uh+1ifLnugq0m4T8SnrwbMWkZbj/+croBF84tOiRnVraRNv2MFEV15drsj6G1tmlBAw7ZKRrzxGaQImLlGnVc096SugFnZWddk+rzvKNicGMn97ZhObqwc/vto55iur38Crv0SeQFjMyh7EWrOqobsrIyLk60F87l5uZGj5Q5CI4kpaVTAuQIcJWBIAe916D4LI95404++nmu+0WTR1uelph00GAyTroZ7aHrpNS4eGL8Yxf4rY5RqFZjxEGxuOo/RUh2tkTppeAcfUrqzG3tOPNloTJhiRsoY5DJatJbuDmG07BqLXE9AL7Y5eyh7IxLFU30XAkBhLeQ++T1wJqmOCLCiBU5wYm0V7sb44JoIYcCSA09d8JCWWwVcG29pggWklLFjGEPqW1R8MESKJf5xmqiRDAHL3oG03PcmItdzwSj+IPvh9ynfDWGh7woVKK7F53z2Z03T9usdniBACBE3oviIxXVNFGyIEIyHEJAxGyo/77nuvrEjGFTllAWMAz+KeRvLP/Vgn5mYFf91PkZYsU9inlGXbjW7DAaFeyUAoAVb/k6W7wVhjm84RnvFsu2Qt+lL4HBeXi4mk+meCgD4loJnDlC+I77t8gMF66H83T+BXVMr3fFRP5mRKcOBRlm/Xsp9SPkZUh7Kl3LI7ZFSLpjqeckIQEo5V454OuWRxVLuT5HScy4spdyYzFnJCkBKuU72nvfJtiNSSvlmsufovYR6qAOWAC8ydm7o4l8IO+JWe1h9JQAAAABJRU5ErkJggg==';
    var materialS = new THREE.SpriteMaterial( { map: texture,  color: 0xffffff, depthTest: false} );
    var sprite = new THREE.Sprite(materialS);
    sprite.scale.set(500,500,500);
    sprite.position.copy(pos);

    var material = new THREE.LineBasicMaterial({ color: 0xffffff, transparent : false}); //side:THREE.DoubleSide, , linewidth: 5
    var geometry = new THREE.Geometry();
    geometry.vertices.push(new THREE.Vector3(), new THREE.Vector3( -pos.x, -pos.y, -pos.z ).divideScalar(1000));
    var line = new THREE.Line( geometry, material );
    sprite.add(line);
    return sprite;
};


FeatureToolBox.prototype.createText = function(pos, prop){

    var sprite = makeTextSprite(prop.name);
    sprite.position.copy(pos);

    var material = new THREE.LineBasicMaterial({ color: 0xffffff, transparent : false}); //side:THREE.DoubleSide, , linewidth: 5,
    var geometry = new THREE.Geometry();
    geometry.vertices.push(new THREE.Vector3(), new THREE.Vector3( -pos.x, -pos.y, -pos.z ).divideScalar(1000));
    var line = new THREE.Line( geometry, material );
    sprite.add(line);

    return sprite;
}

function makeTextSprite( message, parameters )
{
    if ( parameters === undefined ) parameters = {};
    var fontface = parameters.hasOwnProperty("fontface") ?
            parameters["fontface"] : "Arial";

    var fontsize = parameters.hasOwnProperty("fontsize") ?
            parameters["fontsize"] : 22;

    var borderThickness = parameters.hasOwnProperty("borderThickness") ?
            parameters["borderThickness"] : 4;

    var borderColor = parameters.hasOwnProperty("borderColor") ?
            parameters["borderColor"] : { r:0, g:0, b:0, a:1.0 };

    var backgroundColor = parameters.hasOwnProperty("backgroundColor") ?
            parameters["backgroundColor"] : { r:255, g:255, b:255, a:.8 };

    //var spriteAlignment = THREE.SpriteAlignment.topLeft;

    var canvas = document.createElement('canvas');
    var context = canvas.getContext('2d');
    context.font = "Bold " + fontsize + "px " + fontface;

    // get size data (height depends only on font size)
    var metrics = context.measureText( message );
    var textWidth = metrics.width;

    // background color
    context.fillStyle   = "rgba(" + backgroundColor.r + "," + backgroundColor.g + ","
                             + backgroundColor.b + "," + backgroundColor.a + ")";
    // border color
    context.strokeStyle = "rgba(" + borderColor.r + "," + borderColor.g + ","
                             + borderColor.b + "," + borderColor.a + ")";
    context.lineWidth = borderThickness;
    roundRect(context, borderThickness/2, borderThickness/2, textWidth + borderThickness, fontsize * 1.4 + borderThickness, 6);
    // 1.4 is extra height factor for text below baseline: g,j,p,q.
    // text color
    context.fillStyle = "rgba(0, 0, 0, 1.0)";
    context.fillText( message, borderThickness, fontsize + borderThickness);
    // canvas contents will be used for a texture
    var texture = new THREE.Texture(canvas);
    texture.needsUpdate = true;

    var spriteMaterial = new THREE.SpriteMaterial({ map: texture,  color: 0xffffff, depthTest: false} );//, useScreenCoordinates: false} );
    var sprite = new THREE.Sprite( spriteMaterial );
    sprite.scale.set(1000,1000,1000);//100,50,1.0);
    return sprite;
}

// function for drawing rounded rectangles
function roundRect(ctx, x, y, w, h, r)
{
    ctx.beginPath();
    ctx.moveTo(x+r, y);
    ctx.lineTo(x+w-r, y);
    ctx.quadraticCurveTo(x+w, y, x+w, y+r);
    ctx.lineTo(x+w, y+h-r);
    ctx.quadraticCurveTo(x+w, y+h, x+w-r, y+h);
    ctx.lineTo(x+r, y+h);
    ctx.quadraticCurveTo(x, y+h, x, y+h-r);
    ctx.lineTo(x, y+r);
    ctx.quadraticCurveTo(x, y, x+r, y);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();
}

export default FeatureToolBox;
