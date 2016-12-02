/**
 * Generated On: 2015-10-5
 * Class: Layer
 * Description: Le layer est une couche de données. Cette couche peut etre des images ou de l'information 3D. Les requètes de cette couche sont acheminées par une interfaceCommander.
 *
 */

/**
 *
 * @param {type} Node
 * @param {type} InterfaceCommander
 * @param {type} Projection
 * @param {type} NodeMesh
 * @returns {Layer_L15.Layer}
 */
import * as THREE from 'three';
import Node from 'Scene/Node';
import Projection from 'Core/Geographic/Projection';
import NodeMesh from 'Renderer/NodeMesh';

function Layer() {
    // Constructor

    Node.call(this);
    this.descriManager = null;
    this.projection = new Projection();
    this.id = Layer.count++;
    this.name = '';
}

Layer.count = 0;

Layer.prototype = Object.create(Node.prototype);

Layer.prototype.constructor = Layer;

Layer.prototype.setName = function (name) {
    this.name = name;
};

Layer.prototype.getName = function () {
    return this.name;
};

Layer.prototype.removeFeature = function (id) {
    this.children[0].children.splice(id, 1);
};


// Should be plural as it return an array of meshes
Layer.prototype.getMesh = function () {
    var meshs = [];

    for (var i = 0; i < this.children.length; i++) {
        var node = this.children[i];


        if (node instanceof NodeMesh || node instanceof THREE.Mesh || node instanceof THREE.Object3D)
            { meshs.push(node); }
        else if (node instanceof Layer) {
            meshs = meshs.concat(node.getMesh());
        }
    }

    return meshs;
};

export default Layer;
