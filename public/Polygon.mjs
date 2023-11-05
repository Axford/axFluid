
import BoundingBox from "./BoundingBox.mjs";
import Line from "./Line.mjs";
import Vector from "./Vector.mjs";

export default class Polygon {
    constructor() {
        this.edges = [];
        this.bb = new BoundingBox();
    }

    reset() {
        this.edges = [];
        this.bb.reset();
    }

    addEdge(e) {
        this.edges.push(e);
        this.bb.addPoint(e.p0);
        this.bb.addPoint(e.p1);
    }

    finalise() {
        this.bb.finalise();
    }

    intersectsEdge(edge) {
        var res = [];  // array of intersections

        this.edges.forEach((e)=>{
            // check for intersection
            var i = e.intersects(edge);
            if (i.intersects) {
                res.push(i);
            }
        });

        return res;
    }

    contains(p) {
        // calculate crossing number
        var c = 0;

        // construct ray from p to outside upper bound
        var ray = new Line(p, new Vector(p.x, this.bb.topLeft.y));

        this.edges.forEach((e)=>{
            // check for intersection
            var i = ray.intersects(e);
            if (i.intersects) {
                c++;
            }
        });

        return (c % 2) == 1;
    }

}