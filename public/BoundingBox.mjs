
import Vector from "./Vector.mjs";
import Line from "./Line.mjs";

export default class BoundingBox {
    constructor() {
        this.topLeft = new Vector(0,0);
        this.topRight = new Vector(0,0);
        this.bottomRight = new Vector(0,0);
        this.bottomLeft = new Vector(0,0);

        this.reset();
    }

    reset() {
        this.initialised = false;
        this.finalised = false;
    }

    addPoint(p) {
        if (this.initialised) {
            if (p.x < this.topLeft.x) this.topLeft.x = p.x;
            if (p.y < this.topLeft.y) this.topLeft.y = p.y;

            if (p.x > this.bottomRight.x) this.bottomRight.x = p.x;
            if (p.y > this.bottomRight.y) this.bottomRight.y = p.y;

        } else {
            this.topLeft.set(p);
            this.bottomRight.set(p);
        }
        this.initialised = true;
    }

    finalise() {
        // compute all corners
        this.bottomLeft.x = this.topLeft.x;
        this.bottomLeft.y = this.bottomRight.y;
        this.topRight.y = this.topLeft.y;
        this.topRight.x = this.bottomRight.x;

        // ... and Lines for edges
        this.left = new Line(this.bottomLeft, this.topLeft);
        this.top = new Line(this.topLeft, this.topRight);
        this.right = new Line(this.topRight, this.bottomRight);
        this.bottom = new Line(this.bottomRight, this.bottomLeft);
    }

    containsPoint(p) {
        if (!this.initialised) return false;

        return (
            p.x >= this.topLeft.x && 
                p.y >= this.topLeft.y &&
                p.x <= this.bottomRight.x &&
                p.y <= this.bottomRight.y
               );
    }

    containsCoords(x,y) {
        if (!this.initialised) return false;

        return (
            x >= this.topLeft.x && 
                y >= this.topLeft.y &&
                x <= this.bottomRight.x &&
                y <= this.bottomRight.y
               );
    }

    overlaps(b) {
        // compare with another bounding box
        // check each corner to b to see if it's within us
        return  this.containsCoords(b.topLeft.x, b.topLeft.y) ||
                this.containsCoords(b.topLeft.x, b.bottomRight.y) ||
                this.containsCoords(b.bottomRight.x, b.topLeft.y) ||
                this.containsCoords(b.bottomRight.x, b.bottomRight.y);
    }
}