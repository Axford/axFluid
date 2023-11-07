

export default class Line {
    constructor(p0, p1) {
        this.p0 = p0.clone();
        this.p1 = p1.clone();
    }

    intersects(line) {
        // test for intersection with line
        var  res = {
            intersects: false,
            x:0,
            y:0
        }
        const a = this.p0.x,
              b = this.p0.y,
              c = this.p1.x,
              d = this.p1.y;
        const p = line.p0.x,
              q = line.p0.y,
              r = line.p1.x,
              s = line.p1.y;

        var det, gamma, lambda;
        det = (c - a) * (s - q) - (r - p) * (d - b);
        if (det === 0) {
            // parallel or coincident, can't calculate intersection
            
        } else {
            lambda = ((s - q) * (r - a) + (p - r) * (s - b)) / det;
            gamma = ((b - d) * (r - a) + (c - a) * (s - b)) / det;
            if ( (0 <= lambda && lambda <= 1) && (0 <= gamma && gamma <= 1) ) {
                res.intersects = true;
                // calculate point of intersection
                res.x = a + lambda * (c-a);
                res.y = b + lambda * (d-b);
            }
        }
        return res;
    }

    whichSide(p) {
        return (this.p1.x - this.p0.x)*(p.y - this.p0.y) - (this.p1.y - this.p0.y)*(p.x - this.p0.x);
    }

}