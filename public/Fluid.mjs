import TimerStats from "./TimerStats.mjs";
import BoundingBox from "./BoundingBox.mjs";
import Vector from "./Vector.mjs";

export default class Fluid {
  constructor(density, numX, numY, h) {
    this.density = density;
    this.numX = numX + 2;
    this.numY = numY + 2;
    this.numCells = this.numX * this.numY;
    this.h = h;
    this.u = new Float32Array(this.numCells);
    this.v = new Float32Array(this.numCells);
    this.newU = new Float32Array(this.numCells);
    this.newV = new Float32Array(this.numCells);
    this.p = new Float32Array(this.numCells);
    this.s = new Float32Array(this.numCells);
    this.m = new Float32Array(this.numCells);
    this.newM = new Float32Array(this.numCells);
    this.m.fill(1.0);
    var num = numX * numY;

    // compute bounds for each cell
    this.cellBounds = new Array();
    for (var i = 0; i < this.numX; i++) {
      for (var j = 0; j < this.numY; j++) {
        var index = i * this.numY + j;
        var x1 = (i + 0) * this.h;
        var y1 = (j + 0) * this.h;
        var x2 = (i + 1) * this.h;
        var y2 = (j + 1) * this.h;
        var bb = new BoundingBox();
        bb.addPoint(new Vector(x1,y1));
        bb.addPoint(new Vector(x2,y2));
        bb.finalise();
        this.cellBounds.push(bb);
      }
    }

    console.log('cellSize', h);

    this.cellArea = h * h;

    this.U_FIELD = 0;
    this.V_FIELD = 1;
    this.S_FIELD = 2;

    this.scene = null;

    // simulation timers
    this.timers = {};
    this.timers.integration = new TimerStats("integration");
    this.timers.solveIncompressibility = new TimerStats(
      "solveIncompressibility"
    );
    this.timers.extrapolate = new TimerStats("extrapolate");
    this.timers.advectVel = new TimerStats("advectVel");
    this.timers.advectSmoke = new TimerStats("advectSmoke");
    this.timers.total = new TimerStats("total");
  }

  updateObstacle(obstacle) {
    // update "solid" state
    for (var i = 1; i < this.numX - 2; i++) {
      for (var j = 1; j < this.numY - 2; j++) {
        var ci = i * this.numY + j;
        // default all cells to empty
        this.s[ci] = 1.0;

        var cbb = this.cellBounds[ci];

        cbb.li = null;
        cbb.ti = null;
        cbb.ri = null;
        cbb.bi = null;

        // compare to bounds
        if (obstacle.bb.overlaps(cbb)) {

            // there is some overlap, so calculate how much

            // check each corner of bounding box for containment
            cbb.tlc = obstacle.contains(cbb.topLeft);
            cbb.trc = obstacle.contains(cbb.topRight);
            cbb.blc = obstacle.contains(cbb.bottomLeft);
            cbb.brc = obstacle.contains(cbb.bottomRight);

            // solid (s=0) if entirely contained
            this.s[ci] = (cbb.tlc && cbb.trc && cbb.blc && cbb.brc) ? 0 : 1;

            if (this.s[ci] == 0) continue;

            // check for partial containment
            if (cbb.tlc || cbb.trc || cbb.blc || cbb.brc) {

                // check for degenerate cases 
                if (cbb.tlc && cbb.brc && !cbb.blc && !cbb.trc) {
                    // degenerate 1
                } else if (!cbb.tlc && !cbb.brc && cbb.blc && cbb.trc) {
                    // degenerate 2
                } else {

                    var points = [];

                    // starting topLeft, work round corners and edges adding to points if they are inside or intersections

                    if (cbb.tlc) {
                        points.push(cbb.topLeft);
                    }

                    // test top edge
                    var ia = obstacle.intersectsEdge(cbb.top);
                    if (ia.length > 1 ) {
                        cbb.ti = null;
                        // bad - invalid number of intersections
                    } else if (ia.length == 1) {
                        cbb.ti = ia[0];
                        points.push(new Vector(ia[0].x, ia[0].y));
                    } else {
                        cbb.ti = null;
                    }

                    if (cbb.trc) {
                        points.push(cbb.topRight);
                    }

                    // test right edge
                    var ia = obstacle.intersectsEdge(cbb.right);
                    if (ia.length > 1 ) {
                        cbb.ri = null;
                        // bad - invalid number of intersections
                    } else if (ia.length == 1) {
                        cbb.ri = ia[0];
                        points.push(new Vector(ia[0].x, ia[0].y));
                    } else {
                        cbb.ri = null;
                    }

                    if (cbb.brc) {
                        points.push(cbb.bottomRight);
                    }

                    // test bottom edge
                    var ia = obstacle.intersectsEdge(cbb.bottom);
                    if (ia.length > 1 ) {
                        cbb.bi = null;
                        // bad - invalid number of intersections
                    } else if (ia.length == 1) {
                        cbb.bi = ia[0];
                        points.push(new Vector(ia[0].x, ia[0].y));
                    } else {
                        cbb.bi = null;
                    }

                    if (cbb.blc) {
                        points.push(cbb.bottomLeft);
                    }

                    // test left edge
                    var ia = obstacle.intersectsEdge(cbb.left);
                    if (ia.length > 1 ) {
                        cbb.li = null;
                        // bad - invalid number of intersections
                    } else if (ia.length == 1) {
                        cbb.li = ia[0];
                        points.push(new Vector(ia[0].x, ia[0].y));
                    } else {
                        cbb.li = null;
                    }
                    
                    // calc area of points
                    var area = 0;
                    if (points.length > 2) {
                        for (var k=0; k<points.length-1; k++) {
                            area += (points[k].x * points[k+1].y - points[k+1].x * points[k].y);
                        }
                        // add additional area for last point back to first
                        area += (points[points.length-1].x * points[0].y - points[0].x * points[points.length-1].y);
                    }
                    area = 0.5 * Math.abs(area);

                    console.log(points);
                    

                    // set s based on area unobscured
                    var s = 1 - area / this.cellArea;
                    if (s > 1) s = 1;
                    if (s < 0) s = 0;
                    this.s[ci] = s;

                    console.log(i,j,s);
                }
            }


            // check for intersections on left edge of cell

          //this.s[ci] = 0;
        }
      }
    }
  }

  integrate(dt, gravity) {
    var n = this.numY;
    for (var i = 1; i < this.numX; i++) {
      for (var j = 1; j < this.numY - 1; j++) {
        if (this.s[i * n + j] != 0.0 && this.s[i * n + j - 1] != 0.0)
          this.v[i * n + j] += gravity * dt;
      }
    }
  }

  solveIncompressibility(numIters, dt) {
    var n = this.numY;
    var cp = (this.density * this.h) / dt;

    for (var iter = 0; iter < numIters; iter++) {
      for (var i = 1; i < this.numX - 1; i++) {
        for (var j = 1; j < this.numY - 1; j++) {
          // this cell solid if s=0, so skip as nothing to solve
          if (this.s[i * n + j] == 0.0) continue;

          var sx0 = this.s[(i - 1) * n + j];
          var sx1 = this.s[(i + 1) * n + j];
          var sy0 = this.s[i * n + j - 1];
          var sy1 = this.s[i * n + j + 1];
          var s = sx0 + sx1 + sy0 + sy1;
          // if all neighbour cells are solid (s=0), then also nothing to solve
          if (s == 0.0) continue;

          var div =
            this.u[(i + 1) * n + j] -
            this.u[i * n + j] +
            this.v[i * n + j + 1] -
            this.v[i * n + j];

          var p = -div / s;
          p *= this.scene.overRelaxation;
          this.p[i * n + j] += cp * p;

          this.u[i * n + j] -= sx0 * p;
          this.u[(i + 1) * n + j] += sx1 * p;
          this.v[i * n + j] -= sy0 * p;
          this.v[i * n + j + 1] += sy1 * p;
        }
      }
    }
  }

  extrapolate() {
    var n = this.numY;
    for (var i = 0; i < this.numX; i++) {
      this.u[i * n + 0] = this.u[i * n + 1];
      this.u[i * n + this.numY - 1] = this.u[i * n + this.numY - 2];
    }
    for (var j = 0; j < this.numY; j++) {
      this.v[0 * n + j] = this.v[1 * n + j];
      this.v[(this.numX - 1) * n + j] = this.v[(this.numX - 2) * n + j];
    }
  }

  sampleField(x, y, field) {
    var n = this.numY;
    var h = this.h;
    var h1 = 1.0 / h;
    var h2 = 0.5 * h;

    x = Math.max(Math.min(x, this.numX * h), h);
    y = Math.max(Math.min(y, this.numY * h), h);

    var dx = 0.0;
    var dy = 0.0;

    var f;

    switch (field) {
      case this.U_FIELD:
        f = this.u;
        dy = h2;
        break;
      case this.V_FIELD:
        f = this.v;
        dx = h2;
        break;
      case this.S_FIELD:
        f = this.m;
        dx = h2;
        dy = h2;
        break;
    }

    var x0 = Math.min(Math.floor((x - dx) * h1), this.numX - 1);
    var tx = (x - dx - x0 * h) * h1;
    var x1 = Math.min(x0 + 1, this.numX - 1);

    var y0 = Math.min(Math.floor((y - dy) * h1), this.numY - 1);
    var ty = (y - dy - y0 * h) * h1;
    var y1 = Math.min(y0 + 1, this.numY - 1);

    var sx = 1.0 - tx;
    var sy = 1.0 - ty;

    var val =
      sx * sy * f[x0 * n + y0] +
      tx * sy * f[x1 * n + y0] +
      tx * ty * f[x1 * n + y1] +
      sx * ty * f[x0 * n + y1];

    return val;
  }

  avgU(i, j) {
    var n = this.numY;
    var u =
      (this.u[i * n + j - 1] +
        this.u[i * n + j] +
        this.u[(i + 1) * n + j - 1] +
        this.u[(i + 1) * n + j]) *
      0.25;
    return u;
  }

  avgV(i, j) {
    var n = this.numY;
    var v =
      (this.v[(i - 1) * n + j] +
        this.v[i * n + j] +
        this.v[(i - 1) * n + j + 1] +
        this.v[i * n + j + 1]) *
      0.25;
    return v;
  }

  advectVel(dt) {
    this.newU.set(this.u);
    this.newV.set(this.v);

    var n = this.numY;
    var h = this.h;
    var h2 = 0.5 * h;

    for (var i = 1; i < this.numX; i++) {
      for (var j = 1; j < this.numY; j++) {
        //cnt++;

        // u component
        if (
          this.s[i * n + j] != 0.0 &&
          this.s[(i - 1) * n + j] != 0.0 &&
          j < this.numY - 1
        ) {
          var x = i * h;
          var y = j * h + h2;
          var u = this.u[i * n + j];
          var v = this.avgV(i, j);
          //						var v = this.sampleField(x,y, V_FIELD);
          x = x - dt * u;
          y = y - dt * v;
          u = this.sampleField(x, y, this.U_FIELD);
          this.newU[i * n + j] = u;
        }
        // v component
        if (
          this.s[i * n + j] != 0.0 &&
          this.s[i * n + j - 1] != 0.0 &&
          i < this.numX - 1
        ) {
          var x = i * h + h2;
          var y = j * h;
          var u = this.avgU(i, j);
          //						var u = this.sampleField(x,y, U_FIELD);
          var v = this.v[i * n + j];
          x = x - dt * u;
          y = y - dt * v;
          v = this.sampleField(x, y, this.V_FIELD);
          this.newV[i * n + j] = v;
        }
      }
    }

    this.u.set(this.newU);
    this.v.set(this.newV);
  }

  advectSmoke(dt) {
    this.newM.set(this.m);

    var n = this.numY;
    var h = this.h;
    var h2 = 0.5 * h;

    for (var i = 1; i < this.numX - 1; i++) {
      for (var j = 1; j < this.numY - 1; j++) {
        if (this.s[i * n + j] != 0.0) {
          var u = (this.u[i * n + j] + this.u[(i + 1) * n + j]) * 0.5;
          var v = (this.v[i * n + j] + this.v[i * n + j + 1]) * 0.5;
          var x = i * h + h2 - dt * u;
          var y = j * h + h2 - dt * v;

          this.newM[i * n + j] = this.sampleField(x, y, this.S_FIELD);
        }
      }
    }
    this.m.set(this.newM);
  }

  // ----------------- end of simulator ------------------------------

  simulate(dt, gravity, numIters) {
    var t0 = performance.now();

    if (gravity != 0) this.integrate(dt, gravity);

    var t1 = performance.now();

    this.p.fill(0.0);
    this.solveIncompressibility(numIters, dt);

    var t2 = performance.now();

    this.extrapolate();

    var t3 = performance.now();

    this.advectVel(dt);

    var t4 = performance.now();

    this.advectSmoke(dt);

    var t5 = performance.now();

    // calc timers
    this.timers.integration.update(t0, t1);
    this.timers.solveIncompressibility.update(t1, t2);
    this.timers.extrapolate.update(t2, t3);
    this.timers.advectVel.update(t3, t4);
    this.timers.advectSmoke.update(t4, t5);
    this.timers.total.update(t0, t5);
  }
}
