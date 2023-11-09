import TimerStats from "./TimerStats.mjs";
import BoundingBox from "./BoundingBox.mjs";
import Vector from "./Vector.mjs";


const MAX_LEVELS = 4;

// maintain only one level difference between adjacent cells

class FluidCell {
  constructor(fluid, parentCell, x,y,h,level) {
    this.fluid = fluid;
    this.parentCell = parentCell;
    this.x = x;
    this.y = y;
    this.h = h;
    this.u = 0;
    this.v = 0;
    this.inV = 0;
    this.inU = 0;
    this.newU = 0;
    this.newV = 0;
    // pressure
    this.p = 0;
    // solidity
    this.s = 1;  // 1=empty
    this.incompressibility = 0;
    // smoke
    this.m = 1;
    this.newM = 0;

    this.bb = new BoundingBox();
    this.bb.addPoint(new Vector(x, y));
    this.bb.addPoint(new Vector(x+h, y+h));
    this.bb.finalise();

    this.up = new Vector(x, y + h/2);
    this.vp = new Vector(x+h/2, y);
    this.uvx = this.up.x;
    this.uvy = this.up.y;
    this.vvx = this.vp.x;
    this.vvy = this.vp.y;

    this.cells = null;
    this.level = level;
    this.childLevel = 0;
    this.updateChildLevel(level);
  }

  setU(u) {
    if (this.s == 0) this.u = 0;
    else if (this.inU) this.u = this.inU;
    else this.u = u;
  }

  setV(v) {
    if (this.s == 0) this.v = 0;
    else if (this.inV) this.v = this.inV;
    else this.v = v;
  }

  drawPressureFill(c, cX, cY, getSciColor, minP, maxP) {
    if (this.cells) {
      this.cells.forEach((cell)=>{
        cell.drawPressureFill(c, cX, cY, getSciColor, minP, maxP);
      })
    } else {
      c.fillStyle = '#555';
      if (this.s != 0) {
        var color = getSciColor(this.p, minP, maxP);
        
        c.fillStyle = 'rgb('+color[0]+','+color[1]+','+color[2]+')';
      }
      c.fillRect(cX(this.bb.bottomLeft.x), cY(this.bb.topLeft.y), cX(this.h), cX(this.h));

      /*
      c.fillStyle = '#000';
      c.font = "8px Arial";
      c.fillText(this.childLevel.toFixed(0),  cX(this.bb.bottomLeft.x)+1, cY(this.bb.bottomLeft.y)-1);
      */
    }
    /*
    if (this.level == 1) {
      c.fillStyle = '#000';
    c.font = "8px Arial";
    c.fillText(this.childLevel.toFixed(0),  cX(this.bb.bottomLeft.x)+1, cY(this.bb.bottomLeft.y)-1);
    }
    */
  }

  drawBoundary(c, cX, cY) {
    if (this.cells) {
      this.cells.forEach((cell)=>{
        cell.drawBoundary(c, cX, cY);
      })
    } else {
        // left boundary
        c.lineWidth = 1;
        c.strokeStyle = this.hover ? '#f00' : '#888';
        c.beginPath();
        c.moveTo(cX(this.bb.bottomLeft.x), cY(this.bb.bottomLeft.y));
        c.lineTo(cX(this.bb.topLeft.x), cY(this.bb.topLeft.y));
        c.stroke();


        // bottom boundary
        c.beginPath();
        c.moveTo(cX(this.bb.bottomLeft.x), cY(this.bb.bottomLeft.y));
        c.lineTo(cX(this.bb.bottomRight.x), cY(this.bb.bottomRight.y));
        c.stroke();
    }
  }

  drawVelocities(c, cX, cY, scale) {
    var u = this.u;
    var v = this.v;
    
    c.strokeStyle = "#00f";
    c.beginPath();
    var x0 = cX(this.up.x);
    var x1 = cX(this.up.x + u * scale);
    var y = cY(this.up.y);

    c.moveTo(x0, y);
    c.lineTo(x1, y);

    // now show advect Point
    c.moveTo(x0, y);
    c.lineTo(cX(this.uvx), cY(this.uvy));

    c.stroke();

    var x = cX(this.vp.x);
    var y0 = cY(this.vp.y);
    var y1 = cY(this.vp.y + v * scale);

    c.strokeStyle = "#f00";
    c.beginPath();
    c.moveTo(x, y0);
    c.lineTo(x, y1);

    // now show advect Point
    c.moveTo(x, y0);
    c.lineTo(cX(this.vvx), cY(this.vvy));

    c.stroke();
  }

  updateChildLevel(childLevel) {
    if (childLevel > this.childLevel) this.childLevel = childLevel;
    if (this.parentCell) {
      this.parentCell.updateChildLevel(this.childLevel);
    }
  }

  subdivide(toLevel) {
    if (this.cells) return;
    this.cells = [];
    for (var i = 0; i < 2; i++) {
      for (var j = 0; j < 2; j++) {
        var index = i * this.numY + j;
        var x1 = this.x + (i + 0) * this.h/2;
        var y1 = this.y + (j + 0) * this.h/2;
        var cell = new FluidCell(this, this, x1,y1,this.h/2, this.level+1);
        this.cells.push( cell );

        // sud-divide further?
        if (toLevel > this.level+1) {
          cell.subdivide(toLevel);
        }
      }
    }
  }

  getCellAt(x,y) {
    if (!this.cells) return this;

    var h1 = this.h/2;
    var x0 = Math.max(Math.min(Math.floor((x - this.x) / h1), 1),0);
    var y0 = Math.max(Math.min(Math.floor((y - this.y) / h1), 1), 0);

    var cell = this.cells[x0 * 2 + y0];

    if (cell) {
      return cell.getCellAt(x,y);  
    } 
    console.error(this.x, this.y, x,y, x0, y0);
    return this;
  }

  updateObstacle(obstacle) {
    // default all cells to empty
    this.s = 1;

    var cbb = this.bb;
    // compare to bounds
    if (obstacle.bb.overlaps(cbb)) {
      // check the obstacle for containment
      cbb.tlc = obstacle.contains(cbb.topLeft);
      cbb.trc = obstacle.contains(cbb.topRight);
      cbb.blc = obstacle.contains(cbb.bottomLeft);
      cbb.brc = obstacle.contains(cbb.bottomRight);

      // solid (s=0) if entirely contained
      this.s = cbb.tlc && cbb.trc && cbb.blc && cbb.brc ? 0 : 1;
     
      // check for partial containment
      if (this.s == 1 && (cbb.tlc || cbb.trc || cbb.blc || cbb.brc)) {
        // sub-divide if not too deep
        if (this.level <MAX_LEVELS) {
          this.subdivide();

          this.cells.forEach((cell)=>{
            cell.updateObstacle(obstacle);
          });
        }
      }
    }
  }


  solveIncompressibility(relax, cp, dt) {

    if (this.cells) {

    }

    // this cell solid if s=0, so skip as nothing to solve
    if (this.s == 0.0) {
      // force u and v to zero
      this.setU(0);
      this.setV(0);
      return;
    }

    var offset = this.h/16;

    // get neighbours - sample at approx centre of cell, just over the boundary
    var left = this.fluid.getCellAt(this.x - offset, this.up.y);
    var right = this.fluid.getCellAt(this.x + this.h + offset, this.up.y);
    var above = this.fluid.getCellAt(this.vp.x, this.y + this.h + offset);
    var below = this.fluid.getCellAt(this.vp.x, this.y - offset );
    
    var sx0 = left.s; // cell to left
    var sx1 = right.s;  // cell to right
    var sy0 = below.s;  // cell below
    var sy1 = above.s;  // cell above
    var s = sx0 + sx1 + sy0 + sy1;

    // if all neighbour cells are solid (s=0), then also nothing to solve
    if (s == 0.0) {
      // force u and v to zero
      this.setU(0);
      this.setV(0);
      return;
    }
    

    // total divergence (outflow)
    var div =
      right.u - // right boundary
      this.u + // left boundary
      above.v - // top boundary
      this.v; // bottom boundary

    var p = -div / s;
    p *= relax;
    this.p += cp * p * this.h;

    // adjust velocity vectors to rebalance in vs out flow
    this.u -= sx0 * p;
    right.u += sx1 * p;
    this.v -= sy0 * p;
    above.v += sy1 * p;

    // quality check
    var q = div;
    this.incompressibility = q;
  }


}


/*  ======================================================================
*/
export default class Fluid {
  constructor(density, numX, numY, h) {
    this.density = density;
    this.numX = numX + 2;
    this.numY = numY + 2;
    this.numCells = this.numX * this.numY;
    this.h = h;
    var num = numX * numY;
    this.cells = new Array();

    // compute bounds for each cell
    for (var i = 0; i < this.numX; i++) {
      for (var j = 0; j < this.numY; j++) {
        var index = i * this.numY + j;
        var x1 = (i + 0) * this.h;
        var y1 = (j + 0) * this.h;
        this.cells.push(new FluidCell(this, null, x1,y1,h,1) );
      }
    }

    console.log("cellSize", h);

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
        var cell = this.cells[ci];

        cell.updateObstacle(obstacle);
      }
    }
  }

  updateSubdivision() {
    for (var i = 1; i < this.numX - 2; i++) {
      for (var j = 1; j < this.numY - 2; j++) {
        var ci = i * this.numY + j;
        var cell = this.cells[ci];

        // compare neighbour cells, make sure we're subdivided enough
        var maxLevel = Math.max(
          this.cells[(i-1)*this.numY + (j)].childLevel,
          this.cells[(i+1)*this.numY + (j)].childLevel,
          this.cells[(i)*this.numY + (j-1)].childLevel,
          this.cells[(i)*this.numY + (j+1)].childLevel
        );
        if (maxLevel > cell.childLevel+1 && cell.s == 1) {
          //cell.subdivide(maxLevel-1);
        }

      }
    }
  }

  solveIncompressibility(numIters, dt) {
    var n = this.numY;
    //var cp = (this.density * this.h) / dt;
    var cp = this.density / dt;

    for (var iter = 0; iter < numIters; iter++) {
      // tween relaxation toward 1
      var relax = iter == numIters - 1 ? 1 : this.scene.overRelaxation;

      for (var i = 1; i < this.numX - 1; i++) {
        for (var j = 1; j < this.numY - 1; j++) {

          this.cells[i * n + j].solveIncompressibility(relax, cp, dt);
          
        }
      }
    }
  }

  extrapolate() {
    var n = this.numY;
    for (var i = 0; i < this.numX; i++) {
      this.cells[i * n + 0].setU(this.cells[i * n + 1].u);
      this.cells[i * n + this.numY - 1].setU(this.cells[i * n + this.numY - 2].u);
    }
    for (var j = 0; j < this.numY; j++) {
      this.cells[0 * n + j].setV(this.cells[1 * n + j].v);
      this.cells[(this.numX - 1) * n + j].setV(this.cells[(this.numX - 2) * n + j].v);
    }
  }

  getCellAt(x,y) {
    // ensure x and y coords are within bounds
    var x0 = Math.min(Math.floor(x / this.h), this.numX - 1);
    var y0 = Math.min(Math.floor(y / this.h), this.numY - 1);
    if (x0 < 0) x0 = 0;
    if (y0 < 0) y0 = 0;

    return this.cells[x0 * this.numY + y0].getCellAt(x,y);
  }

  sampleField(x, y, field) {
    var n = this.numY;
    var h = this.h;
    var h1 = 1.0 / h;
    var h2 = 0.5 * h;

    // ensure x and y coords are within bounds
    x = Math.max(Math.min(x, this.numX * h), 0);
    y = Math.max(Math.min(y, this.numY * h), 0);

    var cell = this.getCellAt(x,y);

    var tx = (x - cell.x) / cell.h;
    var ty = (y - cell.y) / cell.h;

    var offset = cell.h / 16;

    // get neighbours
    var above = this.getCellAt(cell.x, cell.y + cell.h + offset);
    var below = this.getCellAt(cell.x, cell.y - offset);
    var right = this.getCellAt(cell.x + cell.h + offset, cell.y);

    var sx = 1.0 - tx;
    var sy = 1.0 - ty;

    var val = 0;
    
    switch (field) {
      case this.U_FIELD:
        try {
        val =
        sx * sy * cell.u +
        tx * sy * right.u +
        tx * ty * below.u +
        sx * ty * above.u;
      } catch(e) {
        console.error(x,y,x0,y0,x1,y1, this.cells[x0 * n + y0]);
      }
        break;
      case this.V_FIELD:
        val =
        sx * sy * cell.v +
        tx * sy * right.v +
        tx * ty * below.v +
        sx * ty * above.v;
        break;
      case this.S_FIELD:
        try {
          val =
        sx * sy * cell.m +
        tx * sy * right.m +
        tx * ty * below.m +
        sx * ty * above.m; 
        } catch(e) {
          console.error(x,y,x0,y0,x1,y1, this.cells[x0 * n + y0]);
        }
        break;
    }

    return val;
  }

  avgU(i, j) {
    var n = this.numY;
    var u = (this.cells[i*n + j-1].u + this.cells[i*n + j].u +
      this.cells[(i+1)*n + j-1].u + this.cells[(i+1)*n + j].u) * 0.25;
    return u;
      
  }

  avgV(i, j) {
    var n = this.numY;
    var v = (this.cells[(i-1)*n + j].v + this.cells[i*n + j].v +
      this.cells[(i-1)*n + j+1].v + this.cells[i*n + j+1].v) * 0.25;
    return v;
  }

  advectVel(dt) {
    for (var i = 1; i < this.numX - 1; i++) {
      for (var j = 1; j < this.numY - 1; j++) {
        this.cells[i*this.numY + j].newU = this.cells[i*this.numY + j].u;
        this.cells[i*this.numY + j].newV = this.cells[i*this.numY + j].v;
      }
    }

    var n = this.numY;
    var h = this.h;
    var h2 = 0.5 * h;

    for (var i = 1; i < this.numX; i++) {
      for (var j = 1; j < this.numY; j++) {
        //cnt++;
        var ci = i * n + j;
        var cell = this.cells[ci];
        var cbb = cell.bb;

        // u component...  
        if (cell.s != 0.0 && this.cells[(i-1)*n+j].s != 0 && j < this.numY - 1) {
          var x = cell.up.x;
          var y = cell.up.y;
          var u = cell.u; // starting u
          //var v = this.avgV(i, j);
          var v = this.sampleField(x,y, this.V_FIELD);
          x = x - dt * u;
          y = y - dt * v;
          
          try {
            cell.uvx = x;
            cell.uvy = y;
            u = this.sampleField(x, y, this.U_FIELD);
          } catch(e) {
            console.error(i,j,x,y,u,v);
          }
          if (isNaN(u)) {
            console.error(i,j,x,y,u,v);
          } else cell.newU = u;
        }

        // v component... 
        if (cell.s != 0.0 && this.cells[i*n+(j-1)].s != 0 && i < this.numX - 1) {
          var x = cell.vp.x
          var y = cell.vp.y;

          //var u = this.avgU(i, j);
          var u = this.sampleField(x,y, this.U_FIELD);
          var v = cell.v;
          x = x - dt * u;
          y = y - dt * v;

          cell.vvx = x;
          cell.vvy = y;
          
          v = this.sampleField(x, y, this.V_FIELD);
          if (!isNaN(v)) cell.newV = v;
        }
      }
    }

    for (var i = 1; i < this.numX - 1; i++) {
      for (var j = 1; j < this.numY - 1; j++) {
        this.cells[i*this.numY + j].setU(this.cells[i*this.numY + j].newU);
        this.cells[i*this.numY + j].setV(this.cells[i*this.numY + j].newV);
      }
    }
  }

  advectSmoke(dt) {
    for (var i = 1; i < this.numX - 1; i++) {
      for (var j = 1; j < this.numY - 1; j++) {
        this.cells[i*this.numY + j].newM = this.cells[i*this.numY + j].m;
      }
    }

    var n = this.numY;
    var h = this.h;
    var h2 = 0.5 * h;

    for (var i = 1; i < this.numX - 1; i++) {
      for (var j = 1; j < this.numY - 1; j++) {
        if (this.cells[i * n + j].s != 0.0) {
          var u = (this.cells[i * n + j].u + this.cells[(i + 1) * n + j].u) * 0.5;
          var v = (this.cells[i * n + j].v + this.cells[i * n + j + 1].v) * 0.5;
          var x = i * h + h2 - dt * u;
          var y = j * h + h2 - dt * v;

          this.cells[i * n + j].newM = this.sampleField(x, y, this.S_FIELD);
        }
      }
    }
    for (var i = 1; i < this.numX - 1; i++) {
      for (var j = 1; j < this.numY - 1; j++) {
        this.cells[i*this.numY + j].m = this.cells[i*this.numY + j].newM;
      }
    }
  }

  // ----------------- end of simulator ------------------------------

  simulate(dt, gravity, numIters) {
    var t0 = performance.now();

    this.updateSubdivision();

    var t1 = performance.now();

    //this.p.fill(0.0);
    for (var i = 1; i < this.numX - 1; i++) {
      for (var j = 1; j < this.numY - 1; j++) {
        this.cells[i*this.numY + j].p = 0;
      }
    }

    this.solveIncompressibility(numIters, dt);

    var t2 = performance.now();

    this.extrapolate();

    var t3 = performance.now();

    this.advectVel(dt);

    var t4 = performance.now();

    //this.advectSmoke(dt);

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
