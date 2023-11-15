import TimerStats from "./TimerStats.mjs";
import BoundingBox from "./BoundingBox.mjs";
import Vector from "./Vector.mjs";


const MAX_LEVELS = 4;

function sqr(x) { return x * x; }

const SQRT2 = Math.sqrt(2);

// maintain only one level difference between adjacent cells

class FluidCell {
  constructor(fluid, parentCell, x,y,h,level) {
    this.fluid = fluid;
    this.parentCell = parentCell;
    this.x = x;
    this.y = y;
    this.cx = x + h/2;
    this.cy = y + h/2;
    this.h = h;
    this.area = Math.sqrt(h*h);
    this.u = 0;
    this.v = 0;
    this.inV = 0;
    this.inU = 0;
    this.newU = 0;
    this.newV = 0;
    this.sink = false; // outer edge, can absorb any influx
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

    this.left = [this];
    this.right = [this];
    this.above = [this];
    this.below = [this];

    this.cells = null;
    this.level = level;
    this.childLevel = 0;
    this.updateChildLevel(level);
  }

  reset() {
    this.cells = null;
  }

  setU(u) {
    //if (this.s == 0) this.u = 0;
    if (isNaN(u)) return;
    this.u = u;
  }

  setV(v) {
    if (isNaN(v)) return;
    //if (this.s == 0) this.v = 0;
    this.v = v;
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
      var w = cX(this.bb.bottomRight.x) - cX(this.bb.bottomLeft.x);
      c.fillRect(cX(this.bb.bottomLeft.x), cY(this.bb.topLeft.y), w, w);

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
        c.lineWidth = this.fluid.hoverCell == this ? 2 : 1;
        c.strokeStyle = this.fluid.hoverCell == this ? '#f00' : 'rgba(80,80,80,0.5)';
        c.beginPath();
        c.moveTo(cX(this.bb.bottomLeft.x), cY(this.bb.bottomLeft.y));
        c.lineTo(cX(this.bb.topLeft.x), cY(this.bb.topLeft.y));
        c.lineTo(cX(this.bb.topRight.x), cY(this.bb.topRight.y));
        if (this.fluid.hoverCell == this) {
          c.lineTo(cX(this.bb.bottomRight.x), cY(this.bb.bottomRight.y));
          c.closePath();
        }
        c.stroke();
        c.lineWidth = 1;
    }
  }

  drawVelocities(c, cX, cY, scale) {
    if (this.cells) {
      this.cells.forEach((cell)=>{
        cell.drawVelocities(c, cX, cY, scale);
      })
    } else {
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
      //c.moveTo(x0, y);
      //c.lineTo(cX(this.uvx), cY(this.uvy));

      c.stroke();

      var x = cX(this.vp.x);
      var y0 = cY(this.vp.y);
      var y1 = cY(this.vp.y + v * scale);

      c.strokeStyle = "#f00";
      c.beginPath();
      c.moveTo(x, y0);
      c.lineTo(x, y1);

      // now show advect Point
      //c.moveTo(x, y0);
      //c.lineTo(cX(this.vvx), cY(this.vvy));

      c.stroke();
    }
  }

  drawDiagnostics(c, cX, cY) {
    if (!this.fluid.hoverCell == this) return;

    // draw lines to neighbours
    c.strokeStyle = "rgba(0,0,0,0.3)";
    c.lineWidth = 2;

    c.beginPath();
    
    var cx = cX(this.x + this.h/2);
    var cy = cY(this.y + this.h/2);

    this.left.forEach((a)=>{
      c.moveTo(cx, cy);
      c.lineTo(cX(a.x + a.h/2), cY(a.y + a.h/2));
    });

    this.above.forEach((a)=>{
      c.moveTo(cx, cy);
      c.lineTo(cX(a.x + a.h/2), cY(a.y + a.h/2));
    });

    this.right.forEach((a)=>{
      c.moveTo(cx, cy);
      c.lineTo(cX(a.x + a.h/2), cY(a.y + a.h/2));
    });

    this.below.forEach((a)=>{
      c.moveTo(cx, cy);
      c.lineTo(cX(a.x + a.h/2), cY(a.y + a.h/2));
    });

    c.stroke();
    c.lineWidth = 1;
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
        var x1 = this.x + (i + 0) * this.h/2;
        var y1 = this.y + (j + 0) * this.h/2;
        var cell = new FluidCell(this.fluid, this, x1,y1,this.h/2, this.level+1);
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
    var x0 = x < this.cx ? 0 : 1;
    var y0 = y < this.cy ? 0 : 1;

    var cell = this.cells[x0 * 2 + y0];

    if (cell) {
      return cell.getCellAt(x,y);  
    } 
    console.error(this.x, this.y, x,y, x0, y0);
    return this;
  }


  boundaryCells(i1,i2) {
    if (!this.cells) return [this];

    var res = [];
    res = this.cells[i1].boundaryCells(i1,i2);
    res = res.concat(this.cells[i2].boundaryCells(i1,i2));
    return res;
  }


  coordinateOverlap(p1, w1, p2, w2) {
    var v = Math.min(p1+w1, p2+w2) - Math.max(p1, p2);
    return v > 0;
  }

  updateNeighbours() {
    // expand bordering cells (of any level)
    if (this.above.length == 1) this.above = this.above[0].boundaryCells(0*2+0, 1*2+0);
    if (this.below.length == 1) this.below = this.below[0].boundaryCells(0*2+1, 1*2+1);
    if (this.left.length == 1) this.left = this.left[0].boundaryCells(1*2+0, 1*2+1);
    if (this.right.length == 1) this.right = this.right[0].boundaryCells(0*2+0, 0*2+1);

    // prune boundary cells for only those that overlap with this
    _.remove(this.above, (cell, index) => {
      return !this.coordinateOverlap(this.x, this.h, cell.x, cell.h);
    });
    _.remove(this.below, (cell, index) => {
      return !this.coordinateOverlap(this.x, this.h, cell.x, cell.h);
    });
    _.remove(this.left, (cell, index) => {
      return !this.coordinateOverlap(this.y, this.h, cell.y, cell.h);
    });
    _.remove(this.right, (cell, index) => {
      return !this.coordinateOverlap(this.y, this.h, cell.y, cell.h);
    });


    if (!this.cells) return;


    // bottom left
    this.cells[0 * 2 + 0].left = _.clone(this.left);
    this.cells[0 * 2 + 0].above = [this.cells[0 * 2 + 1]];
    this.cells[0 * 2 + 0].right = [this.cells[1 * 2 + 0]];
    this.cells[0 * 2 + 0].below = _.clone(this.below);

    // top left
    this.cells[0 * 2 + 1].left = _.clone(this.left);
    this.cells[0 * 2 + 1].above = _.clone(this.above);
    this.cells[0 * 2 + 1].right = [this.cells[1 * 2 + 1]];
    this.cells[0 * 2 + 1].below = [this.cells[0 * 2 + 0]];

    // bottom right
    this.cells[1 * 2 + 0].left = [this.cells[0 * 2 + 0]];
    this.cells[1 * 2 + 0].above = [this.cells[1 * 2 + 1]];
    this.cells[1 * 2 + 0].right = _.clone(this.right);
    this.cells[1 * 2 + 0].below = _.clone(this.below);

    // top right
    this.cells[1 * 2 + 1].left = [this.cells[0 * 2 + 1]]
    this.cells[1 * 2 + 1].above = _.clone(this.above);
    this.cells[1 * 2 + 1].right = _.clone(this.right);
    this.cells[1 * 2 + 1].below = [this.cells[1 * 2 + 0]];
    
    
    // update neighbours of child cells
    for (var i = 0; i < 2; i++) {
      for (var j = 0; j < 2; j++) {
        var index = i * 2 + j;
        this.cells[index].updateNeighbours();
      }
    }
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


  resetPressure() {
    this.p = 0;
    if (this.cells) {
      this.cells.forEach((cell)=>{
        cell.resetPressure();
      });
    }
  }


  solveIncompressibility(relax, cp, dt) {

    if (this.cells) {
      this.cells.forEach((cell)=>{
        cell.solveIncompressibility(relax, cp, dt);
      });

      // aggregate values from children to this level
      //this.u = this.cells[0 * 2 + 0].u + this.cells[0 * 2 + 1].u;
      //this.v = this.cells[0 * 2 + 0].v + this.cells[1 * 2 + 0].v;

      return;
    }

    // this cell solid if s=0, so skip as nothing to solve
    if (this.s == 0.0) {
      // force u and v to zero
      this.setU(0);
      this.setV(0);
      return;
    }

    var offset = this.h/16;

    // get neighbours
    var left = this.left;
    var right = this.right;
    var above = this.above;
    var below = this.below;
    
    this.sx0 = _.reduce(left, (sum, a)=>{ return  sum + a.s * Math.min(a.h / this.h,1) }, 0);
    this.sx1 = _.reduce(right, (sum, a)=>{ return  sum + a.s * Math.min(a.h / this.h,1) }, 0);
    this.sy0 = _.reduce(below, (sum, a)=>{ return  sum + a.s * Math.min(a.h / this.h,1) }, 0);
    this.sy1 = _.reduce(above, (sum, a)=>{ return  sum + a.s * Math.min(a.h / this.h,1) }, 0);
    var s = this.sx0 + this.sx1 + this.sy0 + this.sy1;

    // if all neighbour cells are solid (s=0), then also nothing to solve
    if (s == 0.0) {
      // force u and v to zero
      this.setU(0);
      this.setV(0);
      return;
    }
    
    // total divergence (outflow)
    var div =
      _.reduce(right, (sum, a)=>{ return  sum + a.u * Math.min(this.h / a.h, 1) }, 0) - // right boundary
      this.u + // left boundary
      _.reduce(above, (sum, a)=>{ return  sum + a.v * Math.min(this.h / a.h, 1) }, 0) - // top boundary
      this.v; // bottom boundary


    this.div = div;
    if (this.sink) div = 0;

    var p = -div / s;
    p *= relax;
    this.p += cp * p * this.h;

    // adjust velocity vectors to rebalance in vs out flow
    this.u -= this.sx0 * p;
    right.forEach((a)=>{
      if (a.s == 1) a.u += this.sx1 * p / Math.max(this.h / a.h, 1);
    })
    //right.u += sx1 * p;
    this.v -= this.sy0 * p;
    //above.v += sy1 * p;
    above.forEach((a)=>{
      if (a.s == 1) a.v += this.sy1 * p / Math.max(this.h / a.h, 1);
    })

    // quality check
    var q = div;
    this.incompressibility = q;
  }


  startAdvectVel() {
    this.newU = this.u;
    this.newV = this.v;

    // start children
    if (this.cells) {
      this.cells.forEach((cell)=>{
        cell.startAdvectVel();
      });
    }
  }

  advectVel(dt) {

    if (this.cells) {
      this.cells.forEach((cell)=>{
        cell.advectVel(dt);
      });

      return;
    }

    // u component...  
    if (this.s != 0.0 && this.sx0 != 0 && !this.sink) {
      var x = this.up.x;
      var y = this.up.y;
      var u = this.u; // starting u
      //var v = this.avgV(i, j);
      var v = this.fluid.sampleVField(x,y);
      var x1 = x - dt * u;
      var y1 = y - dt * v;

      if (isNaN(v)) {
        console.error('blergh v',x,y,x1,y1);
        this.fluid.errorPoints.push(new Vector(x,y));
      }
      
      try {
        if (!isNaN(x1)) this.uvx = x1;
        if (!isNaN(y1)) this.uvy = y1;
        u = this.fluid.sampleUField(x1, y1);
      } catch(e) {
        console.error(this,x,y,x1,y1,u,v);
      }
      if (isNaN(u)) {
        console.error(this,x,y,x1,y1,u,v);
      } else if (!isNaN(u)) this.newU = u;
    }

    // v component... 
    //if (cell.s != 0.0 && cell.below.s != 0 && i < this.numX - 1) {
    if (this.s != 0.0 && this.sy0 != 0 && !this.sink) {
      var x = this.vp.x;
      var y = this.vp.y;

      //var u = this.avgU(i, j);
      var u = this.fluid.sampleUField(x,y);
      
      if (isNaN(u)) {
        console.error('blergh u', x,y);
        this.fluid.errorPoints.push(new Vector(x,y));
      }

      var v = this.v;
      var x1 = x - dt * u;
      var y1 = y - dt * v;

      try {
        if (!isNaN(x1)) this.vvx = x1;
        if (!isNaN(x1)) this.vvy = y1;
        v = this.fluid.sampleVField(x1, y1);
      } catch(e) {
        console.error(this,x,y,u,v);
      }
      
      if (!isNaN(v)) this.newV = v;
    }
  }

  endAdvectVel() {
    this.setU(this.newU);
    this.setV(this.newV);

    // start children
    if (this.cells) {
      this.cells.forEach((cell)=>{
        cell.endAdvectVel();
      });
    }
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

    this.U_FIELD = 'u';
    this.V_FIELD = 'v';
    this.S_FIELD = 'm';

    this.scene = null;

    this.updateNeighbours();

    this.errorPoints = [];

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

  reset() {
    // reset the quadtree grid
    for (var i = 0; i < this.numX; i++) {
      for (var j = 0; j < this.numY; j++) {
        var ci = i * this.numY + j;
        var cell = this.cells[ci];

        cell.reset();
      }
    }

    this.updateNeighbours();
  }

  updateNeighbours() {
    // update left, right, above, below references
    // blanket set all top level references
    for (var i = 0; i < this.numX; i++) {
      for (var j = 0; j < this.numY; j++) {
        var ci = i * this.numY + j;
        var cell = this.cells[ci];

        if (i>0) cell.left = [this.cells[(i-1) * this.numY + j]];
        if (i<this.numX-1) cell.right = [this.cells[(i+1) * this.numY + j]];

        if (j>0) cell.below = [this.cells[(i) * this.numY + j - 1]];
        if (j<this.numY-1) cell.above = [this.cells[(i) * this.numY + j + 1]];
      }
    }

    // then recurse:
    for (var i = 0; i < this.numX; i++) {
      for (var j = 0; j < this.numY; j++) {
        var ci = i * this.numY + j;
        var cell = this.cells[ci];
        cell.updateNeighbours();
      }
    }
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

    this.updateNeighbours();
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
          cell.subdivide(maxLevel-1);
        }

      }
    }

    this.updateNeighbours();
  }

  solveIncompressibility(numIters, dt) {
    var n = this.numY;
    //var cp = (this.density * this.h) / dt;
    var cp = this.density / dt;

    for (var iter = 0; iter < numIters; iter++) {
      // tween relaxation toward 1
      var relax = iter == numIters - 1 ? 1 : this.scene.overRelaxation;

      for (var i = 1; i < this.numX; i++) {
        for (var j = 1; j < this.numY - 1; j++) {

          this.cells[i * n + j].solveIncompressibility(relax, cp, dt);
          
        }
      }
    }
  }

  extrapolate() {
    var n = this.numY;
    // maintain flow around the boundary
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

    var cell;
    try {
      cell = this.getCellAt(x,y);
    } catch(e) {
      console.error(x,y,cell);
    }

    var tx = (x - cell.x) / cell.h;
    var ty = (y - cell.y) / cell.h;

    // get neighbours
    var above = cell.above[0];
    var below = cell.below[0];
    var right = cell.right[0];

    var sx = 1.0 - tx;
    var sy = 1.0 - ty;

    var val = 0;
    
    if (!cell || !above || !below || !right ) {
      console.error('smeg');
      return 0;
    }

    val =
        sx * sy * cell[field] +
        tx * sy * right[field] +
        tx * ty * below[field] +
        sx * ty * above[field];

    if (isNaN(val)) {
      console.error(tx,ty,sx,sy,field, cell[field], above[field], below[field], right[field]);
    }

    return val;
  }

  sampleUField(x, y) {
    var n = this.numY;
    var h = this.h;
    var h1 = 1.0 / h;
    var h2 = 0.5 * h;

    // sample a linear weighted average of values from around the point (x,y)
    // if y is above midpoint of cell then use the above, right and above right cells as reference
    // otherwise use below, right and below right as reference

    // 
    var bl, br, tl, tr;

    bl = this.getCellAt(x, y);
    // if y above midpoint
    if (y >= bl.up.y) {
      tl = this.getCellAt(x, bl.y + bl.h);
      br = this.getCellAt(bl.x + bl.h, y);
      tr = this.getCellAt(bl.x + bl.h, bl.y + bl.h);

    } else {
      tl = bl;
      tr = this.getCellAt(tl.x + tl.h, y);
      bl = this.getCellAt(x, tl.y - Number.EPSILON);
      br = this.getCellAt(tl.x + tl.h, tl.y - Number.EPSILON);
    }

    var tlv = tl.u;
    var trv = tr.u;
    var blv = bl.u;
    var brv = br.u;
  
    var val = 0;
  
    // lerp factors ... combined lerps
    var brbl = (br.up.x - bl.up.x);
    var tx1 = brbl > 0 ? (x - bl.up.x) / brbl : 0;
    var trtl = (tr.up.x - tl.up.x);
    var tx2 = trtl > 0 ? (x - tl.up.x) / trtl : 0;
    var tx = (tx1 + tx2) / 2;

    var tlbl = (tl.up.y - bl.up.y);
    var ty1 = tlbl > 0 ? (y - bl.up.y) / tlbl : 0;
    var trbr = (tr.up.y - br.up.y);
    var ty2 = trbr > 0 ? (y - br.up.y) / trbr : 0;
    var ty = (ty1 + ty2) / 2;
  
    var sx = 1.0 - tx;
    var sy = 1.0 - ty;

    val =
        sx * sy * blv +
        tx * sy * brv +
        tx * ty * trv +
        sx * ty * tlv;

    if (isNaN(val)) {
      console.error(tx,ty,sx,sy, bl, br, tl, tr);
    }

    return val;
  }

  sampleVField(x, y) {
    var n = this.numY;
    var h = this.h;
    var h1 = 1.0 / h;
    var h2 = 0.5 * h;

    // sample a linear weighted average of values from around the point (x,y)
    // if x is left of midpoint of cell then use the left, above and above left cells as reference
    // otherwise use right, above and above right as reference

    // 
    var bl, br, tl, tr;

    bl = this.getCellAt(x, y);
    // if x to right of midpoint
    if (x >= bl.vp.x) {
      tl = this.getCellAt(x, bl.y + bl.h);
      br = this.getCellAt(bl.x + bl.h, y);
      tr = this.getCellAt(bl.x + bl.h, bl.y + bl.h);

    } else {
      br = bl;
      bl = this.getCellAt(br.x - Number.EPSILON, y);
      tl = this.getCellAt(br.x - Number.EPSILON, br.y + br.h);
      tr = this.getCellAt(x, br.y + br.h);
    }

    var tlv = tl.v;
    var trv = tr.v;
    var blv = bl.v;
    var brv = br.v;
  
    var val = 0;
  
    // lerp factors ... combined lerps
    var brbl = (br.vp.x - bl.vp.x);
    var tx1 = brbl > 0 ? (x - bl.vp.x) / brbl : 0;
    var trtl = (tr.vp.x - tl.vp.x);
    var tx2 = trtl > 0 ? (x - tl.vp.x) / trtl : 0;
    var tx = (tx1 + tx2) / 2;

    var tlbl = (tl.vp.y - bl.vp.y);
    var ty1 = tlbl > 0 ? (y - bl.vp.y) / tlbl : 0;
    var trbr = (tr.vp.y - br.vp.y);
    var ty2 = trbr > 0 ? (y - br.vp.y) / trbr : 0;
    var ty = (ty1 + ty2) / 2;
  
    var sx = 1.0 - tx;
    var sy = 1.0 - ty;

    val =
        sx * sy * blv +
        tx * sy * brv +
        tx * ty * trv +
        sx * ty * tlv;
    

    if (isNaN(val)) {
      console.error(tx,ty,sx,sy, bl, br, tl, tr);
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
    for (var i = 0; i < this.numX; i++) {
      for (var j = 0; j < this.numY; j++) {
        this.cells[i*this.numY + j].startAdvectVel();
      }
    }

    for (var i = 1; i < this.numX; i++) {
      for (var j = 1; j < this.numY-1; j++) {
        //cnt++;
        var ci = i * this.numY + j;
        var cell = this.cells[ci];
        
        cell.advectVel(dt);
      }
    }

    for (var i = 0; i < this.numX; i++) {
      for (var j = 0; j < this.numY; j++) {
        this.cells[i*this.numY + j].endAdvectVel();
      }
    }
  }

  advectSmoke(dt) {
    for (var i = 0; i < this.numX; i++) {
      for (var j = 0; j < this.numY; j++) {
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
    for (var i = 0; i < this.numX; i++) {
      for (var j = 0; j < this.numY; j++) {
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
        this.cells[i*this.numY + j].resetPressure();
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
