
import TimerStats from './TimerStats.mjs';

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
            if (this.s[i * n + j] == 0.0) continue;
  
            var s = this.s[i * n + j];
            var sx0 = this.s[(i - 1) * n + j];
            var sx1 = this.s[(i + 1) * n + j];
            var sy0 = this.s[i * n + j - 1];
            var sy1 = this.s[i * n + j + 1];
            var s = sx0 + sx1 + sy0 + sy1;
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
  
      this.integrate(dt, gravity);
  
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