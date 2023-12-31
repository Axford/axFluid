import Vector from "./Vector.mjs";
import Line from "./Line.mjs";
import TimerStats from "./TimerStats.mjs";
import Fluid from "./Fluid.mjs";
import Polygon from "./Polygon.mjs";

var canvas = document.getElementById("myCanvas");
var c = canvas.getContext("2d", { willReadFrequently: true });
canvas.width = window.innerWidth - 20;
canvas.height = window.innerHeight - 100;

canvas.focus();

var gui;

var f;

var simHeight = 1.2;
var cScale = canvas.height / simHeight;
var simWidth = canvas.width / cScale;

var hoverX = 0;
var hoverY = 0;

var airfoilSteps = 200;
var airfoilUpperProfile = [];
var airfoilLowerProfile = [];

// obstacle is a polygon
var obstacle = new Polygon();

var frame = 0;
var fps = 0;
var fpsTimer = 0;

function cX(x) {
  return x * cScale;
}

function cY(y) {
  return canvas.height - y * cScale;
}

function setupAirfoil() {
  airfoilUpperProfile = [];
  airfoilLowerProfile = [];

  var at = 0.15; // 18% thickness
  var m = 0.04; // camber
  var p = 0.4; // point of max camber

  for (var i = 0; i < airfoilSteps; i++) {
    var ax = i / (airfoilSteps - 1);

    var yc;
    if (ax <= p) {
      yc = (m / (p * p)) * (2 * p * ax - ax * ax);
    } else {
      yc = (m / Math.pow(1 - p, 2)) * (1 - 2 * p + 2 * p * ax - ax * ax);
    }

    var yt =
      5 *
      at *
      (0.2969 * Math.sqrt(ax) -
        0.126 * ax -
        0.3516 * Math.pow(ax, 2) +
        0.2843 * Math.pow(ax, 3) -
        0.1015 * Math.pow(ax, 4));

    var theta = Math.atan2(yc, ax);

    var xu = ax - yt * Math.sin(theta);
    var yu = yc + yt * Math.cos(theta);

    var xl = ax + yt * Math.sin(theta);
    var yl = yc - yt * Math.cos(theta);

    airfoilUpperProfile.push(new Vector(xu, yu));
    airfoilLowerProfile.push(new Vector(xl, yl));
  }
}

// ----------------- start of simulator ------------------------------

var scene = {
  gravity: -9.81,
  dt: 1.0 / 120.0,
  numIters: 100,
  frameNr: 0,
  overRelaxation: 1.9,
  obstacleX: 0.0,
  obstacleY: 0.0,
  obstacleRadius: 0.15,
  paused: false,
  sceneNr: 0,
  showAirfoilBoundary: false,
  showStreamlines: false,
  showVelocities: true,
  showPressure: false,
  showSmoke: true,
  showPressureDistribution:true,
  showCells:true,
  fluid: null,
  airfoilAngle: 0
};

function setupScene(sceneNr = 0) {
  console.log("setupScene", sceneNr);

  scene.sceneNr = sceneNr;
  scene.obstacleRadius = 0.15;
  scene.overRelaxation = 1.9;

  scene.dt = 1.0 / 60.0;
  //scene.numIters = 40;
  scene.numIters = 100;

  var res = 100;

  if (sceneNr == 0) res = 50;
  else if (sceneNr == 1) res = 20;
  else if (sceneNr == 1) res = 200;
  else if (sceneNr == 3) res = 400;

  var domainHeight = 1.0;
  var domainWidth = (domainHeight / simHeight) * simWidth;
  var h = domainHeight / res; // height of a grid cell

  var numX = Math.floor(domainWidth / h);
  var numY = Math.floor(domainHeight / h);

  console.log('numX', numX);
  console.log('numY', numY);

  var density = 5000.0;

  f = scene.fluid = new Fluid(density, numX, numY, h);
  f.scene = scene;

  var n = f.numY;

  if (sceneNr >=1 && sceneNr <= 3) {
    // vortex shedding
    // wind tunnels

    // wind speed
    var inVel = 1;
    for (var i = 0; i < f.numX; i++) {
      for (var j = 0; j < f.numY; j++) {
        var s = 1.0; // fluid
        if (i == 0 || j == 0 || j == f.numY - 1) s = 0.0; // solid
        f.s[i * n + j] = s;
        f.cellBounds[i * n + j].ua = s;
        f.cellBounds[i * n + j].va = s;

        if (i == 1) {
          f.u[i * n + j] = inVel;
          f.cellBounds[i * n + j].ua = 0;
        }
      }
    }

    // smoke injection
    var pipeH = 0.3 * f.numY;
    var minJ = Math.floor(0.5 * f.numY - 0.5 * pipeH);
    var maxJ = Math.floor(0.5 * f.numY + 0.5 * pipeH);

    //for (var j = minJ; j < maxJ; j++) f.m[j] = 0.0;

    for (var j = 0; j < f.numY; j++) {
      if (j % 20 < 4) f.m[j] = 0;
    }

    setObstacle(0.401, 0.501, true);

    scene.overRelaxation = 1.9;

    scene.gravity = 0;
    scene.showPressure = true;
    scene.showSmoke = false;
    scene.showStreamlines = false;
    scene.showVelocities = true;

    if (sceneNr > 1 ) {
      scene.dt = 1.0 / 120.0;
      scene.numIters = 100;
      scene.showPressure = true;
      scene.showVelocities = false;
    }
  }
}

// draw -------------------------------------------------------

function setColor(r, g, b) {
  c.fillStyle = `rgb(
			${Math.floor(255 * r)},
			${Math.floor(255 * g)},
			${Math.floor(255 * b)})`;
  c.strokeStyle = `rgb(
			${Math.floor(255 * r)},
			${Math.floor(255 * g)},
			${Math.floor(255 * b)})`;
}

function getSciColor(val, minVal, maxVal) {
  val = Math.min(Math.max(val, minVal), maxVal - 0.0001);
  var d = maxVal - minVal;
  val = d == 0.0 ? 0.5 : (val - minVal) / d;
  var m = 0.25;
  var num = Math.floor(val / m);
  var s = (val - num * m) / m;
  var r, g, b;

  switch (num) {
    case 0:
      r = 0.0;
      g = s;
      b = 1.0;
      break;
    case 1:
      r = 0.0;
      g = 1.0;
      b = 1.0 - s;
      break;
    case 2:
      r = s;
      g = 1.0;
      b = 0.0;
      break;
    case 3:
      r = 1.0;
      g = 1.0 - s;
      b = 0.0;
      break;
  }

  return [255 * r, 255 * g, 255 * b, 255];
}

function draw() {
  c.clearRect(0, 0, canvas.width, canvas.height);

  c.fillStyle = "#FF0000";
  f = scene.fluid;
  var n = f.numY;
  var cellScale = 1;

  var h = f.h;

  var minP = f.p[0];
  var maxP = f.p[0];

  var minQ = f.incompressibility[0];
  var maxQ = f.incompressibility[0];

  for (var i = 0; i < f.numCells; i++) {
    minP = Math.min(minP, f.p[i]);
    maxP = Math.max(maxP, f.p[i]);

    minQ = Math.min(minQ, f.incompressibility[i]);
    maxQ = Math.max(maxQ, f.incompressibility[i]);
  }

  //minQ = -2;
  //maxQ = 2;

  var id = c.getImageData(0, 0, canvas.width, canvas.height);

  var color = [255, 255, 255, 255];

  for (var i = 0; i < f.numX; i++) {
    for (var j = 0; j < f.numY; j++) {
      if (f.s[i * n + j] == 0) {
        color[0] = 80;
        color[1] = color[0];
        color[2] = color[0];
      } else if (scene.showPressure) {
        var p = f.p[i * n + j];
        var s = f.m[i * n + j];
        color = getSciColor(p, minP, maxP);
        //color[0] *= f.s[i * n + j];
        //color[1] *= f.s[i * n + j];
        //color[2] *= f.s[i * n + j];
        if (scene.showSmoke) {
          color[0] = Math.max(0.0, color[0] - 255 * s);
          color[1] = Math.max(0.0, color[1] - 255 * s);
          color[2] = Math.max(0.0, color[2] - 255 * s);
        }
      } else {
        var s = f.m[i * n + j];
        var q = f.incompressibility[i * n + j];
        color = getSciColor(q, minQ, maxQ);
        //color[0] *= f.s[i * n + j];
        //color[1] *= f.s[i * n + j];
        //color[2] *= f.s[i * n + j];
        if (scene.showSmoke) {
          color[0] = Math.max(0.0, color[0] - 255 * s);
          color[1] = Math.max(0.0, color[1] - 255 * s);
          color[2] = Math.max(0.0, color[2] - 255 * s);
        }
      }

      var x = Math.floor(cX(i * h));
      var y = Math.floor(cY((j + 1) * h));
      var cx = Math.floor(cScale * cellScale * h) + 1;
      var cy = Math.floor(cScale * cellScale * h) + 1;

      var r = color[0];
      var g = color[1];
      var b = color[2];

      for (var yi = y; yi < y + cy; yi++) {
        var p = 4 * (yi * canvas.width + x);

        for (var xi = 0; xi < cx; xi++) {
          id.data[p++] = r;
          id.data[p++] = g;
          id.data[p++] = b;
          id.data[p++] = 255;
        }
      }
    }
  }

  c.putImageData(id, 0, 0);

  
  // show grid and ua / va boundaries
  if (scene.showCells) {
    c.lineWidth = 1.0;

    for (var i = 0; i < f.numX; i++) {
      for (var j = 0; j < f.numY; j++) {
        var ci = i * n + j;
        var bb = f.cellBounds[ci];

        // left boundary
      c.strokeStyle = (hoverX == i && hoverY == j) ? '#0f0' : (bb.ua < 1 ? "#f00" : '#aaa');
      // draw intersection volume
      c.beginPath();
      c.moveTo(cX(bb.bottomLeft.x), cY(bb.bottomLeft.y));
      c.lineTo(cX(bb.topLeft.x), cY(bb.topLeft.y));
      c.stroke();
      

      // bottom boundary
      c.strokeStyle = (hoverX == i && hoverY == j) ? '#0f0' : (bb.va < 1 ? "#f00" : '#aaa');
      // draw intersection volume
      c.beginPath();
      c.moveTo(cX(bb.bottomLeft.x), cY(bb.bottomLeft.y));
      c.lineTo(cX(bb.bottomRight.x), cY(bb.bottomRight.y));
      c.stroke();
    
        
      }
    }

    c.lineWidth = 1.0;
  }
  

  if (scene.showVelocities) {
    var scale = 0.02;

    for (var i = 0; i < f.numX; i++) {
      for (var j = 0; j < f.numY; j++) {
        var u = f.u[i * n + j];
        var v = f.v[i * n + j];
        //var u = f.cellBounds[i * n + j].ua;
        //var v = f.cellBounds[i * n + j].va;

        var ci = i * n + j;
        var bb = f.cellBounds[ci];

        c.strokeStyle = "#000";
        c.beginPath();
        var x0 = cX(bb.up.x);
        var x1 = cX(bb.up.x + u * scale);
        var y = cY(bb.up.y);

        c.moveTo(x0, y);
        c.lineTo(x1, y);

        // now show advect Point
        if (bb.ua > 0 && bb.uvx != 0) {
          c.moveTo(x0, y);
          c.lineTo(cX(bb.uvx), cY(bb.uvy));
        }
        c.stroke();

        var x = cX(bb.vp.x);
        var y0 = cY(bb.vp.y);
        var y1 = cY(bb.vp.y + v * scale);

        c.strokeStyle = "#f00";
        c.beginPath();
        c.moveTo(x, y0);
        c.lineTo(x, y1);

        // now show advect Point
        if (bb.va > 0 && bb.vvy != 0) {
          c.moveTo(x, y0);
          c.lineTo(cX(bb.vvx), cY(bb.vvy));
        }
        c.stroke();

        
      }
    }
  }

  if (scene.showStreamlines) {
    var segLen = f.h * 0.2;
    var xStepSize = Math.max(Math.floor(f.numX / 15), 1);
    var yStepSize = Math.max(Math.floor(f.numY / 30),1);
    var numSegs = Math.floor(xStepSize * 5);

    console.log(xStepSize, yStepSize, numSegs);

    
    for (var i = 1; i < f.numX - 2; i += xStepSize) {
      for (var j = 1; j < f.numY - 2; j += yStepSize) {
        var x = (i + 0.5) * f.h;
        var y = (j + 0.5) * f.h;

        c.strokeStyle = "#000000";
        c.lineWidth = 1.0;
        c.beginPath();
        c.moveTo(cX(x), cY(y));

        for (var k = 0; k < numSegs; k++) {
          // check if outside bounds
          if (x < 0 || x > f.numX * f.h || y < 0 || y > f.numY * f.h) break;

          var u = f.sampleField(x, y, f.U_FIELD);
          var v = f.sampleField(x, y, f.V_FIELD);
          var l = Math.sqrt(u * u + v * v);
          x += (u / l) * segLen;
          y += (v / l) * segLen;



          // check if we're inside a solid region
          if (
            obstacle.bb.containsCoords(x, y) &&
            obstacle.contains(new Vector(x, y))
          ) {
            c.strokeStyle = "#f00";
            break;
          }

          c.lineTo(cX(x), cY(y));
        }
        c.stroke();

      }
    }
    
  }

  // airfoil profile and pressure distribution
  if (scene.showAirfoilBoundary) {
    c.lineWidth = 1.0;
    c.strokeStyle = "#fff";

    c.beginPath();
    obstacle.edges.forEach((line) => {
      // draw airfoil profile
      c.moveTo(cX(line.p0.x), cY(line.p0.y));
      c.lineTo(cX(line.p1.x), cY(line.p1.y));  
    });
    c.stroke();
    c.lineWidth = 1.0;
  }

  if (scene.showPressureDistribution) {
    c.lineWidth = 1.0;

    obstacle.edges.forEach((line) => {
      
      // calc normal vector
      var v = line.p1.clone();
      v.subtract(line.p0);
      // scale it
      v.normalize();
      v.rotate90();
      v.multiply(1.1 * f.h);

      // sample pressure at p0 offset by v
      var i = Math.floor((line.p0.x+v.x)/f.h);
      var j = Math.floor((line.p0.y+v.y)/f.h);
      var p = f.p[i*f.numY + j];

      // scale pressure line based on maxQ
      var p = 0.3 * (p / maxP);
      var pSign = p > 0;
      p = Math.abs(p);

      v.multiply(p/f.h);
      //if (!pSign) v.negative();
      
      //draw it 
      c.strokeStyle =  pSign ? '#f00' : '#00f';
      c.beginPath();
      c.moveTo(cX(line.p0.x), cY(line.p0.y));
      c.lineTo(cX(line.p0.x + v.x), cY(line.p0.y + v.y));
      c.stroke();
      
    });
    c.lineWidth = 1.0;
  }

  if (scene.showPressure) {
    var s = "pressure: " + minP.toFixed(0) + " ... " + maxP.toFixed(0) + " N/m";
    c.fillStyle = "#000000";
    c.font = "16px Arial";
    c.fillText(s, 10, 35);
  }

  var s = "incompressibility: " + minQ.toFixed(2) + " ... " + maxQ.toFixed(2);
  c.fillStyle = "#000000";
  c.font = "16px Arial";
  c.fillText(s, 300, 35);

  // draw colour scale

  // draw cell bounds
  if (scene.showAirfoilBoundary) {
    c.lineWidth = 1.0;

    for (var i = 0; i < f.numX; i++) {
      for (var j = 0; j < f.numY; j++) {
        var ci = i * n + j;
        var bb = f.cellBounds[ci];
        if (f.s[ci] < 1 && bb.points && bb.points.length > 0) {
          c.strokeStyle = "#fff";
          // draw intersection volume
          c.beginPath();
          c.moveTo(cX(bb.points[0].x), cY(bb.points[0].y));
          for (var k = 1; k < bb.points.length; k++) {
            c.lineTo(cX(bb.points[k].x), cY(bb.points[k].y));
          }
          c.closePath();
          c.stroke();
        }

        // see if we have any intersections to draw
        if (bb.li) {
          c.fillStyle = "#f00";
          c.beginPath();
          c.arc(cX(bb.li.x), cY(bb.li.y), 3, 0, 2 * Math.PI);
          c.fill();
        }

        if (bb.ti) {
          c.fillStyle = "#f00";
          c.beginPath();
          c.arc(cX(bb.ti.x), cY(bb.ti.y), 3, 0, 2 * Math.PI);
          c.fill();
        }

        if (bb.ri) {
          c.fillStyle = "#00f";
          c.beginPath();
          c.arc(cX(bb.ri.x), cY(bb.ri.y), 2, 0, 2 * Math.PI);
          c.fill();
        }

        if (bb.bi) {
          c.fillStyle = "#00f";
          c.beginPath();
          c.arc(cX(bb.bi.x), cY(bb.bi.y), 2, 0, 2 * Math.PI);
          c.fill();
        }
      }
    }

    c.lineWidth = 1.0;
  }

  var y = 130;
  for (const [key, t] of Object.entries(scene.fluid.timers)) {
    t.draw(c, 20, y);
    y += 20;
  }

  // draw fps
  c.fillStyle = "#00f";
  c.font = "16px Arial";
  c.fillText("FPS: " + fps.toFixed(2), 20, y + 30);

  c.fillText("Angle: " + scene.airfoilAngle.toFixed(1), 20, y + 60);
}

function setObstacle(x, y, reset) {
  scene.obstacleX = x;
  scene.obstacleY = y;

  setAirfoilObstacle(x, y, reset);

  scene.fluid.updateObstacle(obstacle);

  scene.showAirfoilBoundary = true;
}

function setAirfoilObstacle(x, y, reset) {
  var vx = 0.0;
  var vy = 0.0;

  obstacle.reset();

  if (!reset) {
    //vx = (x - scene.obstacleX) / scene.dt;
    //vy = (y - scene.obstacleY) / scene.dt;
  }

  var r = scene.obstacleRadius;
  var f = scene.fluid;
  var n = f.numY;
  var cd = Math.sqrt(2) * f.h;

  // set airfoil dimensions
  var ac = 0.6; // chord

  var ang = (scene.airfoilAngle * Math.PI) / 180;
  var cs = Math.cos(ang);
  var sn = Math.sin(ang);

  // generate line segments for top of airfoil
  var lp = new Vector(x, y);
  for (var i = 0; i < airfoilSteps; i++) {
    var xt = ac * airfoilUpperProfile[i].x;
    var yt = ac * airfoilUpperProfile[i].y;

    var rx = xt * cs - yt * sn;
    var ry = xt * sn + yt * cs;

    var x1 = x + rx;
    var y1 = y + ry;
    var p = new Vector(x1, y1);

    var line = new Line(lp, p);
    obstacle.addEdge(line);

    lp.x = x1;
    lp.y = y1;
  }
  // generate line segments for bottom of airfoil
  for (var i = 0; i < airfoilSteps; i++) {
    var xt = ac * airfoilLowerProfile[airfoilSteps - i - 1].x;
    var yt = ac * airfoilLowerProfile[airfoilSteps - i - 1].y;

    var rx = xt * cs - yt * sn;
    var ry = xt * sn + yt * cs;

    var x1 = x + rx;
    var y1 = y + ry;
    var p = new Vector(x1, y1);

    var line = new Line(lp, p);
    obstacle.addEdge(line);

    lp.x = x1;
    lp.y = y1;
  }

  obstacle.finalise();

  // for each cell in the fluid grid
  /*
  for (var i = 1; i < f.numX - 2; i++) {
    for (var j = 1; j < f.numY - 2; j++) {
      f.s[i * n + j] = 1.0;

      // determine coords relative to centre of obstacle
      var dx = (i + 0.5) * f.h - x;
      var dy = (j + 0.5) * f.h - y;

      // rotate
      var rx = dx * cs - dy * sn;
      var ry = dx * sn + dy * cs;

      // are we in the x-region that is potentially within the airfoil
      if (rx >= 0 && rx <= ac) {
        var ax = rx / ac;

        // lookup yt
        var ai = Math.round(airfoilSteps * ax);
        var yt = ac * airfoilProfile[ai];

        if (Math.abs(ry) <= yt) {
          // set s for obstacle
          f.s[i * n + j] = 0.0;
          f.m[i * n + j] = 1.0;
          // update u and v
          f.u[i * n + j] = vx;
          f.u[(i + 1) * n + j] = vx;
          f.v[i * n + j] = vy;
          f.v[i * n + j + 1] = vy;
        }
      }
    }
  }
  */
}

// interaction -------------------------------------------------------

var mouseDown = false;

function startDrag(x, y) {
  let bounds = canvas.getBoundingClientRect();

  let mx = x - bounds.left - canvas.clientLeft;
  let my = y - bounds.top - canvas.clientTop;
  mouseDown = true;

  x = mx / cScale;
  y = (canvas.height - my) / cScale;

  
  hoverX = Math.floor(x / f.h);
  hoverY = Math.floor(y / f.h);

  // DIAGNOSTICS
  if (hoverX >= 0 && hoverX < f.numX && hoverY >=0 && hoverY < f.numY) {
    var ci = hoverX * f.numY + hoverY;
    console.log('cellBounds',f.cellBounds[ci]);
    console.log('cellBounds right',f.cellBounds[(hoverX+1) * f.numY + hoverY]);
    console.log('cellBounds above',f.cellBounds[hoverX * f.numY + hoverY+1]);
    console.log('u',f.u[ci]);
    console.log('v',f.v[ci]);
    console.log('u+1',f.u[(hoverX+1) * f.numY + hoverY]);
    console.log('v+1',f.v[hoverX * f.numY + hoverY+1]);
    console.log('incompressibility',f.incompressibility[ci]);
  }

  //setObstacle(x,y, true);
}

function drag(x, y) {
  let bounds = canvas.getBoundingClientRect();
  let mx = x - bounds.left - canvas.clientLeft;
  let my = y - bounds.top - canvas.clientTop;
  x = mx / cScale;
  y = (canvas.height - my) / cScale;
  
  if (mouseDown) {
    //setObstacle(x,y, false);

    /*
    var dy = y - 0.5;
    airfoilAngle = 60 * dy;

    setObstacle(0.401, 0.501, true);
    */
  } 
}

function endDrag() {
  mouseDown = false;
}

canvas.addEventListener("mousedown", (event) => {
  startDrag(event.x, event.y);
});

canvas.addEventListener("mouseup", (event) => {
  endDrag();
});

canvas.addEventListener("mousemove", (event) => {
  drag(event.x, event.y);
});

canvas.addEventListener("touchstart", (event) => {
  startDrag(event.touches[0].clientX, event.touches[0].clientY);
});

canvas.addEventListener("touchend", (event) => {
  endDrag();
});

canvas.addEventListener(
  "touchmove",
  (event) => {
    event.preventDefault();
    event.stopImmediatePropagation();
    drag(event.touches[0].clientX, event.touches[0].clientY);
  },
  { passive: false }
);

document.addEventListener("keydown", (event) => {
  switch (event.key) {
    case "p":
      scene.paused = !scene.paused;
      break;
    case "m":
      scene.paused = false;
      simulate();
      scene.paused = true;
      break;
  }
});

function toggleStart() {
  var button = document.getElementById("startButton");
  if (scene.paused) button.innerHTML = "Stop";
  else button.innerHTML = "Start";
  scene.paused = !scene.paused;
}

// main -------------------------------------------------------

function simulate() {
  if (!scene.paused)
    scene.fluid.simulate(scene.dt, scene.gravity, scene.numIters);
  scene.frameNr++;
}

function update() {
  frame++;
  var loopTime = performance.now();
  if (loopTime > fpsTimer + 1000) {
    var dt = loopTime - fpsTimer;
    fps = frame / (fpsTimer / 1000);
    fpsTimer = loopTime;
  }
  simulate();
  draw();
  requestAnimationFrame(update);
}

function init() {

  setupAirfoil();
  setupScene(1);

  // setup dat.gui
  // Creating a GUI with options.
  gui = new dat.GUI({name: 'GUI'});

  var setupFolder = gui.addFolder('Setup');

  setupFolder.add(scene, 'airfoilAngle', -90, 90, 1).onChange(()=>{
    setObstacle(0.401, 0.501, true);
  });


  var visFolder = gui.addFolder('Visualisation');

  visFolder.add(scene, 'sceneNr', 1, 3, 1).onFinishChange(()=>{
    setupScene(scene.sceneNr);
  });

  visFolder.add(scene, 'showCells');

  visFolder.add(scene, 'showAirfoilBoundary');

  visFolder.add(scene, 'showStreamlines');

  visFolder.add(scene, 'showVelocities');

  visFolder.add(scene, 'showPressure');

  visFolder.add(scene, 'showPressureDistribution');

  visFolder.add(scene, 'showSmoke');

  visFolder.add(scene, 'overRelaxation', 0.5, 1.9, 0.1);

  update();
}

$(document).ready(function () {
  init();
});
