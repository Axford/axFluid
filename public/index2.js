import Vector from "./Vector.mjs";
import Line from "./Line.mjs";
import TimerStats from "./TimerStats.mjs";
import Fluid from "./Fluid2.mjs";
import Polygon from "./Polygon.mjs";

var canvas = document.getElementById("myCanvas");
var c = canvas.getContext("2d", { willReadFrequently: true });
canvas.width = window.innerWidth - 20;
canvas.height = window.innerHeight - 100;

canvas.focus();

var gui;

var f;

var airfoilSteps = 200;
var airfoilUpperProfile = [];
var airfoilLowerProfile = [];

// obstacle is a polygon
var obstacle = new Polygon();

var frame = 0;
var fps = 0;
var fpsTimer = 0;

var mousePos = new Vector(0,0);  
var mouseFlow = new Vector(0,0); // u, v

var scene = {
  simHeight: 1.2,
  gravity: 0,
  dt: 1.0 / 120.0,
  numIters: 100,
  frameNr: 0,
  overRelaxation: 1.9,
  obstacleX: 0.0,
  obstacleY: 0.0,
  obstacleRadius: 0.15,
  paused: true,
  sceneNr: 0,
  showAirfoilBoundary: false,
  showStreamlines: false,
  showVelocities: true,
  showPressure: false,
  showSmoke: true,
  showPressureDistribution:false,
  showCells:true,
  fluid: null,
  airfoilAngle: 0,
  cOffsetX: 0,
  cOffsetY: 0
};
scene.cScale = canvas.height / scene.simHeight;
scene.simWidth = canvas.width / scene.cScale;



function cX(x) {
  return (x + scene.cOffsetX) * scene.cScale;
}

function cY(y) {
  return canvas.height - (y + scene.cOffsetY) * scene.cScale;
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

function setupScene(sceneNr = 0) {
  console.log("setupScene", sceneNr);

  scene.sceneNr = sceneNr;
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
  var domainWidth = (domainHeight / scene.simHeight) * scene.simWidth;
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
        var cell = f.cells[i*n + j];
        var s = 1.0; // fluid
        if (i == 0 || j == 0 || j == f.numY - 1) s = 0.0; // solid

        if (i==4 && j == 3) s = 0;
        cell.s = s;

        if (i==1) {
          cell.inU = inVel;
          cell.u = inVel;
        }

        if (i==f.numX-1) {
          cell.sink = true;
        }
      }
    }

    // smoke injection
    var pipeH = 0.3 * f.numY;
    var minJ = Math.floor(0.5 * f.numY - 0.5 * pipeH);
    var maxJ = Math.floor(0.5 * f.numY + 0.5 * pipeH);

    for (var j = 0; j < f.numY; j++) {
      if (j % 20 < 4) f.cells[j].m = 0;
    }

    setObstacle(0.401, 0.501, true);


    if (sceneNr > 1 ) {
      scene.dt = 1.0 / 120.0;
      scene.numIters = 100;
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

  var minP = f.cells[0].p;
  var maxP = f.cells[0].p;

  var minQ = f.cells[0].incompressibility;
  var maxQ = f.cells[0].incompressibility;

  for (var i = 0; i < f.numCells; i++) {
    minP = Math.min(minP, f.cells[i].p);
    maxP = Math.max(maxP, f.cells[i].p);

    minQ = Math.min(minQ, f.cells[i].incompressibility);
    maxQ = Math.max(maxQ, f.cells[i].incompressibility);
  }

  if (minP > -1) minP = -1;
  if (maxP < 1) maxP = 1;

  //minQ = -2;
  //maxQ = 2;

  /*
  var id = c.getImageData(0, 0, canvas.width, canvas.height);

  var color = [255, 255, 255, 255];

  for (var i = 0; i < f.numX; i++) {
    for (var j = 0; j < f.numY; j++) {
      if (f.cells[i * n + j].s == 0) {
        color[0] = 80;
        color[1] = color[0];
        color[2] = color[0];
      } else if (scene.showPressure) {
        var p = f.cells[i * n + j].p;
        var s = f.cells[i * n + j].m;
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
        var s = f.cells[i * n + j].s;
        var q = f.cells[i * n + j].incompressibility;
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
  */

  for (var i = 0; i < f.numX; i++) {
    for (var j = 0; j < f.numY; j++) {
      var ci = i * n + j;
      var cell = f.cells[ci];

      cell.drawPressureFill(c, cX, cY, getSciColor, minP, maxP);
    }
  }

  
  // show grid and ua / va boundaries
  if (scene.showCells) {
    for (var i = 0; i < f.numX; i++) {
      for (var j = 0; j < f.numY; j++) {
        var ci = i * n + j;
        var cell = f.cells[ci];

        cell.drawBoundary(c, cX, cY);
      }
    }
  }
  

  if (scene.showVelocities) {
    var scale = 0.02;

    for (var i = 0; i < f.numX; i++) {
      for (var j = 0; j < f.numY; j++) {
        var ci = i * n + j;
        f.cells[ci].drawVelocities(c, cX, cY, scale);
      }
    }
  }
  

  if (scene.showStreamlines) {
    var segLen = f.h * 0.2;
    var xStepSize = Math.max(Math.floor(f.numX / 15), 1);
    var yStepSize = Math.max(Math.floor(f.numY / 30), 1);
    var numSegs = Math.floor(xStepSize * 5);
    
    for (var i = 1; i < f.numX - 1; i += xStepSize) {
      for (var j = 1; j < f.numY - 1; j += yStepSize) {
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
          if (l > 0) {
            x += (u / l) * segLen;
            y += (v / l) * segLen;
          }



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
      var l = f.h / 4;
      v.multiply(l);

      // sample pressure at p0 offset by v
      var cell = f.getCellAt(line.p0.x+v.x, line.p0.y+v.y);
      var p = cell.p

      // scale pressure line based on maxQ
      var p = 0.3 * (p / maxP);
      var pSign = p > 0;
      p = Math.abs(p);

      v.multiply(p/l);
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


  if (f.hoverCell) f.hoverCell.drawDiagnostics(c, cX, cY);

  var s = "pressure: " + minP.toFixed(0) + " ... " + maxP.toFixed(0) + " N/m";
  c.fillStyle = "#000000";
  c.font = "16px Arial";
  c.fillText(s, 10, 35);


  var s = "incompressibility: " + minQ.toFixed(2) + " ... " + maxQ.toFixed(2);
  c.fillStyle = "#000000";
  c.font = "16px Arial";
  c.fillText(s, 300, 35);


  var y = 130;
  for (const [key, t] of Object.entries(scene.fluid.timers)) {
    //t.draw(c, 20, y);
    y += 20;
  }

  // draw fps
  c.fillStyle = "#00f";
  c.font = "16px Arial";
  c.fillText("FPS: " + fps.toFixed(2), 20, y + 30);

  c.fillText("Angle: " + scene.airfoilAngle.toFixed(1), 20, y + 60);

  // draw error points
  c.fillStyle = "#f00";
  c.strokeStyle = '#fff';
  f.errorPoints.forEach((p)=>{
    c.beginPath();
    c.arc(cX(p.x), cY(p.y), 3, 0, 2*Math.PI);
    c.fill();
    c.stroke();
  });

  // sample mouse flow
  mouseFlow.x = f.sampleField(mousePos.x, mousePos.y, 'u');
  mouseFlow.y = f.sampleField(mousePos.x, mousePos.y, 'v');

  // draw mouse flow
  c.strokeStyle = "#000";
  c.lineWidth = 2;
  c.beginPath();
  var x0 = cX(mousePos.x);
  var x1 = cX(mousePos.x + mouseFlow.x * scale);
  var y0 = cY(mousePos.y);
  var y1 = cY(mousePos.y + mouseFlow.y * scale);

  // u
  c.moveTo(x0, y0);
  c.lineTo(x1, y0);
  // v
  c.moveTo(x0, y0);
  c.lineTo(x0, y1);

  c.stroke();
}

function setObstacle(x, y, reset) {
  scene.obstacleX = x;
  scene.obstacleY = y;

  setAirfoilObstacle(x, y, reset);

  scene.fluid.updateObstacle(obstacle);
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

  x = mx / scene.cScale - scene.cOffsetX;
  y = (canvas.height - my) / scene.cScale - scene.cOffsetY;

  // DIAGNOSTICS
  var cell = f.getCellAt(x,y);
  if (cell) {
    f.hoverCell = cell;
    console.log('cell',x,y, cell);
  }

  //setObstacle(x,y, true);
}

function drag(x, y) {
  let bounds = canvas.getBoundingClientRect();
  let mx = x - bounds.left - canvas.clientLeft;
  let my = y - bounds.top - canvas.clientTop;
  mousePos.x = mx / scene.cScale - scene.cOffsetX;
  mousePos.y = (canvas.height - my) / scene.cScale - scene.cOffsetY;
  
  if (mouseDown) {
    //setObstacle(x,y, false);

    /*
    var dy = y - 0.5;
    airfoilAngle = 60 * dy;

    setObstacle(0.401, 0.501, true);
    */
  } else {
    // analyse what's under the cursor
    try {
      var cell = f.getCellAt(mousePos.x, mousePos.y);
      
      //console.log(u, cell);

    } catch(e) {
      console.error(e, x, y, mousePos, mouseFlow);
    }
    

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
    case 's':  
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
    f.reset();
    setObstacle(0.401, 0.501, true);
  });

  setupFolder.add(scene, 'sceneNr', 1, 3, 1).onFinishChange(()=>{
    setupScene(scene.sceneNr);
  });


  var visFolder = gui.addFolder('Visualisation');

  visFolder.add(scene, 'showCells');

  visFolder.add(scene, 'showAirfoilBoundary');

  visFolder.add(scene, 'showStreamlines');

  visFolder.add(scene, 'showVelocities');

  visFolder.add(scene, 'showPressure');

  visFolder.add(scene, 'showPressureDistribution');

  visFolder.add(scene, 'showSmoke');

  visFolder.add(scene, 'overRelaxation', 0.5, 1.9, 0.1);

  var viewFolder = gui.addFolder('View');
  viewFolder.add(scene, 'cOffsetX',-2,2,0.01);
  viewFolder.add(scene, 'cOffsetY',-2,2,0.01);
  viewFolder.add(scene, 'cScale',400,5000,10);

  update();
}

$(document).ready(function () {
  init();
});
