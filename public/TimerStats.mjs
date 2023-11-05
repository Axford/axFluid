export default class TimerStats {
    
    constructor(label) {
      this.v = 0;
      this.label = label;
      this.samples = 10;
    }
  
    update(t0, t1) {
      var t = t1 - t0;
      this.v = (this.v * (this.samples - 1) + t) / this.samples;
    }
  
    draw(c, x, y) {
      c.font = "14px serif";
      c.fillStyle = "#00f";
      c.fillText(this.label + ": " + this.v.toFixed(4) + " ms", x, y);
    }
  }