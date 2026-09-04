// Replay fig1_spec.json as native, editable PowerPoint shapes. No images.
const fs = require('fs');
const P  = require('pptxgenjs');

const spec = JSON.parse(fs.readFileSync(process.argv[2] || 'fig1_spec.json', 'utf8'));
const OUT  = process.argv[3] || 'fig1_lean_editable.pptx';
const MM   = 25.4;
const inch = v => v / MM;
const hex  = c => c.replace('#', '').toUpperCase();

// Blend a colour toward white; PowerPoint line colours cannot carry alpha, so a
// de-emphasised stroke is baked rather than made translucent.
function fade(c, a) {
  if (a >= 0.999) return hex(c);
  const n = parseInt(hex(c), 16);
  const mix = k => Math.round(((n >> k) & 255) * a + 255 * (1 - a));
  return [mix(16), mix(8), mix(0)]
    .map(v => v.toString(16).padStart(2, '0')).join('').toUpperCase();
}
const bbox = pts => ({
  x0: Math.min(...pts.map(p => p[0])), y0: Math.min(...pts.map(p => p[1])),
  x1: Math.max(...pts.map(p => p[0])), y1: Math.max(...pts.map(p => p[1])),
});

const pres = new P();
pres.defineLayout({ name: 'FIG', width: inch(spec.w_mm), height: inch(spec.h_mm) });
pres.layout = 'FIG';
pres.author = 'Fig. 1 — multiplex hypergraph EBCM';
const s = pres.addSlide();
s.background = { color: 'FFFFFF' };

const line = (a, b, opt) => {
  const [x0, y0] = a, [x1, y1] = b;
  s.addShape(pres.ShapeType.line, {
    x: inch(Math.min(x0, x1)), y: inch(Math.min(y0, y1)),
    w: inch(Math.abs(x1 - x0)), h: inch(Math.abs(y1 - y0)),
    flipH: x1 < x0, flipV: y1 < y0, line: opt,
  });
};

for (const o of spec.shapes) {
  switch (o.kind) {

    case 'blob': {                                   // hyperedge, closed freeform
      const b = bbox(o.pts);
      s.addShape(pres.ShapeType.custGeom, {
        x: inch(b.x0), y: inch(b.y0), w: inch(b.x1 - b.x0), h: inch(b.y1 - b.y0),
        fill: { color: hex(o.fill), transparency: Math.round((1 - o.alpha) * 100) },
        line: { color: fade(o.edge, o.edge_alpha), width: o.lw,
                dashType: o.dashed ? 'dash' : 'solid' },
        points: o.pts.map(p => ({ x: inch(p[0] - b.x0), y: inch(p[1] - b.y0) }))
                     .concat([{ close: true }]),
      });
      break;
    }

    case 'node':
      s.addShape(pres.ShapeType.ellipse, {
        x: inch(o.x - o.r), y: inch(o.y - o.r), w: inch(2 * o.r), h: inch(2 * o.r),
        fill: { color: fade(o.fill, o.alpha) },
        line: { color: fade(o.edge, o.alpha), width: o.lw },
      });
      break;

    case 'arrow': {
      const head = { color: hex(o.color), width: o.lw, endArrowType: 'triangle',
                     dashType: o.dashed ? 'dash' : 'solid' };
      if (o.alpha < 0.999) head.color = fade(o.color, o.alpha);
      if (o.ctl) {                                   // quadratic Bezier
        const pts = [o.p0, o.ctl, o.p1], b = bbox(pts);
        s.addShape(pres.ShapeType.custGeom, {
          x: inch(b.x0), y: inch(b.y0), w: inch(b.x1 - b.x0), h: inch(b.y1 - b.y0),
          line: head,
          points: [
            { x: inch(o.p0[0] - b.x0), y: inch(o.p0[1] - b.y0) },
            { x: inch(o.p1[0] - b.x0), y: inch(o.p1[1] - b.y0),
              curve: { type: 'quadratic',
                       x1: inch(o.ctl[0] - b.x0), y1: inch(o.ctl[1] - b.y0) } },
          ],
        });
      } else line(o.p0, o.p1, head);
      break;
    }

    case 'rule':
      line(o.p0, o.p1, { color: hex(o.color), width: o.lw,
                         dashType: o.dashed ? 'dash' : 'solid' });
      break;

    case 'cross':
      line([o.x - o.s, o.y - o.s], [o.x + o.s, o.y + o.s], { color: hex(o.color), width: o.lw });
      line([o.x - o.s, o.y + o.s], [o.x + o.s, o.y - o.s], { color: hex(o.color), width: o.lw });
      break;

    case 'roundbox':
      s.addShape(pres.ShapeType.roundRect, {
        x: inch(o.x), y: inch(o.y), w: inch(o.w), h: inch(o.h),
        rectRadius: inch(o.r), fill: { color: hex(o.fill) },
        line: { color: hex(o.edge), width: o.lw },
      });
      break;

    case 'text': {
      const PAD = 1.2, w = o.w + 2 * PAD, h = o.h + 1.4;
      const x = o.ha === 'center' ? o.x - w / 2 : o.ha === 'right' ? o.x - w : o.x;
      const y = o.va === 'top' ? o.y : o.va === 'bottom' ? o.y - h : o.y - h / 2;
      s.addText(o.runs.map(r => ({
        text: r.t,
        options: { subscript: r.sub, superscript: r.sup,
                   italic: r.i || o.style === 'italic',
                   bold: r.b || o.weight === 'bold' },
      })), {
        x: inch(x), y: inch(y), w: inch(w), h: inch(h),
        fontSize: o.size, fontFace: 'Arial', color: hex(o.color),
        align: o.ha, valign: 'middle', margin: 0, isTextBox: true, wrap: false,
      });
      break;
    }
  }
}

pres.writeFile({ fileName: OUT }).then(() =>
  console.log(`wrote ${OUT}  (${spec.shapes.length} native shapes, ${spec.w_mm}x${spec.h_mm} mm)`));
