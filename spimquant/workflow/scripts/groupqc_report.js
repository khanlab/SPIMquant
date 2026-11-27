// Group QC report. DATA is embedded by groupqc_report.py (see build_data there).
(() => {
  const DATA = JSON.parse(document.getElementById("qc-data").textContent);
  const S = DATA.subjects;
  const $ = (sel, el = document) => el.querySelector(sel);
  const h = (tag, attrs = {}, ...kids) => {
    const el = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs)) {
      if (v === undefined || v === null || v === false) continue;
      if (k === "class") el.className = v;
      else if (k === "style") el.style.cssText = v;
      else if (k.startsWith("on")) el.addEventListener(k.slice(2), v);
      else el.setAttribute(k, v === true ? "" : v);
    }
    for (const kid of kids.flat()) {
      if (kid === null || kid === undefined || kid === false) continue;
      el.append(kid instanceof Node ? kid : document.createTextNode(kid));
    }
    return el;
  };
  const store = {
    get(k, d) { try { const v = localStorage.getItem("groupqc:" + k); return v === null ? d : JSON.parse(v); } catch { return d; } },
    set(k, v) { try { localStorage.setItem("groupqc:" + k, JSON.stringify(v)); } catch { /* storage unavailable */ } },
  };

  const planes = Object.keys(DATA.planes);
  const layerByKey = Object.fromEntries(DATA.layers.map((l) => [l.key, l]));
  const metricByKey = Object.fromEntries(DATA.metrics.map((m) => [m.key, m]));
  const regLayer = DATA.layers.find((l) => l.kind === "gray");

  const state = {
    view: "subject",
    cur: 0,
    plane: planes[0],
    slice: Object.fromEntries(planes.map((p) => [p, Math.floor(DATA.planes[p].positions.length / 2)])),
    layer: store.get("layer", regLayer ? regLayer.key : DATA.layers[0].key),
    outlines: store.get("outlines", true),
    shared: store.get("shared", true),
    template: false,
    sort: "label",
    sortDir: 1,
    filter: "",
    groupFilter: "",
    tile: store.get("tile", 210),
    fps: store.get("fps", 4),
    playing: false,
  };
  if (!layerByKey[state.layer]) state.layer = DATA.layers[0].key;

  // ---------- formatting ----------
  const fmt = (v, m) => {
    if (v === undefined || v === null || Number.isNaN(v)) return "–";
    const a = Math.abs(v);
    let s;
    if (a !== 0 && (a >= 1e5 || a < 1e-3)) s = v.toExponential(2);
    else if (m && m.unit === "%") s = v.toFixed(a >= 10 ? 1 : 2);
    else if (a >= 100) s = v.toFixed(0);
    else if (a >= 10) s = v.toFixed(1);
    else s = v.toFixed(m && m.unit === "" && a < 1 ? 3 : 2);
    return m && m.unit ? `${s} ${m.unit}` : s;
  };
  const flagCount = (s) => Object.keys(s.flags).length;

  // ---------- group colors ----------
  const groupCol = DATA.group_by[0] || null;
  const groupVals = groupCol ? [...new Set(S.map((s) => s.cov[groupCol] ?? "n/a"))].sort() : [];
  const PALETTE = ["#2a7fb8", "#d9822b", "#3aa17e", "#b9488f", "#7d6bc4", "#c2a634", "#5d8a9b", "#c4534a"];
  const groupColor = (s) => (groupCol ? PALETTE[groupVals.indexOf(s.cov[groupCol] ?? "n/a") % PALETTE.length] : "var(--dot)");

  // ---------- order (filter + sort) ----------
  let order = [];
  function computeOrder() {
    const q = state.filter.trim().toLowerCase();
    let idx = S.map((_, i) => i).filter((i) => {
      const s = S[i];
      if (state.groupFilter && groupCol && String(s.cov[groupCol] ?? "n/a") !== state.groupFilter) return false;
      if (!q) return true;
      return (s.label + " " + Object.values(s.cov).join(" ")).toLowerCase().includes(q);
    });
    const key = state.sort;
    const val = (i) => {
      const s = S[i];
      if (key === "label") return s.label;
      if (key === "flags") return -flagCount(s);
      if (key.startsWith("cov:")) return s.cov[key.slice(4)] ?? "";
      return s.metrics[key] ?? null;
    };
    idx.sort((a, b) => {
      const va = val(a), vb = val(b);
      if (va === null && vb === null) return 0;
      if (va === null) return 1;
      if (vb === null) return -1;
      const c = typeof va === "number" && typeof vb === "number" ? va - vb : String(va).localeCompare(String(vb), undefined, { numeric: true });
      return c * state.sortDir || S[a].label.localeCompare(S[b].label, undefined, { numeric: true });
    });
    order = idx;
  }

  // ---------- images ----------
  // colormap for value layers: inferno-like, transparent at the low end
  const STOPS = [[0, [0, 0, 4]], [0.25, [87, 16, 110]], [0.5, [188, 55, 84]], [0.75, [249, 142, 9]], [1, [252, 255, 164]]];
  const LUT = new Uint8ClampedArray(256 * 4);
  for (let i = 0; i < 256; i++) {
    const t = i / 255;
    let j = 0;
    while (j < STOPS.length - 2 && t > STOPS[j + 1][0]) j++;
    const [t0, c0] = STOPS[j], [t1, c1] = STOPS[j + 1];
    const f = (t - t0) / (t1 - t0);
    for (let c = 0; c < 3; c++) LUT[i * 4 + c] = c0[c] + f * (c1[c] - c0[c]);
    LUT[i * 4 + 3] = i === 0 ? 0 : Math.min(1, t / 0.12) * 235;
  }
  const lutCss = `linear-gradient(90deg, ${STOPS.map(([t, c]) => `rgb(${c.join(",")}) ${t * 100}%`).join(",")})`;

  const imgCache = new Map();
  function loadImg(src) {
    if (!imgCache.has(src)) {
      imgCache.set(src, new Promise((resolve, reject) => {
        const im = new Image();
        im.onload = () => resolve(im);
        im.onerror = reject;
        im.src = src;
      }));
    }
    return imgCache.get(src);
  }
  const colorCache = new Map();
  function colorize(src, ratio) {
    // ratio: stored scale / display scale, so 255 * ratio is the top of the display range
    const key = src.length + ":" + src.slice(-48) + ":" + ratio.toFixed(4);
    if (!colorCache.has(key)) {
      colorCache.set(key, loadImg(src).then((im) => {
        const c = document.createElement("canvas");
        c.width = im.naturalWidth; c.height = im.naturalHeight;
        const ctx = c.getContext("2d", { willReadFrequently: true });
        ctx.drawImage(im, 0, 0);
        const d = ctx.getImageData(0, 0, c.width, c.height);
        const px = d.data;
        for (let i = 0; i < px.length; i += 4) {
          const v = px[i];
          const k = v === 0 ? 0 : Math.max(1, Math.min(255, Math.round(v * ratio))) * 4;
          px[i] = LUT[k]; px[i + 1] = LUT[k + 1]; px[i + 2] = LUT[k + 2]; px[i + 3] = LUT[k + 3];
        }
        ctx.putImageData(d, 0, 0);
        return c.toDataURL();
      }));
    }
    return colorCache.get(key);
  }
  const displayScale = (layer, s) => (state.shared ? DATA.shared_scale[layer.key] : s.layers[layer.key].scale);

  // One slice of a layer for subject s, as an image stack
  function stack(s, layerKey, plane, k, opts = {}) {
    const layer = layerByKey[layerKey];
    const p = DATA.planes[plane];
    const el = h("div", { class: "stack", style: `aspect-ratio: ${p.aspect}` });
    const L = s && s.layers[layerKey];
    const showTemplate = opts.template || !L;
    let baseSrc;
    if (showTemplate) baseSrc = DATA.template.anat[plane][k];
    else if (layer.kind === "gray") baseSrc = L.slices[plane][k];
    else {
      const under = s.layers["stain:" + layer.stain];
      baseSrc = under ? under.slices[plane][k] : DATA.template.anat[plane][k];
    }
    el.append(h("img", { src: baseSrc, alt: "", decoding: "async", draggable: "false" }));
    if (!showTemplate && layer.kind === "value") {
      const ov = h("img", { class: "ol", alt: "", draggable: "false" });
      el.append(ov);
      colorize(L.slices[plane][k], L.scale / displayScale(layer, s)).then((u) => { ov.src = u; });
    }
    if (state.outlines && DATA.template.outlines) {
      el.append(h("img", { class: "ol", src: DATA.template.outlines[plane][k], alt: "", draggable: "false" }));
    }
    if (opts.caption) el.append(h("span", { class: "cap" }, opts.caption));
    if (showTemplate && !opts.noTag) el.append(h("span", { class: "tag" }, L ? "Template" : "No image"));
    return el;
  }

  // ---------- strip plots ----------
  function dots(key, curIdx) {
    // positions in percent of the width, so dots stay round at any width
    const NS = "http://www.w3.org/2000/svg";
    const svg = document.createElementNS(NS, "svg");
    svg.setAttribute("class", "dots");
    const vals = S.map((s, i) => [i, s.metrics[key]]).filter(([, v]) => v !== undefined && v !== null);
    if (!vals.length) return svg;
    const st = DATA.metric_stats[key];
    let lo = Math.min(...vals.map((x) => x[1])), hi = Math.max(...vals.map((x) => x[1]));
    if (hi === lo) { lo -= 1; hi += 1; }
    const pad = (hi - lo) * 0.04;
    const x = (v) => 2 + ((v - lo + pad) / (hi - lo + 2 * pad)) * 96;
    const pct = (v) => `${v.toFixed(3)}%`;
    const add = (tag, attrs) => {
      const e = document.createElementNS(NS, tag);
      for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
      svg.append(e);
      return e;
    };
    if (st && st.spread > 0) {
      const z = DATA.outlier_z;
      const a = Math.max(lo - pad, st.median - z * st.spread), b = Math.min(hi + pad, st.median + z * st.spread);
      add("rect", { class: "band", x: pct(x(a)), y: 3, width: pct(Math.max(0, x(b) - x(a))), height: 20 });
      add("line", { class: "med", x1: pct(x(st.median)), x2: pct(x(st.median)), y1: 3, y2: 23 });
    }
    add("line", { class: "axis", x1: "0%", x2: "100%", y1: 13, y2: 13 });
    // jitter by subject index so equal values stay visible
    for (const [i, v] of vals) {
      if (i === curIdx) continue;
      const c = add("circle", { class: "d" + (S[i].flags[key] ? " flag" : ""), cx: pct(x(v)), cy: 13 + (((i * 7) % 9) - 4), r: 3.2 });
      const t = document.createElementNS(NS, "title");
      t.textContent = `${S[i].label}: ${fmt(v, metricByKey[key])}`;
      c.append(t);
      c.addEventListener("click", () => go("subject", i));
    }
    const cv = S[curIdx] && S[curIdx].metrics[key];
    if (cv !== undefined && cv !== null) add("circle", { class: "cur", cx: pct(x(cv)), cy: 13, r: 6 });
    return svg;
  }

  function metricRows(keys, s, i) {
    const box = h("div", { class: "mrows" });
    for (const key of keys) {
      const m = metricByKey[key];
      if (!m) continue;
      const v = s.metrics[key];
      const z = s.z[key];
      const cls = s.flags[key] ? "bad" : z !== undefined && Math.abs(z) > DATA.outlier_z * 0.67 ? "warn" : "";
      box.append(
        h("span", { class: "ml", title: m.help }, m.label),
        h("span", { class: "mv " + cls, title: z !== undefined ? `robust z = ${z.toFixed(1)}` : "" }, fmt(v, m)),
        dots(key, i),
      );
    }
    return box;
  }

  // ---------- shared controls ----------
  function planeControl(onchange) {
    const seg = h("div", { class: "seg", role: "group", "aria-label": "Plane" });
    for (const p of planes) {
      seg.append(h("button", {
        type: "button", "aria-pressed": String(p === state.plane),
        onclick: () => { state.plane = p; onchange(); },
      }, p[0].toUpperCase() + p.slice(1)));
    }
    return seg;
  }
  function check(id, label, prop, onchange) {
    const input = h("input", { type: "checkbox", id });
    input.checked = state[prop];
    input.addEventListener("change", () => { state[prop] = input.checked; store.set(prop, input.checked); onchange(); });
    return h("label", { class: "check", for: id }, input, label);
  }
  function layerSelect(id, onchange) {
    const sel = h("select", { id, "aria-label": "Image" });
    const groups = {};
    for (const l of DATA.layers) (groups[l.group] ||= []).push(l);
    for (const [g, ls] of Object.entries(groups)) {
      const og = h("optgroup", { label: g });
      for (const l of ls) og.append(h("option", { value: l.key, selected: l.key === state.layer }, l.label));
      sel.append(og);
    }
    sel.addEventListener("change", () => { state.layer = sel.value; store.set("layer", sel.value); onchange(); });
    return sel;
  }
  function sliceControl(id, onchange) {
    const n = DATA.planes[state.plane].positions.length;
    const input = h("input", { type: "range", id, min: 0, max: n - 1, step: 1, value: state.slice[state.plane] });
    const out = h("output", { for: id }, DATA.planes[state.plane].positions[state.slice[state.plane]]);
    input.addEventListener("input", () => {
      state.slice[state.plane] = +input.value;
      out.textContent = DATA.planes[state.plane].positions[+input.value];
      onchange();
    });
    return h("label", { class: "ctl", for: id }, "Slice", input, out);
  }
  function sortSelect(id, onchange) {
    const sel = h("select", { id, "aria-label": "Order" });
    const opts = [["label", "Name"], ["flags", "Most flags first"]];
    for (const c of DATA.group_by) opts.push(["cov:" + c, c]);
    for (const m of DATA.metrics) opts.push([m.key, m.label + (m.stain ? "" : "")]);
    for (const [v, l] of opts) sel.append(h("option", { value: v, selected: v === state.sort }, "Order: " + l));
    sel.addEventListener("change", () => { state.sort = sel.value; state.sortDir = 1; computeOrder(); onchange(); });
    return sel;
  }
  function groupSelect(id, onchange) {
    if (!groupCol) return null;
    const sel = h("select", { id, "aria-label": groupCol });
    sel.append(h("option", { value: "" }, `All ${groupCol}`));
    for (const g of groupVals) sel.append(h("option", { value: g, selected: g === state.groupFilter }, `${groupCol}: ${g}`));
    sel.addEventListener("change", () => { state.groupFilter = sel.value; computeOrder(); onchange(); });
    return sel;
  }
  function legend(layer) {
    if (layer.kind !== "value") return null;
    const top = state.shared ? fmt(DATA.shared_scale[layer.key]) : "own max";
    return h("span", { class: "legend" }, "0", h("span", { class: "bar", style: `background:${lutCss}` }), `${top} ${layer.unit || ""}`);
  }

  // ---------- subject view ----------
  const sideList = $("#side-list");
  function renderSidebar() {
    sideList.replaceChildren(...order.map((i) => {
      const s = S[i];
      const n = flagCount(s);
      return h("li", {}, h("button", {
        type: "button", "aria-current": String(i === state.cur), onclick: () => go("subject", i),
        title: groupCol ? `${groupCol}: ${s.cov[groupCol] ?? "n/a"}` : "",
      }, h("span", { class: "gdot", style: `background:${groupColor(s)}` }), h("span", { class: "name" }, s.label), n ? h("span", { class: "badge", title: `${n} flagged metric(s)` }, String(n)) : null));
    }));
    const cur = $('[aria-current="true"]', sideList);
    if (cur) cur.scrollIntoView({ block: "nearest" });
  }

  function step(delta) {
    if (!order.length) return;
    let pos = order.indexOf(state.cur);
    pos = pos < 0 ? 0 : (pos + delta + order.length) % order.length;
    state.cur = order[pos];
  }

  function renderSubject() {
    const s = S[state.cur];
    const i = state.cur;
    const pos = order.indexOf(i);
    const root = $("#subject-main");
    const redraw = () => renderSubject();

    const head = h("div", { class: "subject-head" },
      h("h2", {}, s.label),
      h("span", { class: "pos" }, pos >= 0 ? `${pos + 1} of ${order.length}` : "not in the current filter"),
      h("button", { class: "btn", type: "button", onclick: () => { step(-1); update(); }, title: "Previous subject (←)" }, "← Previous"),
      h("button", { class: "btn primary", type: "button", onclick: () => { step(1); update(); }, title: "Next subject (→)" }, "Next →"),
    );
    const covs = h("div", { class: "covs" }, ...Object.entries(s.cov).map(([k, v]) => h("span", { class: "chip" }, `${k} `, h("b", {}, String(v)))));
    const flagged = Object.entries(s.flags);
    const flags = h("div", { class: "flags" },
      flagged.length
        ? flagged.map(([k, why]) => h("span", { class: "pill bad", title: metricByKey[k] ? metricByKey[k].help : "" }, `${metricByKey[k] ? metricByKey[k].label : k}: ${why}`))
        : h("span", { class: "pill good" }, "No metric flagged"));

    // view controls shared by the image sections
    const holdBtn = h("button", { class: "btn", type: "button", title: "Hold to show the template (or hold T)" }, "Hold for template");
    const setTpl = (on) => { if (state.template !== on) { state.template = on; redrawImages(); } };
    holdBtn.addEventListener("pointerdown", () => setTpl(true));
    for (const ev of ["pointerup", "pointerleave", "pointercancel"]) holdBtn.addEventListener(ev, () => setTpl(false));
    const controls = h("div", { class: "toolbar" },
      planeControl(redraw),
      check("s-outlines", "Atlas outlines", "outlines", () => redrawImages()),
      DATA.layers.some((l) => l.kind === "value") ? check("s-shared", "Shared color scale", "shared", () => redrawImages()) : null,
      holdBtn,
      h("label", { class: "ctl", for: "s-tile" }, "Size", (() => {
        const r = h("input", { type: "range", id: "s-tile", min: 120, max: 420, step: 10, value: state.tile });
        r.addEventListener("input", () => { state.tile = +r.value; store.set("tile", state.tile); document.documentElement.style.setProperty("--tile", state.tile + "px"); });
        return r;
      })()),
    );

    const strips = [];
    const strip = (layerKey) => {
      const el = h("div", { class: "strip" });
      strips.push([el, layerKey]);
      return el;
    };
    function redrawImages() {
      const p = DATA.planes[state.plane];
      for (const [el, key] of strips) {
        el.replaceChildren(...p.positions.map((label, k) => stack(s, key, state.plane, k, { caption: label, template: state.template })));
      }
    }

    const cards = [];
    // registration
    const reports = s.figures.filter((f) => f.section === "registration");
    if (regLayer) {
      cards.push(h("section", { class: "card", id: "sec-registration" },
        h("div", { class: "card-head" }, h("h2", {}, "Registration"),
          h("span", { class: "sub" }, `${regLayer.stain} in ${DATA.template_name} space`),
          h("span", { class: "spacer" }),
          ...reports.map((f) => h("a", { href: f.src, target: "_blank", rel: "noopener" }, "Open the subject's registration report ↗"))),
        controls,
        h("div", { style: "height:10px" }),
        strip(regLayer.key),
        metricRows(DATA.metrics.filter((m) => m.section === "registration").map((m) => m.key), s, i),
      ));
    }
    // stains
    const otherStains = DATA.layers.filter((l) => l.kind === "gray" && l !== regLayer);
    const intensityKeys = (stain) => DATA.metrics.filter((m) => m.section === "intensity" && m.stain === stain).map((m) => m.key);
    cards.push(h("section", { class: "card", id: "sec-stains" },
      h("div", { class: "card-head" }, h("h2", {}, "Stains in template space"), h("span", { class: "sub" }, "Each image is windowed on its own; intensities are compared in the metrics.")),
      ...DATA.layers.filter((l) => l.kind === "gray").map((l) => h("div", { class: "block" },
        h("div", { class: "block-title" }, h("h3", {}, l.stain), l === regLayer ? h("span", { class: "sub" }, "registration stain") : null),
        l === regLayer ? null : strip(l.key),
        metricRows(intensityKeys(l.stain), s, i))),
    ));
    if (!otherStains.length && !intensityKeys(regLayer && regLayer.stain).length) cards.pop();
    // segmentation
    const segLayers = DATA.layers.filter((l) => l.kind === "value");
    if (segLayers.length) {
      cards.push(h("section", { class: "card", id: "sec-segmentation" },
        h("div", { class: "card-head" }, h("h2", {}, "Segmentation"), h("span", { class: "sub" }, "Field fraction over the warped stain")),
        ...segLayers.map((l) => h("div", { class: "block" },
          h("div", { class: "block-title" }, h("h3", {}, l.label), h("span", { class: "sub" }, l.desc), legend(l)),
          strip(l.key),
          metricRows(DATA.metrics.filter((m) => m.section === "segmentation" && m.stain === l.stain && (!m.desc || m.desc === l.desc)).map((m) => m.key), s, i))),
      ));
    }
    // figures
    const figs = s.figures.filter((f) => f.section !== "registration");
    if (figs.length) {
      const bySec = {};
      for (const f of figs) (bySec[f.section] ||= []).push(f);
      cards.push(h("section", { class: "card", id: "sec-figures" },
        h("div", { class: "card-head" }, h("h2", {}, "Subject QC figures"), h("span", { class: "sub" }, DATA.embedded ? "Embedded in this report" : "Linked from the output folder: open the report where SPIMquant wrote it")),
        ...Object.entries(bySec).map(([sec, fs]) => h("div", { class: "fig-group" },
          h("h3", {}, sec),
          h("div", { class: "figs" }, ...fs.map(figCard))))));
    }

    root.replaceChildren(head, covs, flags, ...cards);
    redrawImages();
    renderSidebar();
  }

  function figCard(f) {
    const thumb = h("div", { class: "thumb" });
    if (/\.(png|gif|webp|jpe?g)$/i.test(f.path) || f.src.startsWith("data:")) {
      const im = h("img", { src: f.src, alt: f.label, loading: "lazy" });
      im.addEventListener("error", () => thumb.replaceChildren(h("span", { class: "missing" }, "Not found next to this report")));
      thumb.append(im);
    } else thumb.append(h("span", { class: "doc" }, "HTML report ↗"));
    return h("a", { class: "fig", href: f.src, target: "_blank", rel: "noopener", title: f.path },
      thumb, h("div", { class: "meta" }, h("b", {}, f.label), h("span", {}, f.detail)));
  }

  // ---------- montage view ----------
  function renderMontage() {
    const layer = layerByKey[state.layer];
    const controls = $("#montage-controls");
    const redraw = () => renderMontage();
    controls.replaceChildren(
      h("div", { class: "toolbar" },
        layerSelect("m-layer", redraw),
        planeControl(redraw),
        sliceControl("m-slice", () => drawTiles()),
        check("m-outlines", "Atlas outlines", "outlines", () => drawTiles()),
        layer.kind === "value" ? check("m-shared", "Shared color scale", "shared", redraw) : null,
        legend(layer),
      ),
      h("div", { class: "toolbar" },
        sortSelect("m-sort", redraw),
        groupSelect("m-group", redraw),
        h("label", { class: "ctl", for: "m-tile" }, "Tile size", (() => {
          const r = h("input", { type: "range", id: "m-tile", min: 110, max: 420, step: 10, value: state.tile });
          r.addEventListener("input", () => { state.tile = +r.value; store.set("tile", state.tile); document.documentElement.style.setProperty("--tile", state.tile + "px"); });
          return r;
        })()),
        h("span", { class: "ctl" }, `${order.length} of ${S.length} subjects · click a tile to open it`),
      ),
    );
    drawTiles();
  }
  function drawTiles() {
    const layer = layerByKey[state.layer];
    const k = state.slice[state.plane];
    const sortMetric = metricByKey[state.sort];
    $("#montage").replaceChildren(...order.map((i) => {
      const s = S[i];
      const n = flagCount(s);
      return h("button", { type: "button", class: "mtile" + (n ? " flagged" : ""), onclick: () => go("subject", i), title: n ? Object.entries(s.flags).map(([key, why]) => `${metricByKey[key] ? metricByKey[key].label : key}: ${why}`).join("\n") : s.label },
        stack(s, layer.key, state.plane, k, { noTag: false }),
        h("span", { class: "lab" },
          h("span", { class: "gdot", style: `background:${groupColor(s)}` }),
          h("span", { class: "name" }, s.label),
          sortMetric ? h("span", { class: "v" }, fmt(s.metrics[state.sort], sortMetric)) : null,
          n ? h("span", { class: "badge" }, String(n)) : null));
    }));
  }

  // ---------- flip view ----------
  let flipPos = 0, timer = null, withTemplate = store.get("flipTemplate", false);
  function flipFrames() {
    const frames = order.map((i) => ({ i }));
    if (withTemplate) frames.unshift({ i: -1 });
    return frames;
  }
  function renderFlip() {
    const layer = layerByKey[state.layer];
    const redraw = () => renderFlip();
    const tplId = "f-template";
    const tpl = h("input", { type: "checkbox", id: tplId });
    tpl.checked = withTemplate;
    tpl.addEventListener("change", () => { withTemplate = tpl.checked; store.set("flipTemplate", withTemplate); flipPos = 0; redraw(); });
    $("#flip-controls").replaceChildren(
      h("div", { class: "toolbar" },
        layerSelect("f-layer", redraw),
        planeControl(redraw),
        sliceControl("f-slice", () => drawFrame()),
        check("f-outlines", "Atlas outlines", "outlines", () => drawFrame()),
        layer.kind === "value" ? check("f-shared", "Shared color scale", "shared", redraw) : null,
        legend(layer),
      ),
      h("div", { class: "toolbar" },
        sortSelect("f-sort", () => { flipPos = 0; redraw(); }),
        groupSelect("f-group", () => { flipPos = 0; redraw(); }),
        h("label", { class: "check", for: tplId }, tpl, "Start each loop with the template"),
      ),
    );
    const frames = flipFrames();
    flipPos = Math.min(flipPos, Math.max(0, frames.length - 1));
    const scrub = h("input", { type: "range", id: "f-scrub", min: 0, max: Math.max(0, frames.length - 1), value: flipPos, "aria-label": "Subject" });
    scrub.addEventListener("input", () => { flipPos = +scrub.value; drawFrame(); });
    const ticks = h("div", { class: "ticks", "aria-hidden": "true" }, ...frames.map((f, p) => (f.i >= 0 && flagCount(S[f.i]) ? h("i", { style: `left:${frames.length > 1 ? (p / (frames.length - 1)) * 100 : 0}%` }) : null)));
    const fps = h("input", { type: "range", id: "f-fps", min: 1, max: 20, step: 1, value: state.fps });
    const fpsOut = h("output", { for: "f-fps" }, `${state.fps} / s`);
    fps.addEventListener("input", () => { state.fps = +fps.value; store.set("fps", state.fps); fpsOut.textContent = `${state.fps} / s`; if (state.playing) play(true); });
    const playBtn = h("button", { class: "btn primary", type: "button", id: "f-play", onclick: () => play(!state.playing) }, state.playing ? "Pause" : "Play");
    $("#flip-stage").replaceChildren(
      h("div", { id: "flip-frame", style: "width:100%;display:grid;justify-items:center" }),
      h("div", { class: "scrub" }, scrub, ticks),
      h("div", { class: "transport" },
        h("button", { class: "btn", type: "button", onclick: () => { flipPos = (flipPos - 1 + frames.length) % frames.length; drawFrame(); } }, "← Previous"),
        playBtn,
        h("button", { class: "btn", type: "button", onclick: () => { flipPos = (flipPos + 1) % frames.length; drawFrame(); } }, "Next →"),
        h("label", { class: "ctl", for: "f-fps" }, "Speed", fps, fpsOut),
        h("button", { class: "btn", type: "button", onclick: () => { const f = frames[flipPos]; if (f && f.i >= 0) go("subject", f.i); } }, "Open this subject"),
        h("span", { class: "ctl" }, h("kbd", {}, "Space"), " play · ", h("kbd", {}, "←"), h("kbd", {}, "→"), " step"),
      ),
    );
    // decode every frame ahead of playback
    for (const f of frames) if (f.i >= 0) stack(S[f.i], layer.key, state.plane, state.slice[state.plane]);
    drawFrame();
  }
  function drawFrame() {
    const frames = flipFrames();
    const box = $("#flip-frame");
    if (!box || !frames.length) return;
    const f = frames[flipPos];
    const layer = layerByKey[state.layer];
    const k = state.slice[state.plane];
    const s = f.i >= 0 ? S[f.i] : null;
    const st = stack(s, layer.key, state.plane, k, { template: !s, noTag: true });
    const flagsTxt = s && flagCount(s) ? ` · ${flagCount(s)} flagged` : "";
    st.append(h("span", { class: "big" }, s ? s.label : `Template ${DATA.template_name}`,
      h("small", {}, s ? `${flipPos + 1 - (withTemplate ? 1 : 0)} of ${order.length}${groupCol ? ` · ${s.cov[groupCol] ?? "n/a"}` : ""}${flagsTxt}` : "reference")));
    box.replaceChildren(st);
    const scrub = $("#f-scrub");
    if (scrub) scrub.value = flipPos;
  }
  function play(on) {
    state.playing = on;
    clearInterval(timer);
    const btn = $("#f-play");
    if (btn) btn.textContent = on ? "Pause" : "Play";
    if (on) {
      timer = setInterval(() => {
        const n = flipFrames().length;
        if (!n) return;
        flipPos = (flipPos + 1) % n;
        drawFrame();
      }, 1000 / state.fps);
    }
  }

  // ---------- metrics view ----------
  function renderMetrics() {
    const dists = $("#dists");
    dists.replaceChildren(...DATA.metrics.map((m) => {
      const n = S.filter((s) => s.flags[m.key]).length;
      return h("div", { class: "dist", title: m.help },
        h("div", { class: "t" }, h("b", {}, m.label), h("span", {}, n ? `${n} flagged` : `median ${fmt(DATA.metric_stats[m.key] ? DATA.metric_stats[m.key].median : null, m)}`)),
        dots(m.key, -1));
    }));
    const thead = h("tr", {}, h("th", { "aria-sort": state.sort === "label" ? "ascending" : null, onclick: () => sortBy("label") }, "Subject"),
      ...DATA.group_by.map((c) => h("th", { onclick: () => sortBy("cov:" + c), "aria-sort": state.sort === "cov:" + c ? "ascending" : null }, c)),
      h("th", { onclick: () => sortBy("flags") }, "Flags"),
      ...DATA.metrics.map((m) => h("th", { title: m.help, onclick: () => sortBy(m.key), "aria-sort": state.sort === m.key ? (state.sortDir > 0 ? "ascending" : "descending") : null },
        h("span", { class: "sec" }, m.section), m.label)));
    const rows = order.map((i) => {
      const s = S[i];
      return h("tr", { onclick: () => go("subject", i) },
        h("th", {}, s.label),
        ...DATA.group_by.map((c) => h("td", {}, String(s.cov[c] ?? ""))),
        h("td", { class: "num" + (flagCount(s) ? " bad" : "") }, String(flagCount(s))),
        ...DATA.metrics.map((m) => {
          const z = s.z[m.key];
          const cls = s.flags[m.key] ? "bad" : z !== undefined && Math.abs(z) > DATA.outlier_z * 0.67 ? "warn" : "";
          return h("td", { class: "num " + cls, title: z !== undefined ? `robust z = ${z.toFixed(1)}${s.flags[m.key] ? " · " + s.flags[m.key] : ""}` : "" }, fmt(s.metrics[m.key], m));
        }));
    });
    $("#metrics-table").replaceChildren(h("thead", {}, thead), h("tbody", {}, ...rows));
  }
  function sortBy(key) {
    if (state.sort === key) state.sortDir *= -1;
    else { state.sort = key; state.sortDir = 1; }
    computeOrder();
    renderMetrics();
  }

  // ---------- navigation ----------
  const VIEWS = ["subject", "montage", "flip", "metrics"];
  function go(view, i) {
    if (i !== undefined) state.cur = i;
    state.view = view;
    const hash = view === "subject" ? S[state.cur].anchor : view;
    if (location.hash.slice(1) !== hash) history.replaceState(null, "", "#" + hash);
    update();
    window.scrollTo({ top: 0 });
  }
  function update() {
    if (state.view !== "flip") play(false);
    for (const v of VIEWS) $("#view-" + v).hidden = v !== state.view;
    for (const b of document.querySelectorAll("#tabs button")) b.setAttribute("aria-pressed", String(b.dataset.view === state.view));
    if (state.view === "subject") {
      const hash = S[state.cur].anchor;
      if (location.hash.slice(1) !== hash) history.replaceState(null, "", "#" + hash);
      renderSubject();
    } else if (state.view === "montage") renderMontage();
    else if (state.view === "flip") renderFlip();
    else renderMetrics();
  }
  function fromHash() {
    const name = decodeURIComponent(location.hash.slice(1));
    if (VIEWS.includes(name) && name !== "subject") state.view = name;
    else {
      const i = S.findIndex((s) => s.anchor === name);
      state.view = "subject";
      if (i >= 0) state.cur = i;
    }
  }

  // ---------- boot ----------
  document.documentElement.style.setProperty("--tile", state.tile + "px");
  for (const b of document.querySelectorAll("#tabs button")) b.addEventListener("click", () => go(b.dataset.view));
  const search = $("#side-search");
  search.addEventListener("input", () => { state.filter = search.value; computeOrder(); renderSidebar(); });
  const sideGroup = groupSelect("side-group", () => { renderSidebar(); renderSubject(); });
  if (sideGroup) search.after(sideGroup);
  const sideSort = sortSelect("side-sort", () => { renderSidebar(); renderSubject(); });
  search.after(sideSort);
  $("#theme").addEventListener("click", () => {
    const root = document.documentElement;
    const dark = root.dataset.theme ? root.dataset.theme === "dark" : matchMedia("(prefers-color-scheme: dark)").matches;
    root.dataset.theme = dark ? "light" : "dark";
  });
  document.addEventListener("keydown", (e) => {
    if (e.target.closest("input[type=search], select, textarea") || e.metaKey || e.ctrlKey || e.altKey) return;
    if (e.target.matches("input[type=range]") && (e.key === "ArrowLeft" || e.key === "ArrowRight")) return;
    if (state.view === "subject" && (e.key === "ArrowRight" || e.key === "j")) { step(1); update(); e.preventDefault(); }
    else if (state.view === "subject" && (e.key === "ArrowLeft" || e.key === "k")) { step(-1); update(); e.preventDefault(); }
    else if (state.view === "subject" && e.key === "t" && !e.repeat && !state.template) { state.template = true; renderSubject(); }
    else if (state.view === "flip" && e.key === " ") { play(!state.playing); e.preventDefault(); }
    else if (state.view === "flip" && (e.key === "ArrowRight" || e.key === "ArrowLeft")) {
      const n = flipFrames().length;
      flipPos = (flipPos + (e.key === "ArrowRight" ? 1 : -1) + n) % n;
      drawFrame();
      e.preventDefault();
    }
  });
  document.addEventListener("keyup", (e) => {
    if (state.view === "subject" && e.key === "t" && state.template) { state.template = false; renderSubject(); }
  });
  window.addEventListener("hashchange", () => { fromHash(); update(); });
  computeOrder();
  fromHash();
  update();
})();
