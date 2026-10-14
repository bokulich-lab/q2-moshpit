"use strict";

const tfa = JSON.parse(document.getElementById("tfa-data").textContent);
const showUnclassified = document.getElementById("show-unclassified");
const levelSelect = document.getElementById("level");
const taxonSearch = document.getElementById("taxon-search");
const selectedTaxonomy = document.getElementById("selected-taxonomy");
const taxonSelect = document.getElementById("taxon");
const functionSelect = document.getElementById("function");
const groupSelect = document.getElementById("group");
const metricSelect = document.getElementById("metric");
const limitSelect = document.getElementById("limit");
const status = document.getElementById("status");
const schema = "https://vega.github.io/schema/vega-lite/v5.json";
let activeViews = [];
let renderQueue = Promise.resolve();

function addOptions(select, labels, firstLabel) {
  const first = document.createElement("option");
  first.value = "";
  first.textContent = firstLabel;
  select.appendChild(first);
  for (const label of labels) {
    const option = document.createElement("option");
    option.value = label;
    option.textContent = label;
    select.appendChild(option);
  }
}

let current = tfa;
const functions = [...new Set(tfa.pairs.map(pair => pair.function))].sort();
let loadsByPair = [];
const depth = tfa.pairs.reduce((depth, pair) => Math.max(depth, pair.lineage.length), 0);
levelSelect.add(new Option("Full taxonomy", "0"));
for (let level = 1; level <= depth; level++) {
  levelSelect.add(new Option(`Level ${level}`, String(level)));
}
levelSelect.value = "0";

function pairVisible(pair) {
  return showUnclassified.checked || !isUnclassifiedTaxon(pair, Number(levelSelect.value));
}

function pairLabel(pair) {
  return pair.taxon_label || latestTaxonLabel(pair.lineage || [], pair.taxon_id);
}

function refreshTaxa() {
  const previous = taxonSelect.value;
  const query = taxonSearch.value.trim().toLocaleLowerCase();
  const taxa = [...new Set(current.pairs
    .filter(pairVisible)
    .filter(pair => !functionSelect.value || pair.function === functionSelect.value)
    .map(pair => pair.taxon))]
    .filter(taxon => taxon.toLocaleLowerCase().includes(query)).sort();
  taxonSelect.replaceChildren();
  addOptions(taxonSelect, taxa, query ? "Select a matching taxon" : "All taxa");
  const labelByTaxon = new Map(current.pairs.map(pair => [pair.taxon, pairLabel(pair)]));
  const terminal = taxon => labelByTaxon.get(taxon);
  const counts = new Map();
  for (const taxon of taxa) counts.set(terminal(taxon), (counts.get(terminal(taxon)) || 0) + 1);
  for (const option of taxonSelect.options) {
    if (!option.value) continue;
    const label = terminal(option.value);
    const parent = option.value.split(";").slice(0, -1).join("; ").trim();
    option.textContent = levelSelect.value === "0" ? option.value : counts.get(label) > 1 ? `${label} [${parent}]` : label;
    option.title = option.value;
  }
  taxonSelect.value = taxa.includes(previous) ? previous : "";
  showSelectedTaxonomy();
}

function showSelectedTaxonomy() {
  taxonSelect.title = taxonSelect.value || "All taxa";
  selectedTaxonomy.textContent = taxonSelect.value
    ? `Selected lineage: ${taxonSelect.value}`
    : taxonSearch.value.trim() ? "Search filters the taxon list. Select a taxon to filter the plots." : "";
}

function changeLevel() {
  current = collapseTaxonomy(tfa, Number(levelSelect.value));
  loadsByPair = current.pairs.map(() => []);
  for (const load of current.loads) loadsByPair[load[0]].push(load);
  taxonSelect.value = "";
  refreshTaxa();
}
changeLevel();
addOptions(functionSelect, functions, "All functions");
addOptions(groupSelect, Object.keys(tfa.groups), "All samples");

function matchedPairs() {
  const indices = [];
  current.pairs.forEach((pair, index) => {
    if (pairVisible(pair) && (!taxonSelect.value || pair.taxon === taxonSelect.value) &&
        (!functionSelect.value || pair.function === functionSelect.value)) {
      indices.push(index);
    }
  });
  return indices;
}

function sampleValues(pairIndices) {
  const values = new Float64Array(tfa.samples.length);
  for (const pairIndex of pairIndices) {
    for (const [, sampleIndex, load] of loadsByPair[pairIndex]) {
      values[sampleIndex] += load;
    }
  }
  const grouping = tfa.groups[groupSelect.value];
  return tfa.samples.map((sample, index) => ({
    sample,
    group: grouping ? grouping[index] : "All samples",
    load: values[index]
  }));
}

function plotWidth(id) {
  const gutter = {heatmap: 170, boxplot: 90, samples: 180}[id];
  return Math.max(240, document.getElementById(id).parentElement.clientWidth - gutter);
}

function heatmapSpec(indices) {
  const field = {sum: "total", mean: "mean", median: "median"}[metricSelect.value];
  const label = {sum: "Sum", mean: "Mean", median: "Median"}[metricSelect.value];
  const top = indices
    .map(index => current.pairs[index])
    .sort((a, b) => b[field] - a[field])
    .slice(0, Number(limitSelect.value));
  const taxaOrder = [...new Set(top.map(pair => pair.taxon))];
  const functionOrder = [...new Set(top.map(pair => pair.function))];
  const axisLabels = Object.fromEntries(top.map(pair => [pair.taxon, pairLabel(pair)]));
  return {
    $schema: schema,
    data: {values: top},
    width: Math.min(plotWidth("heatmap"), Math.max(300, functionOrder.length * 52)),
    height: {step: 29},
    mark: {type: "rect", cursor: "pointer", stroke: "white", strokeWidth: 1},
    encoding: {
      x: {field: "function", type: "nominal", sort: functionOrder, title: "Function"},
      y: {
        field: "taxon", type: "nominal", sort: taxaOrder, title: "Taxon",
        axis: {labelLimit: 0, labelExpr: `${JSON.stringify(axisLabels)}[datum.label]`}
      },
      color: {field, type: "quantitative", title: `${label} load`, scale: {scheme: "blues", zero: true}},
      tooltip: [
        {field: "taxon", type: "nominal", title: "Taxon"},
        {field: "taxon_id", type: "nominal", title: "Taxon ID"},
        {field: "function", type: "nominal", title: "Function"},
        {field, type: "quantitative", title: `${label} load`, format: ".4g"}
      ]
    },
    config: {view: {stroke: null}}
  };
}

function sampleGroupColor(samples) {
  // Use every sample, including categories absent from the top 100 bars.
  return {
    field: "group", type: "nominal", legend: null,
    scale: {domain: [...new Set(samples.map(sample => sample.group))].sort(), scheme: "tableau10"}
  };
}

function boxplotSpec(samples) {
  return {
    $schema: schema,
    data: {values: samples},
    width: plotWidth("boxplot"),
    height: 300,
    mark: {type: "boxplot", extent: 1.5},
    encoding: {
      x: {field: "group", type: "nominal", sort: sampleGroupColor(samples).scale.domain, title: groupSelect.value || "Samples"},
      y: {field: "load", type: "quantitative", title: "TFA load", scale: {zero: true}},
      color: sampleGroupColor(samples)
    },
    config: {view: {stroke: null}}
  };
}

function sampleSpec(samples) {
  const top = [...samples].sort((a, b) => b.load - a.load).slice(0, 100);
  return {
    $schema: schema,
    data: {values: top},
    width: plotWidth("samples"),
    height: 300,
    mark: {type: "bar"},
    encoding: {
      x: {field: "sample", type: "nominal", sort: "-y", title: "Sample", axis: {labels: top.length <= 30}},
      y: {field: "load", type: "quantitative", title: "TFA load", scale: {zero: true}},
      color: sampleGroupColor(samples),
      tooltip: [
        {field: "sample", type: "nominal", title: "Sample"},
        {field: "group", type: "nominal", title: groupSelect.value || "Group"},
        {field: "load", type: "quantitative", title: "TFA load", format: ".4g"}
      ]
    },
    config: {view: {stroke: null}}
  };
}

async function draw() {
  for (const view of activeViews) view.finalize();
  activeViews = [];
  const indices = matchedPairs();
  const samples = sampleValues(indices);
  const method = metricSelect.value;
  status.textContent = `${indices.length} matching taxon–function pairs · ${samples.length} samples`;
  document.getElementById("heatmap-caption").textContent =
    `Each square is one pair at the selected taxonomy level; color shows its ${{sum: "sum", mean: "mean", median: "median"}[method]} load across all samples, including zeros. Links are ranked by this metric. Taxa are summed within each sample before this summary. Click a square to inspect that pair.`;
  document.getElementById("group-caption").textContent = groupSelect.value
    ? `Each box compares sample loads for one ${groupSelect.value} category. Zeros are included.`
    : "One box summarizes all samples. Add categorical metadata to compare groups.";

  if (!indices.length || !samples.length) {
    for (const id of ["heatmap-plot", "boxplot", "samples"]) {
      document.getElementById(id).textContent = "No values to plot for this selection.";
    }
    return;
  }
  if (typeof vegaEmbed !== "function") {
    status.textContent += " · Vega could not load; an internet connection is needed to display the plots.";
    return;
  }
  try {
    const heatmapContainer = document.getElementById("heatmap");
    heatmapContainer.scrollTop = 0;
    heatmapContainer.scrollLeft = 0;
    // Use fresh targets so a previous Vega view cannot hide the next render.
    for (const id of ["heatmap-plot", "boxplot", "samples"]) {
      document.getElementById(id).replaceChildren();
    }
    const options = {actions: false, renderer: "svg"};
    const [heatmap, boxplot, samplePlot] = await Promise.all([
      vegaEmbed("#heatmap-plot", heatmapSpec(indices), options),
      vegaEmbed("#boxplot", boxplotSpec(samples), options),
      vegaEmbed("#samples", sampleSpec(samples), options)
    ]);
    activeViews = [heatmap.view, boxplot.view, samplePlot.view];
    heatmap.view.addEventListener("click", (_, item) => {
      const pair = item && item.datum;
      if (pair && pair.taxon && pair.function) {
        taxonSearch.value = "";
        functionSelect.value = pair.function;
        refreshTaxa();
        taxonSelect.value = pair.taxon;
        update();
      }
    });
  } catch (error) {
    status.textContent = `Could not draw the plots: ${error.message}`;
  }
}

function update() {
  showSelectedTaxonomy();
  renderQueue = renderQueue.then(draw).catch(error => {
    status.textContent = `Could not draw the plots: ${error.message}`;
  });
}

showUnclassified.addEventListener("change", () => {
  refreshTaxa();
  update();
});
levelSelect.addEventListener("change", () => {
  changeLevel();
  update();
});
taxonSearch.addEventListener("input", () => {
  refreshTaxa();
  update();
});
functionSelect.addEventListener("change", () => {
  refreshTaxa();
  update();
});
for (const select of [taxonSelect, groupSelect, metricSelect, limitSelect]) {
  select.addEventListener("change", update);
}
let resizeTimer;
window.addEventListener("resize", () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(update, 150);
});
update();
