"use strict";
const assert = require("node:assert/strict");
const {test} = require("node:test");
const {collapseTaxonomy, latestTaxonLabel, isUnclassifiedTaxon} = require("../../assets/tfa_explore/taxonomy.js");
const data = require("./data/taxonomy-collapse.json");

test("original taxa preserve the input data", () => {
  assert.equal(collapseTaxonomy(data, 0), data);
});

test("collapse sums within samples before computing mean and median", () => {
  const result = collapseTaxonomy(data, 3);
  const index = result.pairs.findIndex(pair => pair.taxon === "d__Bacteria; p__A; g__Shared" && pair.function === "geneA");
  const pair = result.pairs[index];
  assert.equal(pair.total, 8);
  assert.equal(pair.mean, 8 / 3);
  assert.equal(pair.median, 3);
  assert.deepEqual(result.loads.filter(load => load[0] === index), [[index, 0, 3], [index, 1, 5]]);
  assert.equal(result.pairs.filter(pair => pair.function === "geneB").length, 1);
  assert.equal(result.groups, data.groups);
});

test("identical terminal names under different parents remain distinct", () => {
  const pairs = collapseTaxonomy(data, 3).pairs.filter(pair => pair.taxon.endsWith("g__Shared") && pair.function === "geneA");
  assert.equal(pairs.length, 2);
  assert.notEqual(pairs[0].taxon_short, pairs[1].taxon_short);
});

test("missing ranks and empty rank prefixes stay under their known parent", () => {
  const pairs = collapseTaxonomy(data, 3).pairs.filter(pair => pair.taxon.includes("Unclassified"));
  assert.deepEqual(pairs.map(pair => pair.taxon), [
    "d__Bacteria; p__A; Unclassified (level 3)",
    "d__Bacteria; p__B; Unclassified (level 3)"
  ]);
});

test("total load is conserved at each depth and zeros are retained in summaries", () => {
  for (let level = 1; level <= 4; level++) {
    const result = collapseTaxonomy(data, level);
    assert.equal(result.pairs.reduce((total, pair) => total + pair.total, 0), 22);
    assert.equal(result.loads.reduce((total, load) => total + load[2], 0), 22);
  }
  const geneB = collapseTaxonomy(data, 1).pairs.find(pair => pair.function === "geneB");
  assert.equal(geneB.median, 0);
});

test("changing the level refreshes the dropdown and every plot, including heatmap clicks", async () => {
  const fs = require("node:fs");
  const path = require("node:path");
  const vm = require("node:vm");
  const elements = new Map();
  const specs = new Map();
  const clicks = new Map();
  const element = id => {
    if (!elements.has(id)) elements.set(id, {
      value: ({level: "0", metric: "sum", limit: "20"})[id] || "",
      textContent: id === "tfa-data" ? JSON.stringify(data) : "",
      checked: true, parentElement: {clientWidth: 1000}, options: [], listeners: {},
      add(option) { this.options.push(option); },
      appendChild(option) { this.options.push(option); },
      replaceChildren() { this.options = []; this.value = ""; },
      addEventListener(event, fn) { this.listeners[event] = fn; }
    });
    return elements.get(id);
  };
  const context = vm.createContext({
    document: {getElementById: element, createElement: () => ({})},
    window: {addEventListener() {}},
    Option: function(text, value) { this.textContent = text; this.value = value; },
    collapseTaxonomy, latestTaxonLabel, isUnclassifiedTaxon,
    vegaEmbed: async (id, spec) => {
      specs.set(id, spec);
      return {view: {finalize() {}, addEventListener(event, fn) { clicks.set(id, fn); }}};
    }
  });
  vm.runInContext(fs.readFileSync(path.join(__dirname, "../../assets/tfa_explore/explore.js"), "utf8"), context);
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(element("level").options.map(option => [option.textContent, option.value]), [
    ["Full taxonomy", "0"], ["Level 1", "1"], ["Level 2", "2"], ["Level 3", "3"], ["Level 4", "4"]
  ]);
  assert.equal(element("level").value, "0");
  element("level").value = "3";
  element("level").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.equal(element("taxon").options.length, 5); // All taxa plus four lineages.
  element("show-unclassified").checked = false;
  element("show-unclassified").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.equal(element("taxon").options.length, 3);
  assert.ok(specs.get("#heatmap-plot").data.values.every(pair => !pair.unclassified));
  assert.deepEqual(Array.from(specs.get("#boxplot").data.values, row => row.load), [10, 5, 2]);
  element("show-unclassified").checked = true;
  element("show-unclassified").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.equal(element("taxon").options.length, 5);
  // Function filtering clears a taxon that no longer belongs in the list.
  element("taxon").value = "d__Bacteria; p__B; g__Shared";
  element("function").value = "geneB";
  element("function").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.equal(element("taxon").value, "");
  assert.equal(element("taxon").options.length, 2);
  assert.equal(element("taxon").options[1].value, "d__Bacteria; p__A; g__Shared");
  assert.equal(element("taxon").options[1].textContent, "g__Shared");
  // Search matches the full lineage, irrespective of case.
  element("function").value = "";
  element("function").listeners.change();
  element("taxon-search").value = "P__b";
  element("taxon-search").listeners.input();
  await vm.runInContext("renderQueue", context);
  assert.equal(element("taxon").options.length, 3);
  assert.ok(element("taxon").options.slice(1).every(option => option.value.includes("p__B")));
  element("taxon-search").value = "no matching taxon";
  element("taxon-search").listeners.input();
  await vm.runInContext("renderQueue", context);
  assert.equal(element("taxon").options.length, 1);
  const shared = specs.get("#heatmap-plot").data.values.find(pair => pair.taxon === "d__Bacteria; p__A; g__Shared" && pair.function === "geneA");
  clicks.get("#heatmap-plot")({}, {datum: shared});
  await vm.runInContext("renderQueue", context);
  assert.equal(element("taxon").value, shared.taxon);
  assert.equal(element("taxon-search").value, "");
  assert.equal(element("selected-taxonomy").textContent, `Selected lineage: ${shared.taxon}`);
  assert.equal(specs.get("#heatmap-plot").encoding.y.field, "taxon");
  assert.equal(specs.get("#heatmap-plot").encoding.y.axis.labelLimit, 0);
  assert.deepEqual(Array.from(specs.get("#boxplot").data.values, row => row.load), [3, 5, 0]);
  assert.deepEqual(Array.from(specs.get("#samples").data.values, row => row.load), [5, 3, 0]);
  element("metric").value = "median";
  element("metric").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.equal(specs.get("#heatmap-plot").encoding.color.field, "median");
  assert.equal(specs.get("#heatmap-plot").data.values[0].median, 3);
  element("function").value = "";
  element("function").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.equal(specs.get("#heatmap-plot").data.values.length, 2);
  element("level").value = "1";
  element("level").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.equal(element("taxon").value, "");
  assert.equal(element("taxon").options.length, 2);
  element("limit").value = "Infinity";
  element("limit").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.equal(specs.get("#heatmap-plot").data.values.length, 2);
  assert.deepEqual(Array.from(specs.get("#boxplot").data.values, row => row.load), [11, 9, 2]);
  element("group").value = "treatment";
  element("group").listeners.change();
  await vm.runInContext("renderQueue", context);
  const boxColor = specs.get("#boxplot").encoding.color;
  const barColor = specs.get("#samples").encoding.color;
  assert.deepEqual(boxColor, barColor);
  assert.equal(boxColor.legend, null);
  assert.deepEqual(Array.from(boxColor.scale.domain), ["A", "B"]);
  // A low-load category omitted from the top 100 must not shift other colors.
  const limited = JSON.parse(JSON.stringify(vm.runInContext(`(() => {
    const samples = Array.from({length:101}, (_, index) => ({
      sample: String(index), group: index === 100 ? "Absent from bars" : "Visible",
      load: 101 - index
    }));
    return {box: boxplotSpec(samples), bars: sampleSpec(samples)};
  })()`, context)));
  assert.equal(limited.bars.data.values.length, 100);
  assert.ok(limited.bars.data.values.every(row => row.group === "Visible"));
  assert.deepEqual(limited.box.encoding.color, limited.bars.encoding.color);
  assert.deepEqual(limited.bars.encoding.color.scale.domain, ["Absent from bars", "Visible"]);
  assert.equal(limited.bars.encoding.color.legend, null);
});


test("deep levels use the last classified rank as the heatmap label", () => {
  const result = collapseTaxonomy(data, 10);
  assert.equal(result.pairs.find(pair => pair.taxon_id === "T4").taxon_label, "p__A");
  assert.equal(result.pairs.find(pair => pair.taxon_id === "T5").taxon_label, "p__B");
  assert.equal(latestTaxonLabel(["d__Bacteria", "g__", "s__unclassified"], "ID"), "d__Bacteria");
  assert.equal(latestTaxonLabel(["", "g__"], "ID"), "ID");
});


test("unclassified includes ID-only taxa and taxa missing the selected rank", () => {
  assert.equal(isUnclassifiedTaxon({lineage:["T013"], taxon_id:"T013"}, 0), true);
  assert.equal(isUnclassifiedTaxon(data.pairs[4], 0), false);
  assert.equal(isUnclassifiedTaxon(data.pairs[4], 3), true);
  assert.equal(isUnclassifiedTaxon(data.pairs[5], 3), true);
  assert.equal(isUnclassifiedTaxon(data.pairs[0], 3), false);
  assert.equal(isUnclassifiedTaxon({lineage:["d__Bacteria", "g__unclassified"]}, 0), true);
});
