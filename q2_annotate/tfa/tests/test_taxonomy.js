"use strict";
const assert = require("node:assert/strict");
const {test} = require("node:test");
const {collapseTaxonomy, latestTaxonLabel, isUnclassifiedTaxon} = require("../../assets/tfa_explore/taxonomy.js");
const data = require("./data/taxonomy-collapse.json");
// Normalize structures returned from the isolated browser context.
const plain = value => JSON.parse(JSON.stringify(value));
const expectedCollapse = require("./data/expected-taxonomy-collapse.json");

test("original taxa preserve the input data", () => {
  assert.equal(collapseTaxonomy(data, 0), data);
});

test("collapse sums within samples and preserves parents, functions, and missing ranks", () => {
  assert.deepEqual(collapseTaxonomy(data, 3), expectedCollapse["3"]);
});

test("total load is conserved at each depth and zeros are retained in summaries", () => {
  for (let level = 1; level <= 4; level++) {
    const result = collapseTaxonomy(data, level);
    assert.equal(result.pairs.reduce((total, pair) => total + pair.total, 0), 22);
    assert.equal(result.loads.reduce((total, load) => total + load[2], 0), 22);
  }
  assert.deepEqual(collapseTaxonomy(data, 1), expectedCollapse["1"]);
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
      textContent: "",
      checked: true, parentElement: {clientWidth: 1000}, options: [], listeners: {},
      add(option) { this.options.push(option); },
      appendChild(option) { this.options.push(option); },
      replaceChildren() { this.options = []; this.value = ""; },
      addEventListener(event, fn) { this.listeners[event] = fn; }
    });
    return elements.get(id);
  };
  const context = vm.createContext({
    fetch: async url => {
      assert.equal(url, "data.json");
      return {ok:true, json:async () => data};
    },
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
  await vm.runInContext("initialization", context);
  assert.deepEqual(element("level").options.map(option => [option.textContent, option.value]), [
    ["Full taxonomy", "0"], ["Level 1", "1"], ["Level 2", "2"], ["Level 3", "3"], ["Level 4", "4"]
  ]);
  assert.equal(element("level").value, "0");
  element("level").value = "3";
  element("level").listeners.change();
  await vm.runInContext("renderQueue", context);
  const expectedTaxa = ["", ...[
    "d__Bacteria; p__A; Unclassified (level 3)", "d__Bacteria; p__A; g__Shared",
    "d__Bacteria; p__B; Unclassified (level 3)", "d__Bacteria; p__B; g__Shared"
  ]];
  assert.deepEqual(element("taxon").options.map(option => option.value), expectedTaxa);
  element("show-unclassified").checked = false;
  element("show-unclassified").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(element("taxon").options.map(option => option.value), [
    "", "d__Bacteria; p__A; g__Shared", "d__Bacteria; p__B; g__Shared"
  ]);
  assert.deepEqual(plain(specs.get("#heatmap-plot").data.values), [
    expectedCollapse["3"].pairs[0], expectedCollapse["3"].pairs[2], expectedCollapse["3"].pairs[1]
  ]);
  assert.deepEqual(JSON.parse(JSON.stringify(specs.get("#boxplot").data.values)), [
    {sample:"S1", group:"All samples", load:10},
    {sample:"S2", group:"All samples", load:5},
    {sample:"S3", group:"All samples", load:2}
  ]);
  element("show-unclassified").checked = true;
  element("show-unclassified").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(element("taxon").options.map(option => option.value), expectedTaxa);
  // Function filtering clears a taxon that no longer belongs in the list.
  element("taxon").value = "d__Bacteria; p__B; g__Shared";
  element("function").value = "geneB";
  element("function").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual({
    selected:element("taxon").value,
    options:element("taxon").options.map(option => ({value:option.value, label:option.textContent}))
  }, {
    selected:"", options:[
      {value:"", label:"All taxa"},
      {value:"d__Bacteria; p__A; g__Shared", label:"g__Shared"}
    ]
  });
  // Search matches the full lineage, irrespective of case.
  element("function").value = "";
  element("function").listeners.change();
  element("taxon-search").value = "P__b";
  element("taxon-search").listeners.input();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(element("taxon").options.map(option => option.value), [
    "", "d__Bacteria; p__B; Unclassified (level 3)", "d__Bacteria; p__B; g__Shared"
  ]);
  element("taxon-search").value = "no matching taxon";
  element("taxon-search").listeners.input();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(element("taxon").options.map(option => option.value), [""]);
  const shared = specs.get("#heatmap-plot").data.values.find(pair => pair.taxon === "d__Bacteria; p__A; g__Shared" && pair.function === "geneA");
  clicks.get("#heatmap-plot")({}, {datum: shared});
  await vm.runInContext("renderQueue", context);
  assert.deepEqual({
    taxon:element("taxon").value, search:element("taxon-search").value,
    lineage:element("selected-taxonomy").textContent
  }, {taxon:shared.taxon, search:"", lineage:`Selected lineage: ${shared.taxon}`});
  assert.equal(specs.get("#heatmap-plot").encoding.y.field, "taxon");
  assert.equal(specs.get("#heatmap-plot").encoding.y.axis.labelLimit, 0);
  const expectedSamples = [
    {sample:"S1", group:"All samples", load:3},
    {sample:"S2", group:"All samples", load:5},
    {sample:"S3", group:"All samples", load:0}
  ];
  assert.deepEqual(JSON.parse(JSON.stringify(specs.get("#boxplot").data.values)), expectedSamples);
  assert.deepEqual(JSON.parse(JSON.stringify(specs.get("#samples").data.values)), [
    expectedSamples[1], expectedSamples[0], expectedSamples[2]
  ]);
  element("metric").value = "median";
  element("metric").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(JSON.parse(JSON.stringify(specs.get("#heatmap-plot").encoding.color)), {
    field:"median", type:"quantitative", title:"Median load", scale:{scheme:"blues", zero:true}
  });
  assert.deepEqual(plain(specs.get("#heatmap-plot").data.values), [expectedCollapse["3"].pairs[0]]);
  element("function").value = "";
  element("function").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(plain(specs.get("#heatmap-plot").data.values), expectedCollapse["3"].pairs.slice(0,2));
  element("level").value = "1";
  element("level").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual({selected:element("taxon").value, options:element("taxon").options.map(option => option.value)}, {
    selected:"", options:["", "d__Bacteria"]
  });
  element("limit").value = "Infinity";
  element("limit").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(plain(specs.get("#heatmap-plot").data.values), expectedCollapse["1"].pairs);
  assert.deepEqual(JSON.parse(JSON.stringify(specs.get("#boxplot").data.values)), [
    {sample:"S1", group:"All samples", load:11},
    {sample:"S2", group:"All samples", load:9},
    {sample:"S3", group:"All samples", load:2}
  ]);
  element("group").value = "treatment";
  element("group").listeners.change();
  await vm.runInContext("renderQueue", context);
  const boxColor = specs.get("#boxplot").encoding.color;
  const barColor = specs.get("#samples").encoding.color;
  const expectedColor = {
    field:"group", type:"nominal", legend:null,
    scale:{domain:["A", "B"], scheme:"tableau10"}
  };
  assert.deepEqual(JSON.parse(JSON.stringify({box:boxColor, bars:barColor})), {
    box:expectedColor, bars:expectedColor
  });
  // A low-load category omitted from the top 100 must not shift other colors.
  const limited = JSON.parse(JSON.stringify(vm.runInContext(`(() => {
    const samples = Array.from({length:101}, (_, index) => ({
      sample: String(index), group: index === 100 ? "Absent from bars" : "Visible",
      load: 101 - index
    }));
    return {box: boxplotSpec(samples), bars: sampleSpec(samples)};
  })()`, context)));
  const expectedLimitedColor = {
    field:"group", type:"nominal", legend:null,
    scale:{domain:["Absent from bars", "Visible"], scheme:"tableau10"}
  };
  assert.deepEqual({box:limited.box.encoding.color, bars:limited.bars.encoding.color}, {
    box:expectedLimitedColor, bars:expectedLimitedColor
  });
  assert.deepEqual(limited.bars.data.values, Array.from({length:100}, (_, index) => ({
    sample:String(index), group:"Visible", load:101-index
  })));
});


test("deep levels use the last classified rank as the heatmap label", () => {
  const result = collapseTaxonomy(data, 10);
  assert.deepEqual(result.pairs.map(pair => ({id:pair.taxon_id, gene:pair.function, label:pair.taxon_label})), [
    {id:"T1", gene:"geneA", label:"s__one"}, {id:"T1", gene:"geneB", label:"s__one"},
    {id:"T2", gene:"geneA", label:"s__two"}, {id:"T3", gene:"geneA", label:"g__Shared"},
    {id:"T4", gene:"geneA", label:"p__A"}, {id:"T5", gene:"geneA", label:"p__B"}
  ]);
  assert.deepEqual([
    latestTaxonLabel(["d__Bacteria", "g__", "s__unclassified"], "ID"),
    latestTaxonLabel(["", "g__"], "ID")
  ], ["d__Bacteria", "ID"]);
});


test("unclassified includes ID-only taxa and taxa missing the selected rank", () => {
  assert.deepEqual([
    isUnclassifiedTaxon({lineage:["T013"], taxon_id:"T013"}, 0),
    isUnclassifiedTaxon(data.pairs[4], 0), isUnclassifiedTaxon(data.pairs[4], 3),
    isUnclassifiedTaxon(data.pairs[5], 3), isUnclassifiedTaxon(data.pairs[0], 3),
    isUnclassifiedTaxon({lineage:["d__Bacteria", "g__unclassified"]}, 0)
  ], [true, false, true, true, false, true]);
});

test("JSON loading failures show a message and keep controls disabled", async () => {
  const fs = require("node:fs");
  const path = require("node:path");
  const vm = require("node:vm");
  const source = fs.readFileSync(path.join(__dirname, "../../assets/tfa_explore/explore.js"), "utf8");
  for (const [fetch, message] of [
    [async () => ({ok:false, status:404}), "HTTP 404"],
    [async () => { throw new Error("Network unavailable"); }, "Network unavailable"],
    [async () => ({ok:true, json:async () => { throw new Error("Invalid JSON"); }}), "Invalid JSON"]
  ]) {
    const elements = new Map();
    const getElementById = id => {
      if (!elements.has(id)) elements.set(id, {value:"", addEventListener() {}});
      return elements.get(id);
    };
    const context = vm.createContext({fetch, document:{getElementById}, window:{addEventListener() {}}});
    vm.runInContext(source, context);
    await vm.runInContext("initialization", context);
    assert.equal(getElementById("status").textContent, `Could not load TFA data: ${message}`);
    assert.equal(getElementById("taxon").disabled, true);
  }
});
