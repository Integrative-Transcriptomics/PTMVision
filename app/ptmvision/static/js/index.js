var __url = null;
var __overviewTable = null;
var __overviewChart = null;
var __dashboardChart = null;
var __dashboardContent = null;

class OverviewTable {
  tabulator = null;
  availableIdentifiers = [];
  availableModifications = [];

  constructor(id) {
    this.tabulator = new Tabulator("#" + id, {
      selectableRows: 1,
      columns: [
        {
          title: "ID",
          field: "id",
          sorter: "string",
          width: "5%",
        },
        {
          title: "Name",
          field: "name",
          sorter: "string",
          width: "25%",
        },
        {
          title: "Primary Sequence Length",
          field: "length",
          sorter: "number",
          width: "20%",
        },
        {
          title: "No. modified positions",
          field: "modified_positions",
          sorter: "number",
          width: "25%",
        },
        {
          title: "No. distinct modifications",
          field: "unique_modifications",
          sorter: "number",
          width: "25%",
        },
        {
          title: "Modifications",
          field: "modifications",
        },
      ],
    });
  }

  registerSelectionAction(action) {
    this.tabulator.on(
      "rowSelectionChanged",
      function (data, rows, selected, deselected) {
        action(data);
      }
    );
  }

  setFilters(idValues, modValues) {
    this.tabulator.clearFilter(true);
    const filters = [];
    for (let _ of modValues) {
      filters.push({
        field: "modifications",
        type: "like",
        value: _,
      });
    }
    for (let _ of idValues) {
      let idValueFields = _.split("$");
      filters.push({
        field: idValueFields[0] == "id" ? "id" : "name",
        type: "like",
        value: idValueFields[1],
      });
    }
    this.tabulator.setFilter(filters);
  }

  setData(data) {
    this.tabulator.setData(data);
    this.availableModifications = new Set();
    this.availableIdentifiers = new Set();
    for (let entry of data) {
      entry.modifications
        .split("$")
        .forEach((m) => this.availableModifications.add(m));
      this.availableIdentifiers.add(entry.id + "$" + entry.name);
    }
    this.availableModifications = [...this.availableModifications];
    this.availableIdentifiers = [...this.availableIdentifiers];
  }

  getSelection() {
    let _ = this.tabulator.getSelectedData();
    if (_.length > 0) return this.tabulator.getSelectedData()[0].id;
    else return null;
  }
}

class Chart {
  instance = null;
  #instanceDomId = null;
  #resizeObserver = null;

  constructor(id) {
    this.#instanceDomId = id;
    this.instance = echarts.init($("#" + this.#instanceDomId)[0], {
      devicePixelRatio: 2,
      renderer: "canvas",
      width: "auto",
      height: "auto",
    });
    this.#resizeObserver = new ResizeObserver((entries) => {
      this.instance.resize({
        width: entries[0].width,
        height: entries[0].height,
      });
    });
    this.#resizeObserver.observe($("#" + this.#instanceDomId)[0]);
  }

  setOption(option, replace) {
    this.instance.setOption(option, replace);
  }
}

class OverviewChart {
  chart = null;
  #data;
  #option;
  #axisStyle = {
    nameLocation: "center",
    nameTextStyle: {
      fontWeight: "bold",
      fontSize: 11,
    },
    axisTick: {
      alignWithLabel: true,
      interval: 0,
    },
  };
  #axisLabelStyle = {
    fontWeight: "lighter",
    fontSize: 9,
  };
  #titleStyle = {
    textStyle: {
      fontSize: 12,
      fontWeight: "bold",
    },
  };
  #tooltipStyle = {
    backgroundColor: "#fbfbfbe6",
    borderColor: "#fbfbfb",
    textStyle: {
      color: "#111111",
      fontSize: 13,
    },
  };
  #axisPointerLabelStyle = {
    show: true,
    fontWeight: "bold",
    fontSize: 11,
    color: "#333333",
    padding: [2, 4, 2, 4],
    backgroundColor: "#fbfbfbe6",
    borderColor: "#fbfbfb",
    margin: 1,
  };
  #sortingIndex = 1;

  constructor(id) {
    this.chart = new Chart(id);
  }

  restoreZoom() {
    if (this.#option == undefined) return;
    [0, 1].forEach((_) =>
      this.chart.instance.dispatchAction({
        type: "dataZoom",
        dataZoomIndex: _,
        start: 0,
        end: 100,
      })
    );
  }

  highlight(indices) {
    if (this.#option == undefined) return;
    let markLineOption = {
      silent: true,
      symbol: [null, null],
      lineStyle: {
        color: "#333333",
      },
      label: {
        formatter: (params) => {
          let name =
            this.#data.modificationNamesSorted[this.#sortingIndex][
              params.data.value
            ];
          return this.#data.modifications[name].display_name;
        },
        fontWeight: "lighter",
        fontSize: 9,
        position: "insideEndBottom",
      },
    };
    /* Deprecated; Error level is very low and not well visualized.
    let markAreaOption = {
      silent: true,
      itemStyle: {
        color: "#ff6663",
        opacity: 1.0,
      },
    };
    */
    this.#option.series[0].markLine = {
      ...markLineOption,
      data: indices
        .map((i) => {
          return { yAxis: i };
        })
        .concat(
          indices.map((i) => {
            return { xAxis: i };
          })
        ),
    };
    this.#option.series[1].markLine = {
      ...markLineOption,
      data: indices.map((i) => {
        return { yAxis: i };
      }),
    };

    /* Deprecated; Error level is very low and not well visualized.
    this.#option.series[1].markArea = {
      ...markAreaOption,
      data: indices.map((i) => {
        let name = this.#data.modificationNamesSorted[this.#sortingIndex][i];
        let massShift = this.#data.modifications[name].mass_shift;
        return [
          { coord: [massShift - 0.02, 0] },
          { coord: [massShift + 0.02, this.#option.yAxis[1].data.length - 1] },
        ];
      }),
    };
    */
    this.#option.series[2].markLine = {
      ...markLineOption,
      data: indices.map((i) => {
        return { yAxis: i };
      }),
    };
    this.#updateOption(false);
  }

  resort() {
    if (this.#option == undefined) return;
    if (this.#sortingIndex == 1) {
      this.#sortingIndex = 0;
      this.fill();
    } else if (this.#sortingIndex == 0) {
      this.#sortingIndex = 1;
      this.fill();
    }
  }

  getDataNames() {
    if (this.#data != undefined)
      return this.#data.modificationNamesSorted[this.#sortingIndex];
    else return [];
  }

  getDataUrl() {
    if (this.#option == undefined) return;
    return this.chart.instance.getDataURL({
      pixelRatio: 8,
      backgroundColor: "#fff",
    });
  }

  fill(data) {
    if (data != undefined)
      this.#data = {
        modifications: data[0],
        modificationNamesSorted: data[1],
        coOccurrence: data[2],
        classCounts: data[3],
      };
    let r = this.chart.instance.getWidth() / this.chart.instance.getHeight();
    // Generate series data from data.
    let dataAxis = this.#data.modificationNamesSorted[this.#sortingIndex];
    let dataMassShift = dataAxis.map((name) => {
      return this.#data.modifications[name]["mass_shift"];
    });
    let dataCount = dataAxis.map((name) => {
      return this.#data.modifications[name]["count"];
    });
    let dataCoOccurrence = [];
    for (let i = 0; i < dataAxis.length; i++) {
      for (let j = 0; j < dataAxis.length; j++) {
        if (i == j) continue;
        let _ = [dataAxis[i], dataAxis[j]].sort().join("@");
        dataCoOccurrence.push([i, j, this.#data.coOccurrence[_]]);
      }
    }
    // Fill in option.
    this.#option = {
      title: [
        {
          text:
            "Total Modifications (e.g. Phospho on Ser249) " +
            Object.values(this.#data.classCounts).reduce((a, b) => a + b, 0) +
            " | Modification Types (e.g. Phospho) " +
            Object.keys(this.#data.modifications).length +
            " | Modification Classes " +
            Object.keys(this.#data.classCounts).length,
          top: "top",
          left: "center",
          ...this.#titleStyle,
        },
        {
          text: "Shared PTM sites between modification types",
          top: 15,
          left: "10%",
          ...this.#titleStyle,
        },
        {
          text: "Mass Shift",
          top: 30,
          left: 10 + 2 + 80 / r + "%",
          ...this.#titleStyle,
        },
        {
          text: "Site Count",
          top: 30,
          left: 10 + 2 * 2 + 12 + 80 / r + "%",
          ...this.#titleStyle,
        },
        {
          text: "PTM Class Counts",
          top: 30,
          left: 10 + 4 * 2 + 2 * 12 + 80 / r + "%",
          ...this.#titleStyle,
        },
      ],
      grid: [
        {
          top: 50,
          left: "10%",
          height: "80%",
          width: 80 / r + "%",
          show: true,
        },
        {
          top: 50,
          left: 10 + 2 + 80 / r + "%",
          height: "80%",
          width: "12%",
          show: true,
        },
        {
          top: 50,
          left: 10 + 2 * 2 + 12 + 80 / r + "%",
          height: "80%",
          width: "12%",
          show: true,
        },
        {
          top: 50,
          bottom: "24%",
          left: 10 + 4 * 2 + 2 * 12 + 80 / r + "%",
          right: "4%",
          height: "auto",
          width: "auto",
          show: true,
        },
      ],
      xAxis: [
        {
          gridIndex: 0,
          type: "category",
          name: "Modification",
          nameLocation: "center",
          ...this.#axisStyle,
          nameGap: 35,
          data: dataAxis,
          axisTick: {
            alignWithLabel: true,
            interval: 0,
            length: 3,
          },
          axisLabel: {
            show: true,
            formatter: (i) => {
              if (this.#data.modifications[i]["display_name"].length > 4) {
                return (
                  this.#data.modifications[i]["display_name"].substring(0, 5) +
                  "..."
                );
              } else {
                return this.#data.modifications[i]["display_name"];
              }
            },
            rotate: 45,
            ...this.#axisLabelStyle,
          },
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
        },
        {
          gridIndex: 1,
          type: "value",
          name: "Mass Shift [Da]",
          ...this.#axisStyle,
          nameGap: 35,
          axisLabel: {
            show: true,
            interval: 0,
            rotate: 45,
            ...this.#axisLabelStyle,
          },
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
        },
        {
          gridIndex: 2,
          type: "value",
          name: "Count",
          ...this.#axisStyle,
          nameGap: 35,
          axisLabel: {
            show: true,
            interval: 0,
            rotate: 45,
            ...this.#axisLabelStyle,
          },
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
        },
        {
          gridIndex: 3,
          type: "category",
          name: "Modification Class",
          ...this.#axisStyle,
          nameGap: 95,
          data: Object.keys(this.#data.classCounts),
          axisTick: {
            alignWithLabel: true,
            interval: 0,
            length: 4,
          },
          axisLabel: {
            show: true,
            interval: 0,
            rotate: 45,
            ...this.#axisLabelStyle,
          },
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
          },
        },
      ],
      yAxis: [
        {
          gridIndex: 0,
          type: "category",
          name: "Modification",
          ...this.#axisStyle,
          nameGap: 140,
          data: dataAxis,
          inverse: true,
          axisTick: {
            alignWithLabel: true,
            interval: 0,
            length: 2,
          },
          axisLabel: {
            show: true,
            formatter: (i) => {
              if (this.#data.modifications[i]["display_name"].length > 20) {
                return (
                  this.#data.modifications[i]["display_name"].substring(0, 20) +
                  "..."
                );
              } else {
                return this.#data.modifications[i]["display_name"];
              }
            },
            ...this.#axisLabelStyle,
          },
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
        },
        {
          gridIndex: 1,
          type: "category",
          data: dataAxis,
          show: false,
          inverse: true,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
          },
        },
        {
          gridIndex: 2,
          type: "category",
          data: dataAxis,
          show: false,
          inverse: true,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
          },
        },
        {
          gridIndex: 3,
          type: "value",
          name: "Count",
          ...this.#axisStyle,
          nameGap: 40,
          axisLabel: {
            show: true,
            ...this.#axisLabelStyle,
          },
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
        },
      ],
      tooltip: {
        ...this.#tooltipStyle,
        formatter: (params) => {
          if (Array.isArray(params)) params = params[0];
          if (params.seriesIndex == 0) {
            let yModName =
              this.#data.modificationNamesSorted[this.#sortingIndex][
                params.data[1]
              ];
            let xModName =
              this.#data.modificationNamesSorted[this.#sortingIndex][
                params.data[0]
              ];
            return (
              "<code>" +
              params.data[2] +
              "</code> sites have been modified by both <code>" +
              this.#data.modifications[yModName].display_name +
              "</code> (" +
              parseFloat(
                String(this.#data.modifications[yModName].mass_shift)
              ).toFixed(2) +
              " Da) and <code>" +
              this.#data.modifications[xModName].display_name +
              "</code> (" +
              parseFloat(
                String(this.#data.modifications[xModName].mass_shift)
              ).toFixed(2) +
              " Da)."
            );
          }
          if (params.seriesIndex == 1) {
            return (
              "Modification <code>" +
              this.#data.modifications[params.name].display_name +
              "</code> assigned mass shift is <code>" +
              params.data +
              "</code> Da."
            );
          }
          if (params.seriesIndex == 2)
            return (
              "Modification <code>" +
              this.#data.modifications[params.name].display_name +
              "</code> observed on <code>" +
              params.data +
              "</code> different sites, spread across <code>" +
              this.#data.modifications[params.name].frequency +
              "%</code> of proteins in data."
            );
          if (params.seriesIndex == 3)
            return (
              "Modification class <code>" +
              params.name +
              "</code> assigned <code>" +
              params.data +
              "</code> times."
            );
        },
      },
      visualMap: [
        {
          type: "continuous",
          seriesIndex: [0],
          inRange: {
            color: [
              "#dddddd",
              "#cccccc",
              "#888888",
              "#666666",
              "#444444",
              "#222222",
              "#000000",
            ],
          },
          outOfRange: {
            color: ["#444444"],
          },
          min: 1,
          max: Math.max(...Object.values(this.#data.coOccurrence)),
          orient: "horizontal",
          top: 30,
          left: "10%",
          itemHeight: 50,
          itemWidth: 11,
          text: [
            Math.max(...Object.values(this.#data.coOccurrence)),
            "No. shared sites 1",
          ],
          textStyle: { fontWeight: "lighter", fontSize: 11 },
        },
      ],
      dataZoom: [
        {
          type: "inside",
          yAxisIndex: [0, 1, 2],
          brushSelect: false,
          throttle: 0,
        },
        {
          type: "inside",
          xAxisIndex: [0],
          throttle: 0,
        },
      ],
      series: [
        {
          type: "heatmap",
          xAxisIndex: 0,
          yAxisIndex: 0,
          progressive: 1000,
          progressiveThreshold: 1500,
          animation: false,
          itemStyle: {
            borderWidth: 0.1,
            borderRadius: 2,
            borderColor: "#fbfbfb",
          },
          data: dataCoOccurrence,
          emphasis: {
            itemStyle: {
              color: "inherit",
              borderColor: "#62a8ac",
              borderWidth: 2,
            },
          },
          markLine: {},
        },
        {
          type: "scatter",
          xAxisIndex: 1,
          yAxisIndex: 1,
          animation: false,
          symbolSize: 8,
          symbol: "diamond",
          itemStyle: {
            color: "#111111",
          },
          data: dataMassShift,
          cursor: "default",
          markLine: {},
          markArea: {},
        },
        {
          type: "bar",
          xAxisIndex: 2,
          yAxisIndex: 2,
          animation: false,
          itemStyle: {
            color: "#111111",
          },
          barWidth: "50%",
          data: dataCount,
          cursor: "default",
          markLine: {},
        },
        {
          type: "bar",
          xAxisIndex: 3,
          yAxisIndex: 3,
          animation: false,
          itemStyle: {
            color: "#111111",
          },
          barWidth: "50%",
          data: Object.values(this.#data.classCounts),
          cursor: "default",
        },
      ],
    };
    // Highlight PTMs that fall within 0.02 Da mass shift tolerance, when sorting by mass shift.
    if (this.#sortingIndex == 0) {
      let segments = [];
      let segment = new Set();
      let areas = [];
      let markAreaOption = {
        silent: true,
        animation: false,
        itemStyle: {
          color: "#ff6663",
          opacity: 0.42,
        },
      };
      for (let i = dataMassShift.length - 1; i >= 0; i--) {
        let j = i - 1;
        if (j < 0) break;
        if (dataMassShift[i] == "null" || dataMassShift[j] == "null") continue;
        if (dataMassShift[i] + 0.02 >= dataMassShift[j]) {
          segment.add(i);
          segment.add(j);
        } else {
          if (segment.size > 0) {
            segments.push(segment);
            segment = new Set();
          }
        }
      }
      if (segments.length > 0) {
        for (let S of segments) {
          let s = Array.from(S);
          let l = s.slice(0)[0];
          let u = s.slice(-1)[0];
          areas.push([{ coord: ["min", u] }, { coord: ["max", l] }]);
        }
        this.#option.series[1].markArea = {
          ...markAreaOption,
          data: areas,
        };
      }
    }
    this.#updateOption(true);
  }

  #updateOption(replace) {
    if (this.#option != undefined) this.chart.setOption(this.#option, replace);
  }
}

class DashboardChart {
  chart = null;
  structure = null;
  #data;
  #option;
  #colors = {
    M: {},
    i: 0,
    C: [
      "#f44336",
      "#00bcd4",
      "#009688",
      "#4caf50",
      "#cddc39",
      "#e81e63",
      "#9c27b0",
      "#3f51b5",
      "#ffc107",
      "#ff9800",
      "#333333",
    ],
    explicit: {
      Helix: "#ee858d",
      Turn: "#edc585",
      "Beta strand": "#85b2ed",
      "Alternative sequence": "#210203",
      "Topological domain": "#525354",
      Transmembrane: "#697268",
      Intramembrane: "#697268",
      Domain: "#525354",
      Repeat: "#090C08",
      "Zinc finger": "#95A3A4",
      "DNA binding": "#95A3A4",
      Region: "#CC5A71",
      "Coiled coil": "#95A3A4",
      Motif: "#CC5A71",
      "Compositional bias": "#34344A",
      contact_unmodified: "#AAAAAA",
      contact_modified: "#000000",
      contact_highlight: "#DC5754",
    },
  };
  #axisStyle = {
    nameLocation: "center",
    nameTextStyle: {
      fontWeight: "bold",
      fontSize: 11,
    },
    axisTick: {
      alignWithLabel: true,
      interval: 0,
    },
  };
  #axisLabelStyle = {
    fontWeight: "lighter",
    fontSize: 9,
  };
  #titleStyle = {
    textStyle: {
      fontSize: 11,
      fontWeight: "bold",
    },
  };
  #tooltipStyle = {
    backgroundColor: "#fbfbfbe6",
    borderColor: "#fbfbfb",
    textStyle: {
      color: "#111111",
      fontSize: 13,
    },
  };
  #axisPointerLabelStyle = {
    show: true,
    fontWeight: "bold",
    fontSize: 11,
    color: "#333333",
    padding: [2, 4, 2, 4],
    backgroundColor: "#fbfbfbe6",
    borderColor: "#fbfbfb",
    margin: 1,
  };
  #aminoAcids = [
    "ALA",
    "CYS",
    "ASP",
    "GLU",
    "PHE",
    "GLY",
    "HIS",
    "ILE",
    "LYS",
    "LEU",
    "MET",
    "ASN",
    "PRO",
    "GLN",
    "ARG",
    "SER",
    "THR",
    "VAL",
    "TRP",
    "TYR",
  ];
  #contentMode = 1;

  constructor(id1, id2) {
    this.chart = new Chart(id1);
    this.structure = new StructureView(id2);
    this.hideStructure();
  }

  showStructure() {
    $("#" + this.structure.domId).show();
  }

  hideStructure() {
    $("#" + this.structure.domId).hide();
  }

  restoreZoom() {
    if (this.#option == undefined) return;
    let I;
    if (this.#contentMode == 1) {
      I = [0, 1, 2];
    } else {
      I = [0, 1, 2];
      this.structure.glviewer.zoomTo();
    }
    I.forEach((_) =>
      this.chart.instance.dispatchAction({
        type: "dataZoom",
        dataZoomIndex: _,
        start: 0,
        end: 100,
      })
    );
  }

  highlight(indices) {
    if (this.#option == undefined) return;
    if (this.#contentMode == 1) {
      let _ = {
        silent: true,
        symbol: [null, null],
        lineStyle: {
          color: "#333333",
        },
        label: {
          formatter: (params) => {
            return this.#data.modifications[params.data.value];
          },
          fontWeight: "lighter",
          fontSize: 10,
          position: "insideEndBottom",
        },
        data: indices.map((i) => {
          return { yAxis: parseInt(i) };
        }),
        animation: false,
      };
      for (let i = 0; i < this.#option.series.length; i++) {
        let series = this.#option.series[i];
        if (
          series.id.startsWith("aacount") ||
          series.id.startsWith("mdcount") ||
          series.id.startsWith("modifications")
        ) {
          this.#option.series[i].markLine = _;
        }
      }
      this.#updateOption(false);
    } else {
      let modifications = indices.map((i) => this.#data.modifications.at(i));
      let positions = [];
      if (modifications.length > 0)
        for (const [position, entry] of Object.entries(this.#data.positions)) {
          if (
            this.#contains(
              modifications,
              entry.modifications.map((_) => _.display_name)
            )
          )
            positions.push(position);
        }
      this.structure.highlightPositions(positions, modifications);
      this.setContactsOption(false, modifications);
      this.#updateOption(true);
    }
  }

  switchContent() {
    if (this.#option == undefined) return;
    if (this.#contentMode == 1) {
      this.showStructure();
      this.setContactsOption(true, null);
      this.#contentMode = 2;
    } else {
      this.hideStructure();
      this.setModificationsOption();
      this.#contentMode = 1;
    }
    this.#updateOption(true);
  }

  getDataNames() {
    if (this.#data != undefined) return this.#data.modifications;
    else return [];
  }

  getDataUrl() {
    if (this.#option == undefined) return;
    return this.chart.instance.getDataURL({
      pixelRatio: 8,
      backgroundColor: "#fff",
    });
  }

  setModificationsOption() {
    // Populate per modification counts series.
    let mdcountSeries = {};
    for (const name of this.#data.modifications) {
      for (const [cls, count] of Object.entries(
        this.#data.modificationCounts[name]
      )) {
        if (!mdcountSeries.hasOwnProperty(cls))
          mdcountSeries[cls] = {
            id: "mdcount@" + cls,
            name: cls,
            stack: "total",
            type: "bar",
            data: [],
            itemStyle: {
              color: this.#getColor(cls),
            },
            xAxisIndex: 4,
            yAxisIndex: 4,
            cursor: "default",
            emphasis: {
              disabled: true,
            },
          };
        mdcountSeries[cls].data.push([
          count,
          this.#data.modifications.indexOf(name),
        ]);
      }
    }
    // Populate per position counts series.
    let pscountSeries = {};
    for (const position of Object.keys(this.#data.positions)) {
      for (const [cls, count] of Object.entries(
        this.#data.positions[position].counts
      )) {
        if (!pscountSeries.hasOwnProperty(cls))
          pscountSeries[cls] = {
            id: "pscount@" + cls,
            name: cls,
            type: "bar",
            stack: "total2",
            data: [],
            itemStyle: {
              color: this.#getColor(cls),
            },
            xAxisIndex: 0,
            yAxisIndex: 0,
            cursor: "default",
            emphasis: {
              disabled: true,
            },
          };
        pscountSeries[cls].data.push([position - 1, count]);
      }
    }
    // Populate aminoacid count series.
    let aacountSeries = {
      id: "aacount",
      name: "Aminoacid Counts",
      type: "heatmap",
      data: [],
      xAxisIndex: 3,
      yAxisIndex: 3,
      itemStyle: {
        borderWidth: 0.2,
        borderRadius: 2,
        borderColor: "#fbfbfb",
      },
      cursor: "default",
      emphasis: {
        disabled: true,
      },
    };
    for (const [aa, _] of Object.entries(this.#data.aminoacidCounts)) {
      for (const [mdname, count] of Object.entries(_)) {
        aacountSeries.data.push([
          this.#aminoAcids.indexOf(aa),
          this.#data.modifications.indexOf(mdname),
          count,
        ]);
      }
    }
    // Populate annotation series.
    let annotationSeries = {};
    let annotationLabels = [];
    for (const position of Object.keys(this.#data.positions)) {
      for (const [cls, entries] of Object.entries(
        this.#data.positions[position].annotations
      )) {
        if (!annotationLabels.includes(cls)) annotationLabels.push(cls);
        if (!annotationSeries.hasOwnProperty(cls))
          annotationSeries[cls] = {
            id: "annotation@" + cls,
            name: cls,
            type: "heatmap",
            data: [],
            xAxisIndex: 2,
            yAxisIndex: 2,
            cursor: "default",
            emphasis: {
              disabled: true,
            },
          };
        annotationSeries[cls].data.push({
          value: [position - 1, annotationLabels.indexOf(cls), 1, entries],
          itemStyle: {
            borderColor: "transparent",
            borderWidth: 0,
            color: this.#colors.explicit.hasOwnProperty(entries[0][0])
              ? this.#colors.explicit[entries[0][0]]
              : "#40434E",
          },
        });
      }
    }
    // Populate modifications map series.
    let modificationsSeries = {};
    for (const position of Object.keys(this.#data.positions)) {
      for (let modification of Object.values(
        this.#data.positions[position].modifications
      )) {
        if (
          !modificationsSeries.hasOwnProperty(
            modification.modification_classification
          )
        )
          modificationsSeries[modification.modification_classification] = {
            id: "modifications@" + modification.modification_classification,
            name: modification.modification_classification,
            type: "heatmap",
            data: [],
            itemStyle: {
              color: this.#getColor(modification.modification_classification),
              borderWidth: 0.2,
              borderRadius: 2,
              borderColor: "#fbfbfb",
            },
            xAxisIndex: 1,
            yAxisIndex: 1,
            cursor: "default",
            emphasis: {
              disabled: true,
            },
            animation: false,
            large: true,
          };
        modificationsSeries[modification.modification_classification].data.push(
          [
            position - 1,
            this.#data.modifications.indexOf(modification.display_name),
            0,
            modification.display_name,
            modification.mass_shift,
            modification.modification_classification,
            this.#data.sequence[position - 1],
          ]
        );
      }
    }
    // Construct option object.
    let n = this.#data.sequence.length;
    this.#option = {
      title: [
        /*{
          top: 4,
          left: 100,
          text: "PTM Classes",
          textStyle: {
            fontSize: 10,
          },
        },*/
        {
          text: "Per Position PTMs and Annotations",
          top: "2%",
          left: "8%",
          ...this.#titleStyle,
        },
        {
          text: "Per Aminoacid PTM Counts",
          top: "12%",
          left: "64%",
          ...this.#titleStyle,
        },
        {
          text: "Per PTM Class Counts",
          top: "14%",
          left: "85%",
          ...this.#titleStyle,
        },
      ],
      grid: [
        {
          // Position counts.
          top: "5%",
          left: "8%",
          height: "11%",
          width: "55%",
          zlevel: 0,
          show: true,
        },
        {
          // Modifications map.
          top: "17%",
          left: "8%",
          height: "65%",
          width: "55%",
          containLabel: false,
          zlevel: 0,
          show: true,
        },
        {
          // Annotations.
          top: "83%",
          left: "8%",
          height: "10%",
          width: "55%",
          containLabel: false,
          zlevel: 0,
          show: true,
        },
        {
          // Amino-acid counts.
          top: "17%",
          left: "64%",
          height: "65%",
          width: "20%",
          containLabel: false,
          zlevel: 0,
          show: true,
        },
        {
          // PTM counts.
          top: "17%",
          left: "85%",
          height: "65%",
          width: "8%",
          containLabel: false,
          zlevel: 0,
          show: true,
        },
      ],
      xAxis: [
        {
          // Position counts.
          data: Object.keys(this.#data.positions),
          show: false,
          gridIndex: 0,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: true,
          },
        },
        {
          // Modifications map.
          data: Object.keys(this.#data.positions),
          show: false,
          gridIndex: 1,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
        },
        {
          // Annotations.
          data: Object.keys(this.#data.positions),
          name: "Protein Position",
          nameGap: 30,
          nameLocation: "center",
          nameTextStyle: {
            fontWeight: "bold",
            fontSize: 11,
          },
          axisTick: {
            interval: "auto",
          },
          axisLabel: {
            ...this.#axisLabelStyle,
            interval: Math.round(n / 50),
          },
          gridIndex: 2,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: true,
          },
        },
        {
          // Amino acids.
          data: this.#aminoAcids,
          name: "Amino Acid",
          nameGap: 30,
          ...this.#axisStyle,
          axisLabel: {
            ...this.#axisLabelStyle,
          },
          gridIndex: 3,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
        },
        {
          // Modification counts.
          type: "value",
          name: "Count",
          nameGap: 30,
          ...this.#axisStyle,
          axisLabel: {
            ...this.#axisLabelStyle,
          },
          gridIndex: 4,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
        },
      ],
      yAxis: [
        {
          // Position counts.
          type: "value",
          name: "Count",
          nameGap: 30,
          ...this.#axisStyle,
          axisLabel: {
            ...this.#axisLabelStyle,
          },
          gridIndex: 0,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
        },
        {
          // Modifications map.
          type: "category",
          data: this.#data.modifications,
          name: "Modification",
          nameGap: 120,
          ...this.#axisStyle,
          axisLabel: {
            ...this.#axisLabelStyle,
            formatter: (value) => {
              return value.length > 20 ? value.substring(0, 21) + "..." : value;
            },
          },
          gridIndex: 1,
          inverse: true,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
        },
        {
          // Annotations.
          // type: "category",
          data: annotationLabels,
          name: "Annotation",
          nameGap: 120,
          ...this.#axisStyle,
          axisLabel: {
            ...this.#axisLabelStyle,
          },
          gridIndex: 2,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
        },
        {
          // Amino acids.
          type: "category",
          data: this.#data.modifications,
          show: false,
          gridIndex: 3,
          inverse: true,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
        },
        {
          // Modification counts.
          type: "category",
          data: this.#data.modifications,
          show: false,
          gridIndex: 4,
          inverse: true,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: true,
          },
        },
      ],
      dataZoom: [
        {
          type: "inside",
          xAxisIndex: [0, 1, 2],
          throttle: 0,
          filterMode: "none",
        },
        {
          type: "inside",
          yAxisIndex: [1, 3, 4],
          throttle: 0,
        },
        {
          type: "inside",
          xAxisIndex: [3],
          throttle: 0,
        },
      ],
      legend: [
        {
          id: "legend",
          top: "top",
          left: "center",
          zlevel: 2,
          icon: "circle",
          itemWidth: 10,
          itemHeight: 10,
          orient: "horizontal",
          textStyle: {
            fontSize: 10,
            fontWeight: "lighter",
          },
          data: [],
          selector: [
            { type: "all", title: "Select all." },
            { type: "inverse", title: "Invert selection." },
          ],
          selectorLabel: {
            fontSize: 10,
            fontWeight: "lighter",
            borderRadius: 2,
          },
        },
      ],
      tooltip: [
        {
          trigger: "item",
          ...this.#tooltipStyle,
        },
      ],
      series: [
        ...Object.values(pscountSeries),
        ...Object.values(mdcountSeries),
        aacountSeries,
        ...Object.values(annotationSeries),
        ...Object.values(modificationsSeries),
      ],
      visualMap: [],
    };
    // Append information from series definitions to option object.
    let legendData = new Set(
      this.#option.series
        .filter(
          (_) => !_.id.startsWith("annotation@") && !_.id.startsWith("aacount")
        )
        .map((_) => {
          return _.name;
        })
    );
    this.#option.legend[0].data = [];
    legendData.forEach((_) => {
      this.#option.legend[0].data.push({ name: _ });
    });
    let aacountSeriesMax = Math.max(...aacountSeries.data.map((_) => _[2]));
    this.#option.visualMap.push({
      seriesIndex: this.#option.series.indexOf(aacountSeries),
      top: "14%",
      left: "64%",
      color: ["#111111", "#d1d1d1"],
      orient: "horizontal",
      itemHeight: 50,
      itemWidth: 11,
      precision: 0,
      textStyle: {
        fontSize: 10,
        fontWeight: "lighter",
      },
      text: [aacountSeriesMax, "No. Occurrence 1"],
      min: 1,
      max: aacountSeriesMax,
    });
    this.#option.tooltip[0].formatter = (params) => {
      let contentHead = ``;
      let contentBody = ``;
      let component = Array.isArray(params) ? params[0] : params;
      let positionIndex;
      let modificationIndex;
      let aminoacidIndex;
      if (component.seriesId.startsWith("annotation"))
        // Interaction on position axis.
        positionIndex = component.data.value[0];
      if (
        component.seriesId.startsWith("pscount") ||
        component.seriesId.startsWith("modifications")
      )
        // Interaction on position axis.
        positionIndex = component.data[0];
      if (
        component.seriesId.startsWith("aacount") ||
        component.seriesId.startsWith("mdcount") ||
        component.seriesId.startsWith("modifications")
      )
        // Interaction on modification axis.
        modificationIndex = component.data[1];
      if (component.seriesId.startsWith("aacount"))
        // Interaction on aminoacid axis.
        aminoacidIndex = component.data[0];
      // Fill content head based on axis.
      if (positionIndex != undefined)
        contentHead +=
          "Position <code>" +
          (positionIndex + 1) +
          "&nbsp;" +
          this.#data.sequence[positionIndex] +
          "</code> ";
      if (modificationIndex != undefined)
        contentHead +=
          "Modification <code>" +
          this.#data.modifications[modificationIndex] +
          "</code> ";
      if (component.seriesId.startsWith("modifications"))
        contentHead += " as <code>" + component.seriesName + "</code>";
      if (aminoacidIndex != undefined)
        contentHead +=
          "has " +
          this.#data.aminoacidCounts[this.#aminoAcids[aminoacidIndex]][
            this.#data.modifications[modificationIndex]
          ] +
          " occurrences on aminoacid <code>" +
          this.#aminoAcids[aminoacidIndex] +
          "</code>";
      // Fill content body based on axis.
      var noData;
      if (positionIndex != undefined) {
        noData = true;
        contentBody += `<hr /><small><b>Position PTM class counts</b></small></br>`;
        for (const [cls, count] of Object.entries(
          this.#data.positions[positionIndex + 1].counts
        )) {
          noData = false;
          contentBody +=
            `<small><code>` + count + `</code>&nbsp;` + cls + `</small></br>`;
        }
        if (noData) contentBody += `<small>No data.</small>`;
        noData = true;
        contentBody += `<hr /><small><b>Position annotations</b></small></br>`;
        for (const [cls, info] of Object.entries(
          this.#data.positions[positionIndex + 1].annotations
        )) {
          for (const i of info) {
            noData = false;
            contentBody +=
              `<small><code>` +
              cls +
              `</code> at position ` +
              i[1] +
              ` to ` +
              i[2] +
              `&nbsp;` +
              i[0] +
              (i[3] != "" ? `&nbsp;` + i[3] : ``) +
              `</br>`;
          }
          contentBody += `</small>`;
        }
        if (noData) contentBody += `<small>No data.</small>`;
      }
      if (modificationIndex != undefined) {
        noData = true;
        contentBody += `<hr /><small><b>PTM class counts</b></small></br>`;
        for (const [cls, count] of Object.entries(
          this.#data.modificationCounts[
            this.#data.modifications[modificationIndex]
          ]
        )) {
          noData = false;
          contentBody +=
            `<small><code>` + count + `</code>&nbsp;` + cls + `</small></br>`;
        }
        if (noData) contentBody += `<small>No data.</small>`;
      }
      return contentHead + contentBody;
    };
    this.chart.instance.on("click", () => {});
  }

  setContactsOption(setStructure, highlightModifications) {
    // Populate per position counts series.
    let pscountSeries = {};
    for (const position of Object.keys(this.#data.positions)) {
      for (const [cls, count] of Object.entries(
        this.#data.positions[position].counts
      )) {
        if (!pscountSeries.hasOwnProperty(cls))
          pscountSeries[cls] = {
            id: "pscount@" + cls,
            name: cls,
            type: "bar",
            stack: "total2",
            data: [],
            itemStyle: {
              color: "#333333",
            },
            xAxisIndex: 0,
            yAxisIndex: 0,
            cursor: "default",
            emphasis: {
              disabled: true,
            },
          };
        pscountSeries[cls].data.push([position - 1, count]);
      }
    }
    // Populate annotation series.
    let annotationSeries = {};
    let annotationLabels = [];
    for (const position of Object.keys(this.#data.positions)) {
      for (const [cls, entries] of Object.entries(
        this.#data.positions[position].annotations
      )) {
        if (!annotationLabels.includes(cls)) annotationLabels.push(cls);
        if (!annotationSeries.hasOwnProperty(cls))
          annotationSeries[cls] = {
            id: "annotation@" + cls,
            name: cls,
            type: "heatmap",
            data: [],
            xAxisIndex: 2,
            yAxisIndex: 2,
            cursor: "default",
            emphasis: {
              disabled: true,
            },
          };
        annotationSeries[cls].data.push({
          value: [position - 1, annotationLabels.indexOf(cls), 1, entries],
          itemStyle: {
            borderColor: "transparent",
            borderWidth: 0,
            color: this.#colors.explicit.hasOwnProperty(entries[0][0])
              ? this.#colors.explicit[entries[0][0]]
              : "#40434E",
          },
        });
      }
    }
    // Populate contact map series.
    let contactsSeries = {};
    var setUnion = (s1, s2) => {
      let _ = new Set();
      s1.forEach((e) => _.add(e));
      s2.forEach((e) => _.add(e));
      return _;
    };
    for (let [index_i, contacts] of Object.entries(this.#data.contacts)) {
      let x = parseInt(index_i);
      for (let contact of contacts) {
        let y = contact[0];
        let iModifications = new Set(
          this.#data.positions[x].modifications.map((_) => _.display_name)
        );
        let jModifications = new Set(
          this.#data.positions[y].modifications.map((_) => _.display_name)
        );
        let unionSize = setUnion(iModifications, jModifications).size;
        let cls;
        let clr;
        if (unionSize == 0) {
          cls = "Unmodified Contact";
          clr = this.#colors.explicit.contact_unmodified;
        } else if (unionSize > 0) {
          cls = "Modified Contact";
          clr = this.#colors.explicit.contact_modified;
          if (
            highlightModifications != undefined &&
            highlightModifications.length > 0 &&
            (this.#contains(highlightModifications, [...iModifications]) ||
              this.#contains(highlightModifications, [...jModifications]))
          ) {
            cls =
              "Modified Contact (Includes " +
              highlightModifications.join(",") +
              ")";
            clr = this.#colors.explicit.contact_highlight;
          }
        } else {
          continue;
        }
        if (!contactsSeries.hasOwnProperty(cls))
          contactsSeries[cls] = {
            id: "contacts@" + cls,
            name: cls,
            type: "heatmap",
            data: [],
            xAxisIndex: 1,
            yAxisIndex: 1,
            large: true,
            cursor: "default",
            emphasis: {
              disabled: true,
            },
            itemStyle: {
              color: clr,
              borderColor: clr,
              borderWidth: 2,
              borderRadius: 1,
            },
          };
        contactsSeries[cls].data.push([x - 1, y - 1, 1]);
      }
    }
    // Initialize structure view.
    if (setStructure)
      this.structure.setStructure(this.#data.structure, this.#data.contacts);
    // Construct option object.
    let r = this.chart.instance.getWidth() / this.chart.instance.getHeight();
    let n = this.#data.sequence.length;
    console.log(n);
    $("#panel-dashboard-structure").css("left", 15 + 75 / r + "%");
    $("#panel-dashboard-structure").css("width", 100 - (20 + 75 / r) + "%");
    this.#option = this.#option = {
      title: [
        /*{
          top: 4,
          left: "left",
          text: "Contact Classes",
          textStyle: {
            fontSize: 10,
          },
        },*/
        {
          text: "Per Position PTM Counts, Residue Contacts and Annotations",
          top: "4%",
          left: "8%",
          ...this.#titleStyle,
        },
        {
          text: "Protein Structure and Residue Contact PTM Details",
          top: "4%",
          left: 15 + 75 / r + "%",
          ...this.#titleStyle,
        },
        {
          top: "80%",
          left: 15 + 75 / r + "%",
          text: "No PTM details to be seen. Click on a cell representing modified residues in the heatmap of residue contacts to populate this chart.",
          textStyle: {
            fontSize: 13,
            fontWeight: "lighter",
          },
          backgroundColor: "#FAFAFA",
          show: true,
        },
      ],
      grid: [
        {
          // Position counts.
          top: "7%",
          left: "8%",
          height: "5%",
          width: 75 / r + "%",
          zlevel: 0,
          show: true,
        },
        {
          // Contact map.
          top: "13%",
          left: "8%",
          height: "75%",
          width: 75 / r + "%",
          containLabel: false,
          zlevel: 0,
          show: true,
          backgroundColor: "#FFFFFF",
        },
        {
          // Annotations.
          top: "89%",
          left: "8%",
          height: "5%",
          width: 75 / r + "%",
          containLabel: false,
          zlevel: 0,
          show: true,
        },
        {
          // Contact detail.
          top: "80%",
          left: 15 + 75 / r + "%",
          width: 100 - (20 + 75 / r) + "%",
          height: "13%",
          containLabel: false,
          zlevel: 0,
          show: false,
        },
      ],
      xAxis: [
        {
          // Position counts.
          type: "category",
          data: Object.keys(this.#data.positions),
          show: false,
          gridIndex: 0,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: true,
          },
        },
        {
          // Contact map.
          type: "category",
          data: Object.keys(this.#data.positions),
          gridIndex: 1,
          show: false,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
        },
        {
          // Annotations.
          type: "category",
          data: Object.keys(this.#data.positions),
          name: "Protein Position",
          nameGap: 30,
          nameLocation: "center",
          nameTextStyle: {
            fontWeight: "bold",
            fontSize: 11,
          },
          axisTick: {
            interval: "auto",
          },
          axisLabel: {
            ...this.#axisLabelStyle,
            interval: Math.round(n / 50),
            rotate: 45,
          },
          gridIndex: 2,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: true,
          },
        },
        {
          // Contact detail.
          type: "value",
          name: "Mass Shift [Da]",
          nameGap: 30,
          ...this.#axisStyle,
          axisLabel: {
            ...this.#axisLabelStyle,
          },
          gridIndex: 3,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
        },
      ],
      yAxis: [
        {
          // Position counts.
          type: "value",
          name: "Count",
          nameGap: 30,
          nameLocation: "center",
          nameTextStyle: {
            fontWeight: "bold",
            fontSize: 11,
          },
          axisTick: {
            alignWithLabel: true,
            interval: "auto",
          },
          axisLabel: {
            ...this.#axisLabelStyle,
          },
          gridIndex: 0,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
        },
        {
          // Contact map.
          type: "category",
          data: Object.keys(this.#data.positions),
          name: "Position",
          nameGap: 30,
          nameLocation: "center",
          nameTextStyle: {
            fontWeight: "bold",
            fontSize: 11,
          },
          axisTick: {
            interval: "auto",
          },
          axisLabel: {
            ...this.#axisLabelStyle,
            interval: Math.round(n / 50),
          },
          gridIndex: 1,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
          inverse: true,
        },
        {
          // Annotations.
          type: "category",
          data: annotationLabels,
          name: "Annotation",
          nameGap: 120,
          ...this.#axisStyle,
          axisLabel: {
            ...this.#axisLabelStyle,
          },
          gridIndex: 2,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
        },
        {
          // Contact detail.
          type: "category",
          data: ["", ""],
          name: "Residue Index",
          nameGap: 30,
          ...this.#axisStyle,
          axisLabel: {
            ...this.#axisLabelStyle,
          },
          gridIndex: 3,
          axisPointer: {
            show: true,
            label: { ...this.#axisPointerLabelStyle },
            triggerEmphasis: false,
            triggerTooltip: false,
          },
        },
      ],
      dataZoom: [
        {
          type: "inside",
          xAxisIndex: [0, 1, 2],
          throttle: 0,
        },
        {
          type: "inside",
          yAxisIndex: [1],
          throttle: 0,
        },
        {
          type: "inside",
          xAxisIndex: [3],
          throttle: 0,
        },
      ],
      legend: [
        {
          id: "legend",
          top: "top",
          left: "center",
          zlevel: 2,
          icon: "circle",
          itemWidth: 10,
          itemHeight: 10,
          orient: "horizontal",
          textStyle: {
            fontSize: 10,
            fontWeight: "lighter",
          },
          data: [],
          selector: [
            { type: "all", title: "Select all." },
            { type: "inverse", title: "Invert selection." },
          ],
          selectorLabel: {
            fontSize: 10,
            fontWeight: "lighter",
            borderRadius: 2,
          },
        },
      ],
      tooltip: [
        {
          trigger: "item",
          ...this.#tooltipStyle,
        },
      ],
      series: [
        ...Object.values(pscountSeries),
        ...Object.values(contactsSeries),
        ...Object.values(annotationSeries),
        {
          id: "contact_detail",
          name: "Contact detail",
          type: "scatter",
          data: [],
          xAxisIndex: 3,
          yAxisIndex: 3,
          cursor: "default",
          emphasis: {
            disabled: true,
          },
          markLine: {
            silent: true,
            label: {
              show: false,
            },
            lineStyle: {
              type: "solid",
              color: "#333333",
            },
            animation: false,
            symbol: [],
            data: [{ yAxis: 0 }, { yAxis: 1 }],
          },
        },
      ],
      visualMap: [],
    };
    // Append information from series definitions to option object.
    let legendData = new Set(
      this.#option.series
        .filter((_) => _.id.startsWith("contacts@"))
        .map((_) => {
          return _.name;
        })
    );
    this.#option.legend[0].data = [];
    legendData.forEach((_) => {
      this.#option.legend[0].data.push({ name: _ });
    });
    this.#option.tooltip[0].formatter = (params) => {
      let content = ``;
      let component = Array.isArray(params) ? params[0] : params;
      let positionIndexX;
      let positionIndexY;
      if (component.seriesId.startsWith("annotation"))
        // Interaction on position axis.
        positionIndexX = component.data.value[0];
      if (component.seriesId.startsWith("pscount"))
        // Interaction on position axis.
        positionIndexX = component.data[0];
      if (component.seriesId.startsWith("contacts")) {
        positionIndexX = component.data[0];
        positionIndexY = component.data[1];
      }
      // Fill content.
      let fillContent = (index) => {
        var noData;
        content +=
          "Position <code>" +
          (index + 1) +
          "&nbsp;" +
          this.#data.sequence[index] +
          "</code> ";
        content += `</br><small><b>PTM class counts</b></small></br>`;
        noData = true;
        for (const [cls, count] of Object.entries(
          this.#data.positions[index + 1].counts
        )) {
          noData = false;
          content +=
            `<small><code>` + count + `</code>&nbsp;` + cls + `</small></br>`;
        }
        if (noData) content += `<small>No data.</small></br>`;
        content += `<hr /><small><b>Position annotations</b></small></br>`;
        noData = true;
        for (const [cls, info] of Object.entries(
          this.#data.positions[index + 1].annotations
        )) {
          for (const i of info) {
            noData = false;
            content +=
              `<small><code>` +
              cls +
              `</code> at position ` +
              i[1] +
              ` to ` +
              i[2] +
              `&nbsp;` +
              i[0] +
              (i[3] != "" ? `&nbsp;` + i[3] : ``) +
              `</br>`;
          }
          content += `</small>`;
        }
        if (noData) content += `<small>No data.</small></br>`;
      };
      if (positionIndexX != undefined) fillContent(positionIndexX);
      content += "</br>";
      if (positionIndexY != undefined) fillContent(positionIndexY);
      if (positionIndexX != undefined && positionIndexY != undefined)
        content +=
          `<hr/><small>` +
          component.seriesName +
          `</small> <code>Click to highlight.</code>`;
      return content;
    };
    this.chart.instance.on("click", (params) => {
      let component = Array.isArray(params) ? params[0] : params;
      if (component.seriesId.startsWith("contacts")) {
        this.structure.highlightContacts(component.data[0] + 1);
        this.#showContactDetail(component.data[0] + 1, component.data[1] + 1);
      }
    });
  }

  fill(data) {
    this.#data = data;
    this.#postprocess();
    this.setModificationsOption();
    this.hideStructure();
    this.#updateOption(true);
  }

  #getColor(key) {
    if (!this.#colors.M.hasOwnProperty(key)) {
      this.#colors.M[key] = this.#colors.i;
      this.#colors.i += 1;
      if (this.#colors.i > this.#colors.C.length) {
        this.#colors.i = 0;
      }
    }
    return this.#colors.C[this.#colors.M[key]];
  }

  #postprocess() {
    if (this.#data == undefined) return;
    let _modifications = [];
    let _modificationsMassShift = [];
    this.#data.aminoacidCounts = {};
    this.#aminoAcids.forEach((aa) => (this.#data.aminoacidCounts[aa] = {}));
    this.#data.modificationCounts = {};
    // Add missing position entries.
    for (let p of [...Array(this.#data.sequence.length).keys()].map(
      (x) => x + 1
    )) {
      if (this.#data.positions.hasOwnProperty(p)) {
        this.#data.positions[p].counts = {};
        this.#data.positions[p].annotations = {};
      } else {
        this.#data.positions[p] = {
          modifications: [],
          counts: {},
          annotations: {},
        };
      }
    }
    // Extract count information from data.
    for (const [position, info] of Object.entries(this.#data.positions)) {
      const aa = this.#data.sequence[position - 1];
      for (const modification of Object.values(info.modifications)) {
        const modificationName = modification.display_name;
        const modificationClass = modification.modification_classification;
        _modifications.push(modificationName);
        if (!this.#data.aminoacidCounts[aa].hasOwnProperty(modificationName)) {
          this.#data.aminoacidCounts[aa][modificationName] = 1;
        } else {
          this.#data.aminoacidCounts[aa][modificationName] += 1;
        }
        if (!_modificationsMassShift.hasOwnProperty(modificationName))
          _modificationsMassShift[modificationName] = modification.mass_shift;
        if (!this.#data.modificationCounts.hasOwnProperty(modificationName))
          this.#data.modificationCounts[modificationName] = {};
        if (
          !this.#data.modificationCounts[modificationName].hasOwnProperty(
            modificationClass
          )
        ) {
          this.#data.modificationCounts[modificationName][
            modificationClass
          ] = 1;
        } else {
          this.#data.modificationCounts[modificationName][
            modificationClass
          ] += 1;
        }
        if (
          !this.#data.positions[position].counts.hasOwnProperty(
            modificationClass
          )
        ) {
          this.#data.positions[position].counts[modificationClass] = 1;
        } else {
          this.#data.positions[position].counts[modificationClass] += 1;
        }
      }
    }
    // Sort modification display names by mass shift and store in data.
    _modifications = [...new Set(_modifications)].sort((m1, m2) => {
      if (_modificationsMassShift[m1] < _modificationsMassShift[m2]) {
        return -1;
      } else if (_modificationsMassShift[m1] > _modificationsMassShift[m2]) {
        return 1;
      } else {
        return 0;
      }
    });
    this.#data.modifications = _modifications;
    // Extract annotation information and store in data.
    let annotationTypeToClass = {
      "Initiator methionine": "Molecule processing",
      Signal: "Molecule processing",
      "Transit peptide": "Molecule processing",
      Propeptide: "Molecule processing",
      Chain: "Molecule processing",
      Peptide: "Molecule processing",
      "Topological domain": "Region",
      Transmembrane: "Region",
      Intramembrane: "Region",
      Domain: "Region",
      Repeat: "Region",
      "Zinc finger": "Region",
      "DNA binding": "Region",
      Region: "Region",
      "Coiled coil": "Region",
      Motif: "Region",
      "Compositional bias": "Region",
      "Active site": "Site",
      "Binding site": "Site",
      Site: "Site",
      "Non-standard residue": "Amino acid modifications",
      "Modified residue": "Amino acid modifications",
      Lipidation: "Amino acid modifications",
      Glycosylation: "Amino acid modifications",
      "Disulfide bond": "Amino acid modifications",
      "Cross-link": "Amino acid modifications",
      "Alternative sequence": "Natural variations",
      "Natural variant": "Natural variations",
      Mutagenesis: "Experimental info",
      "Sequence uncertainty": "Experimental info",
      "Sequence conflict": "Experimental info",
      "Non-adjacent residues": "Experimental info",
      "Non-terminal residue": "Experimental info",
      Helix: "Secondary structure",
      Turn: "Secondary structure",
      "Beta strand": "Secondary structure",
    };
    for (let entry of Object.values(this.#data.annotation.features)) {
      if (entry.type == "Chain") continue;
      for (
        let position = parseInt(entry.location.start.value);
        position <= parseInt(entry.location.end.value);
        position++
      ) {
        let annotationClass = annotationTypeToClass[entry.type];
        if (
          !this.#data.positions[position].annotations.hasOwnProperty(
            annotationClass
          )
        )
          this.#data.positions[position].annotations[annotationClass] = [];
        this.#data.positions[position].annotations[annotationClass].push([
          entry.type,
          entry.location.start.value,
          entry.location.end.value,
          entry.description,
        ]);
      }
    }
  }

  #updateOption(replace) {
    if (this.#option != undefined) this.chart.setOption(this.#option, replace);
  }

  #showContactDetail(x, y) {
    if (this.#contentMode == 2) {
      let xModifications = {};
      this.#data.positions[x].modifications.forEach((o) => {
        xModifications[o.display_name] = o;
      });
      let yModifications = {};
      this.#data.positions[y].modifications.forEach((o) => {
        yModifications[o.display_name] = o;
      });
      if (
        Object.keys(xModifications).length == 0 &&
        Object.keys(yModifications).length == 0
      ) {
        this.#option.title[2].show = true;
      } else {
        this.#option.title[2].show = false;
      }
      var xData = {};
      for (const [name, _] of Object.entries(xModifications)) {
        let clr = this.#colors.explicit.contact_modified;
        let ms = _.mass_shift == "null" ? 0.0 : _.mass_shift;
        if (xData.hasOwnProperty(ms)) {
          xData[ms].name += "\n\t" + name;
          xData[ms].itemStyle.color = clr;
        } else {
          xData[ms] = {
            name: name,
            value: [ms, 0],
            itemStyle: { color: clr },
            symbol: "pin",
            symbolSize: 12,
            symbolRotate: 180,
            label: {
              show: true,
              ...this.#axisLabelStyle,
              position: [5, -5],
              rotate: 45,
              formatter: (params) => {
                return params.name;
              },
              overflow: "truncate",
            },
          };
        }
      }
      var yData = {};
      for (const [name, _] of Object.entries(yModifications)) {
        let clr = this.#colors.explicit.contact_modified;
        let ms = _.mass_shift == "null" ? 0.0 : _.mass_shift;
        if (yData.hasOwnProperty(ms)) {
          yData[ms].name += "\n\t" + name;
          yData[ms].itemStyle.color = clr;
        } else {
          yData[ms] = {
            name: name,
            value: [ms, 1],
            itemStyle: { color: clr },
            symbol: "pin",
            symbolSize: 12,
            label: {
              show: true,
              ...this.#axisLabelStyle,
              position: [5, -5],
              rotate: 45,
              formatter: (params) => {
                return params.name;
              },
              overflow: "truncate",
            },
          };
        }
      }
      this.#option.series.filter((_) => _.id == "contact_detail")[0].data = [
        ...Object.values(xData),
        ...Object.values(yData),
      ];
      this.#option.yAxis[3].data = [x, y];
      this.#updateOption();
    }
  }

  #contains(listA, listB) {
    const isContained = (entryA) => listB.includes(entryA);
    return listA.every(isContained);
  }
}

class StructureView {
  glviewer = null;
  domId = null;
  #contacts = null;
  #highlightPositions;
  #highlightModifications;
  #highlightContacts;
  #colors = {
    contact_unmodified: "#AAAAAA",
    contact_modified: "#000000",
    contact_highlight: "#DC5754",
  };

  constructor(id) {
    this.domId = id;
    this.glviewer = $3Dmol.createViewer($("#" + id), {
      backgroundColor: "#fbfbfb",
      antialias: true,
      cartoonQuality: 6,
    });
  }

  setStructure(pdbString, contacts) {
    this.#highlightPositions = null;
    this.#highlightModifications = null;
    this.#highlightContacts = null;
    this.glviewer.clear();
    this.glviewer.addModel(pdbString, "pdb");
    this.glviewer.zoomTo();
    this.glviewer.addSurface(
      "SAS",
      {
        color: "#d4d4d4",
        opacity: 0.4,
      },
      {}
    );
    this.style();
    this.glviewer.render();
    this.#contacts = contacts;
  }

  style() {
    this.glviewer.removeAllLabels();
    this.glviewer.removeAllShapes();
    this.glviewer.setStyle(
      {},
      {
        cartoon: {
          color: "#D4D4D4",
        },
      }
    );
    // Add style for contact highlight, if set.
    if (this.#highlightContacts != undefined) {
      let targetIndex = this.#highlightContacts;
      let targetCaAtom = this.glviewer.getAtomsFromSel({
        resi: [targetIndex],
        atom: "CA",
      })[0];
      this.glviewer.addStyle(
        {
          resi: [targetIndex],
        },
        {
          stick: {
            color: "#dd9ac2",
            radius: 0.7,
          },
        }
      );
      this.glviewer.addResLabels(
        {
          resi: [targetIndex],
        },
        {
          backgroundColor: "rgb(51, 51, 51)",
          backgroundOpacity: 0.7,
          fontColor: "#f0f5f5",
          fontSize: 10,
        }
      );
      if (this.#contacts.hasOwnProperty(targetIndex)) {
        for (let entry of this.#contacts[targetIndex]) {
          let contactIndex = entry[0];
          this.glviewer.addStyle(
            {
              resi: [contactIndex],
            },
            {
              stick: {
                color: "#b486ab",
                radius: 0.4,
              },
            }
          );
          this.glviewer.addResLabels(
            {
              resi: [contactIndex],
            },
            {
              backgroundColor: "rgb(51, 51, 51)",
              backgroundOpacity: 0.7,
              fontColor: "#f0f5f5",
              fontSize: 10,
            }
          );
          let contactCaAtom = this.glviewer.getAtomsFromSel({
            resi: [contactIndex],
            atom: "CA",
          })[0];
          this.glviewer.addLine({
            color: "#000000",
            hidden: false,
            dashed: false,
            start: {
              x: targetCaAtom.x,
              y: targetCaAtom.y,
              z: targetCaAtom.z,
            },
            end: {
              x: contactCaAtom.x,
              y: contactCaAtom.y,
              z: contactCaAtom.z,
            },
          });
          this.glviewer.addLabel(entry[1].toFixed(2) + " Å", {
            backgroundColor: "rgb(51, 51, 51)",
            backgroundOpacity: 0.4,
            fontColor: "#f0f5f5",
            fontSize: 10,
            position: {
              x: (targetCaAtom.x + contactCaAtom.x) / 2,
              y: (targetCaAtom.y + contactCaAtom.y) / 2,
              z: (targetCaAtom.z + contactCaAtom.z) / 2,
            },
          });
        }
      }
    }
    // Add style fo r position highlight, if set.
    if (
      this.#highlightPositions != undefined &&
      this.#highlightPositions.length > 0
    ) {
      this.glviewer.addStyle(
        { resi: this.#highlightPositions },
        {
          cartoon: {
            color: this.#colors.contact_highlight,
          },
        }
      );
      $("#panel-dashboard-structure-meta").html(
        `<i class="fa-solid fa-circle fa-sm"></i> Highlight positions with modifications: ` +
          this.#highlightModifications
            .map((_) => "<code>" + _ + "</code>")
            .join(", ")
      );
      $("#panel-dashboard-structure-meta").show();
    } else {
      $("#panel-dashboard-structure-meta").hide();
    }
    this.glviewer.render();
  }

  highlightContacts(targetIndex) {
    this.#highlightContacts = targetIndex;
    this.style();
  }

  highlightPositions(indices, modifications) {
    this.#highlightPositions = indices;
    this.#highlightModifications = modifications;
    this.style();
  }
}

/**
 * Inizializes all client side elements of the PTMVision application.
 */
function init() {
  __url = window.location.origin;
  $("#menu")[0].style.display = "flex";
  $("#panel")[0].style.display = "block";
  __overviewTable = new OverviewTable("panel-table-tabulator");
  __overviewTable.registerSelectionAction((data) => {
    if (data.length > 0) {
      $("#dashboard-display-button")[0].disabled = false;
    } else {
      $("#dashboard-display-button")[0].disabled = true;
    }
  });
  __overviewChart = new OverviewChart("panel-overview-chart");
  __dashboardChart = new DashboardChart(
    "panel-dashboard-chart",
    "panel-dashboard-structure"
  );
}

function startExampleSession() {
  displayNotification("Initializing example session.");
  axios
    .get(__url + "/example_session")
    .then((_) => {
      overviewTableInitialize(); // Init. table.
      overviewChartInitialize(); // Init.overview chart.
      togglePanel("panel-inputs");
    })
    .catch((error) => {
      console.error(error);
      displayAlert(error.message);
    })
    .finally(() => {
      removeNotification();
    });
}

async function startExistingSession() {
  displayNotification("Re-initializing session.");
  request = null;
  if ($("#session-input-form")[0].files.length == 0) {
    displayAlert("No session data was supplied.");
    $("body").css("cursor", "auto");
    removeNotification();
    return;
  }
  await readFile($("#session-input-form")[0].files[0]).then((response) => {
    request = response;
  });
  axios
    .post(__url + "/restart_session", request, {
      headers: {
        "Content-Type": "text",
      },
    })
    .then((_) => {
      overviewTableInitialize(); // Init. table.
      overviewChartInitialize(); // Init.overview chart.
      togglePanel("panel-inputs");
    })
    .catch((error) => {
      console.error(error);
      displayAlert(error.message);
    })
    .finally(() => {
      removeNotification();
    });
}

/**
 * Sends the specified search enginge output data to the PTMVision backend and loads the results in the overview table.
 */
async function startSession() {
  displayNotification("Transfer and process entered data.");
  request = {
    contentType: null,
    content: null,
  };
  if ($("#data-input-form")[0].files.length == 0) {
    displayAlert("No search engine output data was supplied.");
    $("body").css("cursor", "auto");
    removeNotification();
    return;
  }
  await readFile($("#data-input-form")[0].files[0]).then((response) => {
    request.content = response;
  });
  request.contentType = $("#data-type-form")[0].value;
  axios
    .post(
      __url + "/process_search_engine_output",
      pako.deflate(JSON.stringify(request)),
      {
        headers: {
          "Content-Type": "application/octet-stream",
          "Content-Encoding": "zlib",
        },
      }
    )
    .then((_) => {
      overviewTableInitialize(overviewChartInitialize); // Init. table and chart.
    })
    .catch((error) => {
      console.error(error);
      removeNotification();
      displayAlert(error.message);
    });
}

/**
 * Downloads the content of a blob to a client file.
 *
 * @param {Blob} blob File-like blob data to download.
 * @param {String} name The file name to save the downloaded content to.
 */
function downloadBlob(blob, name) {
  var download_link = document.createElement("a");
  download_link.href = window.URL.createObjectURL(new Blob([blob]));
  download_link.download = name;
  download_link.click();
  download_link.remove();
}

function downloadSessionData() {
  axios.get(__url + "/download_session").then((response) => {
    const D = new Date();
    downloadBlob(
      response.data,
      "ptmvis-" +
        [D.getFullYear(), D.getMonth() + 1, Date.now()].join("-") +
        ".zlib"
    );
  });
}

/**
 * Reads a file-like blob object to its String content.
 *
 * @param {Blob} file File-like blob data to read.
 * @returns String content of the file-like input blob.
 */
function readFile(file) {
  return new Promise((resolve, reject) => {
    var fileReader = new FileReader();
    fileReader.onload = (event) => {
      resolve(event.target.result);
    };
    fileReader.onerror = (error) => reject(error);
    fileReader.readAsText(file);
  });
}

/**
 * Redirects the browser to the about page of this project.
 */
function redirectAbout() {
  window.open(__url + "/about", "_blank");
}

function togglePanel(id) {
  $("#" + id).css("height", "50px");
  $("#" + id).css("overflow", "hidden");
  $("#" + id).css("color", "#d4d4d4");
  $("#" + id).css("cursor", "pointer");
  var minH = $("#" + id).css("min-height");
  $("#" + id).css("min-height", 0);
  if (id == "panel-overview") $("#panel-overview-chart").css("display", "none");
  document.getElementById(id).onclick = function () {
    $("#" + id).css("height", "auto");
    $("#" + id).css("overflow", "inherit");
    $("#" + id).css("color", "#333333");
    $("#" + id).css("cursor", "default");
    $("#" + id).css("min-height", minH);
    if (id == "panel-overview")
      $("#panel-overview-chart").css("display", "inherit");
  };
}

function displayNotification(text) {
  $("#menu").append(
    `<div class='notification'><i class="fa-duotone fa-spinner-third fa-spin fa-2xl"></i> ` +
      text +
      `</div>`
  );
}

function removeNotification() {
  $(".notification").remove();
}

/**
 * Displays a toast element with a custom error message to the user.
 *
 * @param {String} text The message to display.
 */
function displayAlert(text) {
  $("#menu").append(
    `<div class='alert'><i class="fa-duotone fa-circle-exclamation"></i> ` +
      text +
      `<button class='button-no-decoration float-right' onclick='$(".alert").remove()'><i class="fa-solid fa-x"></i></button></div>`
  );
}

function overviewTableInitialize(afterResponse) {
  axios
    .get(__url + "/get_available_proteins")
    .then((response) => {
      __overviewTable.setData(response.data);
      modifications_data_string = "";
      __overviewTable.availableModifications.forEach((m) => {
        modifications_data_string +=
          `<option value="` + m + `">` + m + `</option>`;
      });
      Metro.getPlugin("#panel-table-filter-modification", "select").data(
        modifications_data_string
      );
      identifiers_data_string = "";
      __overviewTable.availableIdentifiers.forEach((i) => {
        let entry = i.split("$");
        identifiers_data_string +=
          `<option value="id$` +
          entry[0] +
          `">` +
          entry[0] +
          `</option><option value="name$` +
          entry[1] +
          `">` +
          entry[1] +
          `</option>`;
      });
      Metro.getPlugin("#panel-table-filter-id", "select").data(
        identifiers_data_string
      );
      $("#panel-table-title").html(
        "Select single protein of interest (" +
          response.data.length +
          " available)"
      );
      if (afterResponse != undefined) afterResponse();
    })
    .catch((error) => {
      console.error(error);
      removeNotification();
      displayAlert(error.message);
    });
}

function overviewTableSetFilters() {
  __overviewTable.setFilters(
    Metro.getPlugin("#panel-table-filter-id", "select").val(),
    Metro.getPlugin("#panel-table-filter-modification", "select").val()
  );
}

function overviewChartInitialize(afterResponse) {
  axios
    .get(__url + "/get_overview_data")
    .then((response) => {
      __overviewChart.fill(response.data);
      window.scrollTo({
        top: $("#panel-overview").get(0).offsetTop + 40,
        behavior: "smooth",
      });
      if (afterResponse != undefined) afterResponse();
    })
    .catch((error) => {
      console.error(error);
      removeNotification();
      displayAlert(error.message);
    })
    .finally(() => {
      removeNotification();
      togglePanel("panel-inputs");
    });
}

function overviewChartDownloadImage() {
  const _ = document.createElement("a");
  document.body.appendChild(_);
  _.setAttribute("download", "overview.png");
  _.href = __overviewChart.getDataUrl();
  if (!_.href.endsWith("undefined")) _.click();
  _.remove();
}

function overviewChartRestoreZoom() {
  __overviewChart.restoreZoom();
}

function overviewChartHighlight() {
  let options = [];
  let names = __overviewChart.getDataNames();
  for (let i = 0; i < names.length; i++) {
    options.push(`<option value="` + i + `">` + names[i] + `</option>`);
  }
  Swal.fire({
    backdrop: false,
    confirmButtonColor: "#62a8ac",
    width: "auto",
    padding: "1em",
    position: "center",
    html:
      `<h4><small>Select PTMs to highlight:</small></h4>
    <select id="tmpSelect" data-role="select" class="input-small" multiple>` +
      options.join("") +
      `</select>`,
  }).then((result) => {
    if (result.isConfirmed) {
      __overviewChart.highlight(
        Metro.getPlugin("#tmpSelect", "select")
          .val()
          .map((_) => parseInt(_))
      );
    }
  });
}

function overviewChartSort() {
  __overviewChart.resort();
}

function dashboardChartInitialize(cutoff_value, pdb_text_value) {
  displayNotification("Initializing dashboard.");
  if (cutoff_value == undefined) cutoff_value = 4.69;
  if (pdb_text_value == undefined) pdb_text_value = null;
  let sel_protein_id = __overviewTable.getSelection();
  if (sel_protein_id == null) {
    removeNotification();
    displayAlert(
      `No protein was selected from panel <i class="fa-duotone fa-circle-3"></i>`
    );
    return;
  }
  request = {
    uniprot_pa: sel_protein_id,
    pdb_text: pdb_text_value,
    cutoff: cutoff_value,
  };
  axios
    .post(__url + "/get_protein_data", pako.deflate(JSON.stringify(request)), {
      headers: {
        "Content-Type": "application/octet-stream",
        "Content-Encoding": "zlib",
      },
    })
    .then((response) => {
      if (response.status == 303) {
        displayAlert("Failed to fetch protein structure.");
      } else {
        var protein_data = response.data;
        let primaryAccession = protein_data.annotation.primaryAccession;
        let proteinName = protein_data.annotation.hasOwnProperty(
          "proteinDescription"
        )
          ? protein_data.annotation.proteinDescription.hasOwnProperty(
              "recommendedName"
            )
            ? protein_data.annotation.proteinDescription.recommendedName
                .fullName.value
            : Object.values(protein_data.annotation.proteinDescription)[0][0]
                .fullName.value
          : "N/A";
        let organismName = protein_data.annotation.hasOwnProperty("organism")
          ? "<em>" + protein_data.annotation.organism.scientificName + "</em>"
          : "N/A";
        $("#panel-dashboard-selection").html(
          [primaryAccession, proteinName, organismName].join(", ")
        );
        $("#panel-dashboard-selection").css("cursor", "pointer");
        $("#panel-dashboard-selection").on("click", () =>
          Swal.fire({
            html: protein_data.annotation.comments
              .filter((_) => {
                return _.hasOwnProperty("texts");
              })
              .map((_) => {
                return (
                  "<code>" +
                  _.commentType +
                  "</code><p class='text-left'>" +
                  _.texts[0].value +
                  "</p>"
                );
              })
              .join("</br>"),
            confirmButtonColor: "#d4d4d4",
            confirmButtonText: "Close Info",
          })
        );
        __dashboardChart.fill(protein_data);
        __dashboardContent = "modifications";
        $("#panel-dashboard-title").html("Explore detail - Modifications view");
        window.scrollTo({
          top: $("#panel-dashboard").get(0).offsetTop + 40,
          behavior: "smooth",
        });
        togglePanel("panel-overview");
        togglePanel("panel-table");
      }
    })
    .catch((error) => {
      console.error(error);
      displayAlert(error.message);
    })
    .finally(() => {
      removeNotification();
    });
}

function dashboardChartDownloadImage() {
  const _ = document.createElement("a");
  document.body.appendChild(_);
  _.setAttribute("download", "dashboard.png");
  _.href = __dashboardChart.getDataUrl();
  if (!_.href.endsWith("undefined")) _.click();
  _.remove();
}

function dashboardChartRestoreZoom() {
  __dashboardChart.restoreZoom();
}

function dashboardChartHighlight() {
  let options = [];
  let names = __dashboardChart.getDataNames();
  for (let i = 0; i < names.length; i++) {
    options.push(`<option value="` + i + `">` + names[i] + `</option>`);
  }
  Swal.fire({
    backdrop: false,
    confirmButtonColor: "#62a8ac",
    width: "auto",
    padding: "1em",
    position: "center",
    html:
      `<h4><small>Select PTMs to highlight:</small></h4>
    <select id="tmpSelect" data-role="select" class="input-small" multiple>` +
      options.join("") +
      `</select>`,
  }).then((result) => {
    if (result.isConfirmed) {
      __dashboardChart.highlight(Metro.getPlugin("#tmpSelect", "select").val());
    }
  });
}

function dashboardChartSwitch() {
  if (__dashboardContent == "modifications") {
    __dashboardContent = "structure";
    $("#panel-dashboard-title").html("Explore detail - Structure view");
  } else if (__dashboardContent == "structure") {
    __dashboardContent = "modifications";
    $("#panel-dashboard-title").html("Explore detail - Modifications view");
  }
  __dashboardChart.switchContent();
}

function to(id) {
  let e = $("#" + id).get(0);
  let offset = Math.round(e.getBoundingClientRect().top) - 80;
  window.scrollTo({
    top: offset,
    behavior: "smooth",
  });
}

function stringToHash(string) {
  let hash = 0;

  if (string.length == 0) return hash;

  for (i = 0; i < string.length; i++) {
    let char = string.charCodeAt(i);
    hash = (hash << 5) - hash + char;
    hash = hash & hash;
  }

  return hash;
}
