/*!
 * PTMVision v1.1 (https://github.com/Integrative-Transcriptomics/PTMVision)
 * Copyright 2023-2024 by Caroline Jachmann and Simon Hackl
 * Licensed under GPL-3.0
 !*/

/**
 * Instance of the OverviewTable class.
 */
var _overviewTable = null;
/**
 * Instance of the OverviewChart class.
 */
var _overviewChart = null;
/**
 * Instance of the DashboardChart class.
 */
var _dashboardChart = null;
/**
 * The currently displayed dashboard content. One of "modifications" or "structure".
 */
var _dashboardContent = null;

/**
 * ECharts option for global axis style.
 */
const STYLE_AXIS = {
  nameLocation: "center",
  nameGap: 40,
  nameTextStyle: {
    fontWeight: "bold",
    fontSize: 14,
  },
  axisTick: {
    alignWithLabel: true,
    length: 4,
    interval: 0,
  },
};
/**
 * ECharts option for global axis label style.
 */
const STYLE_AXIS_LABEL = {
  fontWeight: "normal",
  fontSize: 12,
};
/**
 * ECharts option for global tooltip style.
 */
const STYLE_TOOLTIP = {
  backgroundColor: "#f0f5f5",
  borderColor: "#d4d4d4",
  textStyle: {
    color: "#333333",
    fontSize: 12,
  },
};
/**
 * ECharts option for global title style.
 */
const STYLE_TITLE = {
  textStyle: {
    fontSize: 14,
    fontWeight: "bold",
  },
};
/**
 * ECharts option for global pointer style.
 */
const STYLE_POINTER = {
  label: {
    show: true,
    fontWeight: "bold",
    fontSize: 11,
    color: "#333333",
    padding: [2, 4, 2, 4],
    backgroundColor: "#fbfbfbe6",
    borderColor: "#fbfbfb",
    margin: 1,
  },
};

/**
 * Internal class for handling the overview table component.
 */
class OverviewTable {
  /**
   * Tabulator instance.
   */
  tabulator = null;
  /**
   * Available identifiers for filtering.
   */
  availableIdentifiers = [];
  /**
   * Available modifications for filtering.
   */
  availableModifications = [];

  /**
   * Class constructor.
   *
   * @param {String} id DOM id of the element to bind the component to.
   */
  constructor(id) {
    this.tabulator = new Tabulator("#" + id, {
      selectableRows: 1,
      columns: [
        {
          title: "ID",
          field: "id",
          sorter: "string",
          width: "12%",
        },
        {
          title: "Name",
          field: "name",
          sorter: "string",
          width: "32%",
        },
        {
          title: "Primary Sequence Length",
          field: "length",
          sorter: "number",
          width: "18%",
        },
        {
          title: "No. modified positions",
          field: "modified_positions",
          sorter: "number",
          width: "18%",
        },
        {
          title: "No. distinct modifications",
          field: "unique_modifications",
          sorter: "number",
          width: "18%",
        },
        {
          title: "Modifications",
          field: "modifications",
          visible: false,
          download: false,
        },
      ],
    });
  }

  /**
   * Exposes the Tabulator rowSelectionChanged event handler.
   *
   * @param {Function} action Function to be executed on row selection change.
   */
  registerSelectionAction(action) {
    this.tabulator.on(
      "rowSelectionChanged",
      function (data, rows, selected, deselected) {
        action(data);
      }
    );
  }

  /**
   * Sets the filters for the table.
   *
   * @param {Array} idValues List of identifiers to filter for.
   * @param {Array} modValues List of modifications to filter for.
   */
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

  /**
   * Fills the table with data.
   *
   * @param {Array} data List of Objects containing the data to be displayed.
   */
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

  /**
   * Returns the selected row's ID.
   *
   * @returns The selected row's ID.
   */
  getSelection() {
    let _ = this.tabulator.getSelectedData();
    if (_.length > 0) return this.tabulator.getSelectedData()[0].id;
    else return null;
  }
}

/**
 * Internal class for handling EChart instances.
 */
class Chart {
  /**
   * EChart instance.
   */
  instance = null;
  /**
   * DOM id of the element to bind the chart to.
   */
  #instanceDomId = null;
  /**
   * ResizeObserver instance for the chart.
   */
  #resizeObserver = null;

  /**
   * Class constructor.
   *
   * @param {String} id DOM id of the element to bind the component to.
   */
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

  /**
   *
   * @param {Object} option An ECharts option object (https://echarts.apache.org/en/option.html).
   * @param {Boolean} replace Handles the replace/merge behavior of the setOption method (https://echarts.apache.org/en/api.html#echartsInstance.setOption).
   */
  setOption(option, replace) {
    this.instance.setOption(option, replace);
  }
}

/**
 * Internal class for handling the overview chart component.
 */
class OverviewChart {
  /**
   * Chart instance.
   */
  chart = null;
  /**
   * Internal data object.
   */
  #data;
  /**
   * Internal ECharts option object.
   */
  #option;
  /**
   * Internal sorting index.
   */
  #sortingIndex = 1;

  /**
   * Class constructor.
   *
   * @param {String} id DOM id of the element to bind the component to.
   */
  constructor(id) {
    this.chart = new Chart(id);
  }

  /**
   * Restores the zoom level of the chart.
   *
   * @returns None if no option is set.
   */
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

  /**
   * Highlights the provided indices in the chart.
   *
   * @param {Array} indices
   * @returns None if no option is set.
   */
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
    this.#option.series[2].markLine = {
      ...markLineOption,
      data: indices.map((i) => {
        return { yAxis: i };
      }),
    };
    this.#updateOption(false);
  }

  /**
   * Switches the sorting of the chart.
   *
   * @returns None if no option is set.
   */
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

  /**
   * Returns the modification names of the data wrt. to the current sorting index.
   *
   * @returns The names of the data or an empty list if no data is set.
   */
  getDataNames() {
    if (this.#data != undefined)
      return this.#data.modificationNamesSorted[this.#sortingIndex];
    else return [];
  }

  /**
   * Returns the data URL of the chart used for exporting images.
   *
   * @returns The data URL of the chart or None if no option is set.
   */
  getDataUrl() {
    if (this.#option == undefined) return;
    return this.chart.instance.getDataURL({
      pixelRatio: 4,
      backgroundColor: "#fff",
    });
  }

  /**
   * Fills the chart with data; I.e., generates the ECharts option object from the data.
   *
   * @param {Object} data Contains the data to be displayed.
   */
  fill(data) {
    if (data != undefined)
      this.#data = {
        modifications: data[0],
        modificationNamesSorted: data[1],
        coOccurrence: data[2],
        classCounts: Object.entries(data[3]) // Sort class counts by value.
          .sort(([, v1], [, v2]) => v2 - v1)
          .reduce((acc, [key, value]) => ({ ...acc, [key]: value }), {}),
        meta: data[4],
      };
    let r = this.chart.instance.getWidth() / this.chart.instance.getHeight(); // Approximate quadratic aspect ratio.
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
        if (this.#data.coOccurrence.hasOwnProperty(_))
          dataCoOccurrence.push([i, j, this.#data.coOccurrence[_]]);
      }
    }
    // Fill in option.
    this.#option = {
      title: [
        {
          text:
            "Total modifications: " +
            Object.values(this.#data.classCounts).reduce((a, b) => a + b, 0) +
            " | Distinct modification types: " +
            Object.keys(this.#data.modifications).length +
            " | Distinct UniMod modification classes: " +
            Object.keys(this.#data.classCounts).length,
          top: "top",
          left: "center",
          ...STYLE_TITLE,
        },
        {
          text: "Shared PTM Sites between Modification Types",
          top: 35,
          left: "10%",
          ...STYLE_TITLE,
        },
        {
          text: "Mass shift in Dalton",
          top: 54,
          left: 10 + 1 + 65 / r + "%",
          ...STYLE_TITLE,
        },
        {
          text: "Site count",
          top: 54,
          left: 10 + 1 + 15 + 1 + 65 / r + "%",
          ...STYLE_TITLE,
        },
        {
          text: "Unimod PTM class counts",
          top: 54,
          left: 10 + 1 + 15 + 1 + 15 + 1 + 5 + 65 / r + "%",
          ...STYLE_TITLE,
        },
      ],
      grid: [
        {
          top: 80,
          left: "10%",
          height: "65%",
          width: 65 / r + "%",
          show: true,
        },
        {
          top: 80,
          left: 10 + 1 + 65 / r + "%",
          height: "65%",
          width: "15%",
          show: true,
        },
        {
          top: 80,
          left: 10 + 1 + 15 + 1 + 65 / r + "%",
          height: "65%",
          width: "15%",
          show: true,
        },
        {
          top: 80,
          left: 10 + 1 + 15 + 1 + 15 + 1 + 5 + 65 / r + "%",
          right: "1%",
          height: "60%",
          width: "auto",
          show: true,
        },
      ],
      xAxis: [
        {
          gridIndex: 0,
          type: "category",
          name: "Modification",
          data: dataAxis,
          ...STYLE_AXIS,
          nameGap: 100,
          axisLabel: {
            show: true,
            rotate: 50,
            formatter: (i) => {
              let displayName = this.#data.modifications[i]["display_name"];
              displayName = displayName.replace(
                "Unannotated mass-shift",
                "Mass-shift"
              );
              return displayName.length > 12
                ? displayName.slice(0, 12) + "..."
                : displayName;
            },
            ...STYLE_AXIS_LABEL,
          },
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
          },
        },
        {
          gridIndex: 1,
          type: "value",
          name: "Mass shift [Da]",
          ...STYLE_AXIS,
          nameGap: 100,
          axisLabel: {
            show: true,
            interval: 0,
            rotate: 50,
            ...STYLE_AXIS_LABEL,
          },
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
          },
        },
        {
          gridIndex: 2,
          type: "value",
          name: "Count",
          ...STYLE_AXIS,
          nameGap: 100,
          axisLabel: {
            show: true,
            interval: 0,
            rotate: 50,
            ...STYLE_AXIS_LABEL,
          },
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
          },
        },
        {
          gridIndex: 3,
          type: "category",
          name: "Class",
          ...STYLE_AXIS,
          nameGap: 130,
          data: Object.keys(this.#data.classCounts),
          axisLabel: {
            show: true,
            interval: 0,
            rotate: 50,
            ...STYLE_AXIS_LABEL,
          },
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            ...STYLE_POINTER,
          },
        },
      ],
      yAxis: [
        {
          gridIndex: 0,
          type: "category",
          name: "Modification",
          ...STYLE_AXIS,
          nameGap: 120,
          data: dataAxis,
          inverse: true,
          axisLabel: {
            show: true,
            formatter: (i) => {
              let displayName = this.#data.modifications[i]["display_name"];
              displayName = displayName.replace(
                "Unannotated mass-shift",
                "Mass-shift"
              );
              return displayName.length > 12
                ? displayName.slice(0, 12) + "..."
                : displayName;
            },
            ...STYLE_AXIS_LABEL,
          },
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
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
            triggerEmphasis: false,
            ...STYLE_POINTER,
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
            triggerEmphasis: false,
            ...STYLE_POINTER,
          },
        },
        {
          gridIndex: 3,
          type: "value",
          name: "Count",
          ...STYLE_AXIS,
          axisLabel: {
            show: true,
            ...STYLE_AXIS_LABEL,
          },
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
          },
        },
      ],
      tooltip: {
        ...STYLE_TOOLTIP,
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
        dataCoOccurrence.length > 0
          ? {
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
              top: 54,
              left: "10%",
              itemHeight: 50,
              itemWidth: 11,
              text: [
                Math.max(...Object.values(this.#data.coOccurrence)),
                "No. shared sites 1",
              ],
              textStyle: { fontWeight: "normal", fontSize: 12 },
            }
          : null,
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
            borderWidth: 0.2,
            borderRadius: 2,
            borderColor: "#fbfbfb",
          },
          data: dataCoOccurrence,
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
    // Highlight PTMs that fall within mass shift tolerance, when sorting by mass shift.
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
        if (
          dataMassShift[i] + this.#data.meta.mass_shift_tolerance >=
          dataMassShift[j]
        ) {
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

  /**
   * Clears the chart option.
   */
  clear() {
    this.chart.instance.clear();
  }

  /**
   * Internal method for updating the chart option.
   *
   * @param {Boolean} replace Passed to the setOption method of the Chart instance.
   */
  #updateOption(replace) {
    if (this.#option != undefined)
      this.chart.setOption(this.#option, { notMerge: replace });
  }
}

/**
 * Internal class for handling the modification- and structure view components.
 */
class DashboardChart {
  /**
   * Chart instance.
   */
  chart = null;
  /**
   * StructureView instance.
   */
  structure = null;
  /**
   * Internal data object.
   */
  #data;
  /**
   * Internal ECharts option object.
   */
  #option;
  /**
   * Internal color object.
   */
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
  /**
   * Internal amino acid three letter code list.
   */
  #aminoAcids3 = [
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
  /**
   * Internal amino acid one letter code list.
   */
  #aminoAcids1 = [
    "A", // "ALA",
    "C", // "CYS",
    "D", // "ASP",
    "E", // "GLU",
    "F", // "PHE",
    "G", // "GLY",
    "H", // "HIS",
    "I", // "ILE",
    "K", // "LYS",
    "L", // "LEU",
    "M", // "MET",
    "N", // "ASN",
    "P", // "PRO",
    "Q", // "GLN",
    "R", // "ARG",
    "S", // "SER",
    "T", // "THR",
    "V", // "VAL",
    "W", // "TRP",
    "Y", // "TYR",
  ];
  /**
   * The current content mode. One of 1 (modifications) or 2 (structure).
   */
  #contentMode = 1;

  /**
   * Class constructor.
   *
   * @param {String} id1 DOM id of the element to bind the chart component to.
   * @param {String} id2 DOM id of the element to bind the 3DMol structure view component to.
   */
  constructor(id1, id2) {
    this.chart = new Chart(id1);
    this.structure = new StructureView(id2);
    this.hideStructure();
  }

  /**
   * Shows the 3DMol structure view component.
   */
  showStructure() {
    $("#" + this.structure.domId).show();
  }

  /**
   * Hides the 3DMol structure view component.
   */
  hideStructure() {
    $("#" + this.structure.domId).hide();
  }

  /**
   * Restores the zoom level of the chart.
   *
   * @returns None if no option is set.
   */
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

  /**
   * Highlights the provided indices in the chart.
   *
   * @param {Array} indices
   * @returns None if no option is set.
   */
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
          fontSize: 9,
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

  /**
   * Switches the content mode of the component.
   *
   * @returns None if no option is set.
   */
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

  /**
   * Returns the modification names of the data wrt. to the current sorting index.
   *
   * @returns The names of the data or an empty list if no data is set.
   */

  getDataNames() {
    if (this.#data != undefined) return this.#data.modifications;
    else return [];
  }

  /**
   * Returns the data URL of the chart used for exporting images.
   *
   * @returns The data URL of the chart or None if no option is set.
   */
  getDataUrl() {
    if (this.#option == undefined) return;
    return this.chart.instance.getDataURL({
      pixelRatio: 4,
      backgroundColor: "#fff",
    });
  }

  /**
   * Generates the ECharts option object for the modifications view from the currently set data.
   */
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
        borderRadius: 1,
        borderColor: "#fbfbfb",
      },
      cursor: "default",
      emphasis: {
        disabled: true,
      },
    };
    let minMaxScaleCounts = (count, aa) => {
      let min = Math.min(
        ...Object.values(this.#data.aminoacidCounts[aa]).map((_) => _[0])
      );
      let max = Math.max(
        ...Object.values(this.#data.aminoacidCounts[aa]).map((_) => _[0])
      );
      if (min == max) return 1;
      else return (count - min) / (max - min);
    };
    for (const [aa, _] of Object.entries(this.#data.aminoacidCounts)) {
      for (const [mdname, info] of Object.entries(_)) {
        aacountSeries.data.push([
          this.#aminoAcids3.indexOf(aa),
          this.#data.modifications.indexOf(mdname),
          minMaxScaleCounts(info[0], aa),
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
              borderRadius: 1,
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
        {
          text: "Per Position PTMs and Annotations",
          top: "2%",
          left: "8%",
          ...STYLE_TITLE,
        },
        {
          text: "Per Aminoacid PTM Counts",
          top: "12%",
          left: "64%",
          ...STYLE_TITLE,
        },
        {
          text: "Per PTM Class Counts",
          top: "14%",
          left: "85%",
          ...STYLE_TITLE,
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
          width: "12%",
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
            triggerEmphasis: false,
            triggerTooltip: true,
            ...STYLE_POINTER,
          },
        },
        {
          // Modifications map.
          data: Object.keys(this.#data.positions),
          show: false,
          gridIndex: 1,
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
          },
        },
        {
          // Annotations.
          data: Object.keys(this.#data.positions),
          name: "Protein Position",
          ...STYLE_AXIS,
          nameGap: 35,
          axisLabel: {
            interval: Math.round(n / 40),
          },
          gridIndex: 2,
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: true,
            ...STYLE_POINTER,
          },
        },
        {
          // Amino acids.
          data: this.#aminoAcids3,
          name: "Amino Acid",
          ...STYLE_AXIS,
          nameGap: 35,
          gridIndex: 3,
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
          },
          axisLabel: {
            rotate: 60,
            interval: 0,
          },
        },
        {
          // Modification counts.
          type: "value",
          name: "Count",
          ...STYLE_AXIS,
          nameGap: 35,
          gridIndex: 4,
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
          },
        },
      ],
      yAxis: [
        {
          // Position counts.
          type: "value",
          name: "Count",
          ...STYLE_AXIS,
          gridIndex: 0,
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
          },
        },
        {
          // Modifications map.
          type: "category",
          data: this.#data.modifications,
          name: "Modification",
          ...STYLE_AXIS,
          nameGap: 130,
          axisLabel: {
            formatter: (value) => {
              value = value.replace("Unannotated mass-shift", "Mass-shift");
              return value.length > 17 ? value.slice(0, 17) + "..." : value;
            },
          },
          gridIndex: 1,
          inverse: true,
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
          },
        },
        {
          // Annotations.
          // type: "category",
          data: annotationLabels,
          name: "Annotation",
          ...STYLE_AXIS,
          nameGap: 130,
          gridIndex: 2,
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
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
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
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
            triggerEmphasis: false,
            triggerTooltip: true,
            ...STYLE_POINTER,
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
            fontSize: 11,
            fontWeight: "lighter",
          },
          data: [],
          selector: [
            { type: "all", title: "Select all." },
            { type: "inverse", title: "Invert selection." },
          ],
          selectorLabel: {
            fontSize: 11,
            fontWeight: "lighter",
            borderRadius: 2,
          },
        },
      ],
      tooltip: [
        {
          trigger: "item",
          ...STYLE_TOOLTIP,
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
        fontSize: 11,
        fontWeight: "normal",
      },
      text: ["1", "Per Aminoacid Min-Max Scaled No. Occurrence 0"],
      min: 0,
      max: 1,
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
      if (aminoacidIndex != undefined) {
        let countAndClass =
          this.#data.aminoacidCounts[this.#aminoAcids3[aminoacidIndex]][
            this.#data.modifications[modificationIndex]
          ];
        contentHead +=
          "has " +
          countAndClass[0] +
          " occurrences (<code>" +
          countAndClass[1] +
          "</code>) on aminoacid <code>" +
          this.#aminoAcids3[aminoacidIndex] +
          "</code>";
        return contentHead; // Do not fill content body for aminoacid axis.
      }
      // Fill content body based on axis.
      var noData;
      if (positionIndex != undefined) {
        noData = true;
        contentBody += `<hr /><b>Position PTM class counts</b></br>`;
        for (const [cls, count] of Object.entries(
          this.#data.positions[positionIndex + 1].counts
        )) {
          noData = false;
          contentBody += `<code>` + count + `</code>&nbsp;` + cls + `</br>`;
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
              `<code>` +
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
        }
        if (noData) contentBody += `No data.`;
      }
      if (modificationIndex != undefined) {
        noData = true;
        contentBody += `<hr /><b>PTM class counts</b></br>`;
        for (const [cls, count] of Object.entries(
          this.#data.modificationCounts[
            this.#data.modifications[modificationIndex]
          ]
        )) {
          noData = false;
          contentBody += `<code>` + count + `</code>&nbsp;` + cls + `</br>`;
        }
        if (noData) contentBody += `No data.`;
      }
      return contentHead + contentBody;
    };
    this.chart.instance.on("click", () => {});
  }

  /**
   * Generates the ECharts option object for the structure view from the currently set data.
   *
   * @param {Boolean} setStructure Whether to initialize the structure view.
   * @param {Array} highlightModifications  Modifications to highlight in the components.
   */
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
    for (let [residueI, contactEntries] of Object.entries(
      this.#data.contacts
    )) {
      let x = parseInt(residueI);
      for (let [residueJ, _] of Object.entries(contactEntries)) {
        let y = parseInt(residueJ);
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
        contactsSeries[cls].data.push([x - 1, y - 1, 1]); // -1 to account for 0-based index of series.
      }
    }
    // Initialize structure view.
    if (setStructure)
      this.structure.setStructure(this.#data.structure, this.#data.contacts);
    // Construct option object.
    let r = this.chart.instance.getWidth() / this.chart.instance.getHeight();
    let n = this.#data.sequence.length;
    $("#panel-dashboard-structure").css("left", 12 + 75 / r + "%");
    $("#panel-dashboard-structure").css("width", 100 - (16 + 75 / r) + "%");
    this.#option = this.#option = {
      title: [
        {
          text: "Per Position PTM Counts, Residue Contacts and Annotations",
          top: "4%",
          left: "8%",
          ...STYLE_TITLE,
        },
        {
          text: "Protein Structure and Residue Contact PTM Details",
          top: "4%",
          left: 12 + 75 / r + "%",
          ...STYLE_TITLE,
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
          left: 12 + 75 / r + "%",
          width: 100 - (16 + 75 / r) + "%",
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
            triggerEmphasis: false,
            triggerTooltip: true,
            ...STYLE_POINTER,
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
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
          },
        },
        {
          // Annotations.
          type: "category",
          data: Object.keys(this.#data.positions),
          name: "Protein Position",
          ...STYLE_AXIS,
          nameGap: 35,
          gridIndex: 2,
          axisLabel: {
            interval: Math.round(n / 25),
          },
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: true,
            ...STYLE_POINTER,
          },
        },
        {
          // Contact detail.
          type: "value",
          name: "Mass Shift [Da]",
          ...STYLE_AXIS,
          nameGap: 35,
          gridIndex: 3,
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
          },
        },
      ],
      yAxis: [
        {
          // Position counts.
          type: "value",
          name: "Count",
          ...STYLE_AXIS,
          gridIndex: 0,
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
          },
          intervaL: 4,
        },
        {
          // Contact map.
          type: "category",
          data: Object.keys(this.#data.positions),
          name: "Position",
          ...STYLE_AXIS,
          axisLabel: {
            interval: Math.round(n / 40),
          },
          gridIndex: 1,
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
          },
          inverse: true,
        },
        {
          // Annotations.
          type: "category",
          data: annotationLabels,
          name: "Annotation",
          ...STYLE_AXIS,
          nameGap: 125,
          gridIndex: 2,
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
          },
        },
        {
          // Contact detail.
          type: "category",
          data: ["", ""],
          name: "Residue Index",
          ...STYLE_AXIS,
          gridIndex: 3,
          axisPointer: {
            show: true,
            triggerEmphasis: false,
            triggerTooltip: false,
            ...STYLE_POINTER,
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
            fontSize: 11,
            fontWeight: "normal",
          },
          data: [],
          selector: [
            { type: "all", title: "Select all." },
            { type: "inverse", title: "Invert selection." },
          ],
          selectorLabel: {
            fontSize: 11,
            fontWeight: "lighter",
            borderRadius: 2,
          },
        },
      ],
      tooltip: [
        {
          trigger: "item",
          ...STYLE_TOOLTIP,
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
        content += `</br><b>PTM class counts</b></br>`;
        noData = true;
        for (const [cls, count] of Object.entries(
          this.#data.positions[index + 1].counts
        )) {
          noData = false;
          content += `<code>` + count + `</code>&nbsp;` + cls + `</br>`;
        }
        if (noData) content += `No data.</br>`;
        content += `<hr /><b>Position annotations</b></br>`;
        noData = true;
        for (const [cls, info] of Object.entries(
          this.#data.positions[index + 1].annotations
        )) {
          for (const i of info) {
            noData = false;
            content +=
              `<code>` +
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
        }
        if (noData) content += `No data.</br>`;
      };
      if (positionIndexX != undefined) fillContent(positionIndexX);
      content += "</br>";
      if (positionIndexY != undefined) fillContent(positionIndexY);
      if (positionIndexX != undefined && positionIndexY != undefined)
        content +=
          `<hr/>` + component.seriesName + ` <code>Click to highlight.</code>`;
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

  /**
   * Fills the chart with data; I.e., generates the ECharts option object from the data.
   *
   * @param {Object} data Contains the data to be displayed.
   */
  fill(data) {
    this.#data = data;
    this.#postprocess();
    this.setModificationsOption();
    this.hideStructure();
    this.#updateOption(true);
  }

  /**
   * Clears the chart option.
   */
  clear() {
    this.chart.instance.clear();
    this.hideStructure();
    $("#panel-dashboard-selection").html("No Protein Selected");
    $("#panel-dashboard-selection").css("cursor", "default");
    $("#panel-dashboard-selection").prop("onclick", null).off("click");
  }

  /**
   * Retrieves the set color for a given key. If the key is not set, a new color is assigned.
   *
   * @param {String} key Identifier for any colorable component.
   * @returns The color for the given key.
   */
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

  /**
   * Postprocesses the data to generate additional information for the chart.
   *
   * @returns None if no data is set.
   */
  #postprocess() {
    if (this.#data == undefined) return;
    let _modifications = [];
    let _modificationsMassShift = [];
    this.#data.aminoacidCounts = {};
    this.#aminoAcids3.forEach((aa) => (this.#data.aminoacidCounts[aa] = {}));
    this.#data.modificationCounts = {};
    // Add missing position entries.
    // Check for positions wrt. the sequence length that are not reflected in the data.
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
    // Check for positions wrt. the data that are not reflected in the sequence length.
    for (const [position, _] of Object.entries(this.#data.positions)) {
      if (
        parseInt(position) < 1 ||
        parseInt(position) > this.#data.sequence.length
      ) {
        // Remove from contact information.
        if (this.#data.contacts.hasOwnProperty(position))
          delete this.#data.contacts[position];
        // Remove modifications information.
        delete this.#data.positions[position];
        console.warn(
          "Delete modifications data at position " +
            position +
            "; Position exceeds the sequence length or is smaller than one."
        );
      }
    }
    for (const [position, entries] of Object.entries(this.#data.contacts)) {
      if (parseInt(position) > this.#data.sequence.length) {
        delete this.#data.contacts[position];
        console.warn(
          "Delete contacts data at position " +
            position +
            "; Position exceeds the sequence length."
        );
      }
      for (const [contact, _] of Object.entries(entries)) {
        if (parseInt(contact) > this.#data.sequence.length) {
          delete this.#data.contacts[position][contact];
          console.warn(
            "Delete contact at position " +
              position +
              " to position " +
              contact +
              "; Position exceeds the sequence length."
          );
        }
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
          this.#data.aminoacidCounts[aa][modificationName] = [
            1,
            modificationClass,
          ];
        } else {
          this.#data.aminoacidCounts[aa][modificationName][0] += 1;
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
        return 1;
      } else if (_modificationsMassShift[m1] > _modificationsMassShift[m2]) {
        return -1;
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
      "Non-standard residue": "Modification",
      "Modified residue": "Modification",
      Lipidation: "Modification",
      Glycosylation: "Modification",
      "Disulfide bond": "Modification",
      "Cross-link": "Modification",
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

  /**
   * Internal method for updating the chart option.
   *
   * @param {Boolean} replace Passed to the setOption method of the Chart instance.
   */
  #updateOption(replace) {
    if (this.#option != undefined)
      this.chart.setOption(this.#option, { notMerge: replace });
  }

  /**
   * Highlights a contact in the structure view.
   *
   * @param {Number} x Position index of the first residue.
   * @param {Number} y Position index of the second residue.
   */
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

  /**
   * Checks whether listB contains all elements of listA.
   *
   * @param {Array} listA List of elements to check for.
   * @param {Array} listB List of elements to check in.
   * @returns True, if listB contains all elements of listA.
   */
  #contains(listA, listB) {
    const isContained = (entryA) => listB.includes(entryA);
    return listA.every(isContained);
  }
}

/**
 * Internal class for handling the structure view.
 */
class StructureView {
  /**
   * GLViewer instance for the structure view.
   */
  glviewer = null;
  /**
   * DOM identifier of the element to bind the the structure view to.
   */
  domId = null;
  /**
   * Data object containing the (structure) contact information.
   */
  #contacts = null;
  /**
   * Highlighted positions in the structure view.
   */
  #highlightPositions;
  /**
   * Highlighted modifications in the structure view.
   */
  #highlightModifications;
  /**
   * Highlighted contacts in the structure view.
   */
  #highlightContacts;
  /**
   * Colors for the structure view.
   */
  #colors = {
    contact_unmodified: "#AAAAAA",
    contact_modified: "#000000",
    contact_highlight: "#DC5754",
  };

  /**
   * Class constructor.
   *
   * @param {String} id DOM id of the element to bind the component to.
   */
  constructor(id) {
    this.domId = id;
    this.glviewer = $3Dmol.createViewer($("#" + id), {
      backgroundColor: "#fbfbfb",
      antialias: true,
      cartoonQuality: 6,
    });
  }

  /**
   * Sets the structure view with the passed protein structure and contact information.
   *
   * @param {String} pdbString Protein structure in PDB format.
   * @param {Object} contacts Object containing the contact information of residues within the passsed structure.
   */
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

  /**
   * Applies the style to the structure view based on current settings.
   */
  style() {
    this.glviewer.removeAllLabels();
    this.glviewer.removeAllShapes();
    this.glviewer.setStyle(
      {},
      {
        cartoon: {
          color: "#D4D4D4",
          radius: 0.4,
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
          backgroundOpacity: 0.8,
          fontColor: "#f0f5f5",
          fontSize: 14,
        }
      );
      if (this.#contacts.hasOwnProperty(targetIndex)) {
        for (let [contactIndex, distance] of Object.entries(
          this.#contacts[targetIndex]
        )) {
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
              backgroundOpacity: 0.8,
              fontColor: "#f0f5f5",
              fontSize: 14,
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
          this.glviewer.addLabel(distance.toFixed(2) + " Å", {
            backgroundColor: "rgb(51, 51, 51)",
            backgroundOpacity: 0.8,
            fontColor: "#f0f5f5",
            fontSize: 14,
            position: {
              x: (targetCaAtom.x + contactCaAtom.x) / 2,
              y: (targetCaAtom.y + contactCaAtom.y) / 2,
              z: (targetCaAtom.z + contactCaAtom.z) / 2,
            },
          });
        }
      }
    }
    // Add style for position highlight, if set.
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

  /**
   * Highlights all residues in contact with the specified resude in the structure view.
   *
   * @param {Number} targetIndex Index of the residue to highlight contact information for.
   */
  highlightContacts(targetIndex) {
    this.#highlightContacts = targetIndex;
    this.style();
  }

  /**
   * Highlights all residues with the specified modifications in the structure view.
   *
   * @param {Array} indices The residue indices to highlight.
   * @param {Array} modifications The modifications to highlight.
   */
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
  $("#menu")[0].style.display = "flex";
  $("#panel")[0].style.display = "block";
  _overviewTable = new OverviewTable("panel-table-tabulator");
  _overviewTable.registerSelectionAction((data) => {
    if (data.length > 0) {
      $("#dashboard-display-button")[0].disabled = false;
    } else {
      $("#dashboard-display-button")[0].disabled = true;
    }
  });
  _overviewChart = new OverviewChart("panel-overview-chart");
  _dashboardChart = new DashboardChart(
    "panel-dashboard-chart",
    "panel-dashboard-structure"
  );
  checkSessionState();
}

/**
 * Checks the current session state and updates the UI accordingly.
 */
function checkSessionState() {
  axios.get(window.location.origin + "/session_state").then((response) => {
    if (response.status == 200 && response.data.has_data) {
      overviewTableInitialize(overviewChartInitialize);
      if (response.data.hasOwnProperty("protein_selected"))
        dashboardChartInitialize(
          response.data.protein_selected,
          undefined,
          undefined,
          false
        );
    }
  });
}

/**
 * Starts an example session by making a GET request to the server and initializing the overview table and chart.
 *
 * @param {String} fileIdentifier The file identifier of the example session to start.
 */
function startExampleSession(fileIdentifier) {
  displayNotification("Initializing example session.");
  axios
    .get(
      window.location.origin +
        "/example_session?fileIdentifier=" +
        fileIdentifier
    )
    .then((_) => {
      clearCharts();
      overviewTableInitialize(); // Init. table.
      overviewChartInitialize(); // Init.overview chart.
    })
    .catch((error) => {
      console.error(error);
      displayAlert(error.response.data);
    })
    .finally(() => {
      removeNotification();
    });
}

/**
 * Restarts an existing PTMVision session.
 */
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
    .post(window.location.origin + "/restart_session", request, {
      headers: {
        "Content-Type": "text",
      },
    })
    .then((_) => {
      clearCharts();
      overviewTableInitialize(); // Init. table.
      overviewChartInitialize(); // Init.overview chart.
    })
    .catch((error) => {
      console.error(error);
      removeNotification();
      displayAlert(error.response.data);
    })
    .finally(() => {
      removeNotification();
    });
}

/**
 * Sends the specified search enginge output data to the PTMVision backend to start a new session.
 */
async function startSession() {
  displayNotification("Transfer and process entered data.");
  request = {
    massShiftTolerance: 0.001,
    excludeClasses: [],
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
  request.massShiftTolerance = parseFloat($("#data-tolerance-form")[0].value);
  request.excludeClasses = Metro.getPlugin(
    "#data-excludecls-form",
    "select"
  ).val();
  axios
    .post(
      window.location.origin + "/process_search_engine_output",
      pako.deflate(JSON.stringify(request)),
      {
        headers: {
          "Content-Type": "application/octet-stream",
          "Content-Encoding": "zlib",
        },
      }
    )
    .then((_) => {
      clearCharts();
      overviewTableInitialize(overviewChartInitialize); // Init. table and chart.
    })
    .catch((error) => {
      console.error(error);
      removeNotification();
      displayAlert(error.response.data);
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
 * Initializes the overview table by fetching the available proteins from the server.
 *
 * @param {Function} afterResponse Optional callback function to execute after the response was received.
 */
function overviewTableInitialize(afterResponse) {
  axios
    .get(window.location.origin + "/available_proteins")
    .then((response) => {
      _overviewTable.setData(response.data);
      modifications_data_string = "";
      _overviewTable.availableModifications.forEach((m) => {
        modifications_data_string +=
          `<option value="` + m + `">` + m + `</option>`;
      });
      Metro.getPlugin("#panel-table-filter-modification", "select").data(
        modifications_data_string
      );
      identifiers_data_string = "";
      _overviewTable.availableIdentifiers.forEach((i) => {
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
      displayAlert(error.response.data);
    });
}

/**
 * Sets the filters of the overview table.
 */
function overviewTableSetFilters() {
  _overviewTable.setFilters(
    Metro.getPlugin("#panel-table-filter-id", "select").val(),
    Metro.getPlugin("#panel-table-filter-modification", "select").val()
  );
}

/**
 * Initializes the overview chart by fetching the data from the server.
 *
 * @param {Function} afterResponse Optional callback function to execute after the response was received.
 */
function overviewChartInitialize(afterResponse) {
  displayNotification("Retrieving session data.");
  axios
    .get(window.location.origin + "/overview_data")
    .then((response) => {
      _overviewChart.fill(response.data);
      scroll_to("panel-overview");
      if (afterResponse != undefined) afterResponse();
    })
    .catch((error) => {
      console.error(error);
      removeNotification();
      displayAlert(error.response.data);
    })
    .finally(() => {
      removeNotification();
    });
}

/**
 * Downloads the image of the overview chart.
 */
function overviewChartDownloadImage() {
  const _ = document.createElement("a");
  document.body.appendChild(_);
  _.setAttribute("download", "overview.png");
  _.href = _overviewChart.getDataUrl();
  if (!_.href.endsWith("undefined")) _.click();
  _.remove();
}

/**
 * Restores the zoom of the overview chart.
 */
function overviewChartRestoreZoom() {
  _overviewChart.restoreZoom();
}

/**
 * Highlights the selected PTMs in the overview chart.
 */
function overviewChartHighlight() {
  let options = [];
  let names = _overviewChart.getDataNames();
  for (let i = 0; i < names.length; i++) {
    options.push(`<option value="` + i + `">` + names[i] + `</option>`);
  }
  Swal.fire({
    backdrop: false,
    confirmButtonColor: "#607196",
    width: "35em",
    padding: "1em",
    position: "center",
    html:
      `<h4><small>Select PTMs to highlight:</small></h4>
    <select id="tmpSelect" data-role="select" class="input-small" multiple>` +
      options.join("") +
      `</select>`,
  }).then((result) => {
    if (result.isConfirmed) {
      _overviewChart.highlight(
        Metro.getPlugin("#tmpSelect", "select")
          .val()
          .map((_) => parseInt(_))
      );
    }
  });
}

/**
 * Resorts the overview chart.
 */
function overviewChartSort() {
  _overviewChart.resort();
}

/**
 * Initialies the dashboard chart by fetching the data of one protein from the server.
 *
 * @param {String} uniprotPaValue Optional UniProt accession value to initialize the dashboard chart with. If undefined, the selected protein from the overview table is used.
 * @param {Number} cutoffValue Distance cutoff value to use for defining residue contacts.
 * @param {String} pdbTextValue Optional PDB text value to initialize the dashboard chart with.
 * @param {Boolean} scrollTo Optional flag to indicate whether to scroll to the dashboard panel after initialization. Default is true.
 */
function dashboardChartInitialize(
  uniprotPaValue,
  cutoffValue,
  pdbTextValue,
  scrollTo
) {
  displayNotification("Initializing dashboard.");
  let sel_protein_id = undefined;
  if (uniprotPaValue == undefined)
    sel_protein_id = _overviewTable.getSelection();
  else sel_protein_id = uniprotPaValue;
  if (cutoffValue == undefined) cutoffValue = 4.69;
  if (pdbTextValue == undefined) pdbTextValue = null;
  if (scrollTo == undefined) scrollTo = true;
  if (sel_protein_id == null) {
    removeNotification();
    displayAlert(
      `No protein was selected from panel <i class="fa-duotone fa-circle-3"></i>`
    );
    return;
  }
  request = {
    uniprot_pa: sel_protein_id,
    pdb_text: pdbTextValue,
    cutoff: cutoffValue,
  };
  axios
    .post(
      window.location.origin + "/protein_data",
      pako.deflate(JSON.stringify(request)),
      {
        headers: {
          "Content-Type": "application/octet-stream",
          "Content-Encoding": "zlib",
        },
      }
    )
    .then((response) => {
      if (response.status == 303) {
        displayAlert("Failed to fetch protein structure.");
      } else {
        // TODO: Prune protein data to remove modifications at positions that exceed the protein's length.
        let primaryAccession = response.data.annotation.primaryAccession;
        let proteinName = response.data.annotation.hasOwnProperty(
          "proteinDescription"
        )
          ? response.data.annotation.proteinDescription.hasOwnProperty(
              "recommendedName"
            )
            ? response.data.annotation.proteinDescription.recommendedName
                .fullName.value
            : Object.values(response.data.annotation.proteinDescription)[0][0]
                .fullName.value
          : "N/A";
        let organismName = response.data.annotation.hasOwnProperty("organism")
          ? "<em>" + response.data.annotation.organism.scientificName + "</em>"
          : "N/A";
        $("#panel-dashboard-selection").html(
          "Selected Protein: <u>" +
            [primaryAccession, proteinName, organismName].join(" | ") +
            "</u> (Click for UniProt Information)"
        );
        $("#panel-dashboard-selection").css("cursor", "pointer");
        $("#panel-dashboard-selection").on("click", () =>
          Swal.fire({
            html: response.data.annotation.comments
              .filter((_) => {
                return _.hasOwnProperty("texts");
              })
              .map((_) => {
                return (
                  "<code>" +
                  _.commentType +
                  "</code><p class='text-just' style='font-size: smaller;'>" +
                  _.texts[0].value +
                  "</p>"
                );
              })
              .join("</br>"),
            confirmButtonColor: "#607196",
            confirmButtonText: "Close Info",
            width: "35vw",
            heightAuto: false,
          })
        );
        // console.log(response.data); // Uncomment for local development.
        _dashboardChart.fill(response.data);
        _dashboardContent = "modifications";
        $("#panel-dashboard-title").html("Explore detail - Modifications view");
        if (scrollTo) scroll_to("panel-dashboard");
      }
    })
    .catch((error) => {
      console.error(error);
      removeNotification();
      displayAlert(error.response.data);
    })
    .finally(() => {
      removeNotification();
    });
}

/**
 * Downloads the image of the dashboard chart.
 */
function dashboardChartDownloadImage() {
  const _ = document.createElement("a");
  document.body.appendChild(_);
  _.setAttribute("download", "dashboard.png");
  _.href = _dashboardChart.getDataUrl();
  if (!_.href.endsWith("undefined")) _.click();
  _.remove();
}

/**
 * Restores the zoom of the dashboard chart.
 */
function dashboardChartRestoreZoom() {
  _dashboardChart.restoreZoom();
}

/**
 * Highlights the selected PTMs in the dashboard chart.
 */
function dashboardChartHighlight() {
  let options = [];
  let names = _dashboardChart.getDataNames();
  for (let i = 0; i < names.length; i++) {
    options.push(`<option value="` + i + `">` + names[i] + `</option>`);
  }
  Swal.fire({
    backdrop: false,
    confirmButtonColor: "#607196",
    width: "35em",
    padding: "1em",
    position: "center",
    html:
      `<h4><small>Select PTMs to highlight:</small></h4>
    <select id="tmpSelect" data-role="select" class="input-small" multiple>` +
      options.join("") +
      `</select>`,
  }).then((result) => {
    if (result.isConfirmed) {
      _dashboardChart.highlight(Metro.getPlugin("#tmpSelect", "select").val());
    }
  });
}

/**
 * Switches the content of the dashboard chart between modifications and structure view.
 */
function dashboardChartSwitch() {
  if (_dashboardContent == "modifications") {
    _dashboardContent = "structure";
    $("#panel-dashboard-title").html("Explore detail - Structure view");
  } else if (_dashboardContent == "structure") {
    _dashboardContent = "modifications";
    $("#panel-dashboard-title").html("Explore detail - Modifications view");
  }
  _dashboardChart.switchContent();
}

/**
 * Clears the overview and dashboard charts.
 */
function clearCharts() {
  _overviewChart.clear();
  _dashboardChart.clear();
}
