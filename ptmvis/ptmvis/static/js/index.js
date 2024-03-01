var __url = "http://127.0.0.1:5000/";
var __proteinsOverviewTable = null;
var __overviewChart = null;
var __overviewOption = null;
var __overviewSizeObserver = null;
var __structureView = null;
var __dashboardChart = null;
var __dashboardOption = null;
var __dashboardSizeObserver = null;
var __data = null;
var __modifications = [];

const COLORS = {
  M: {},
  i: 0,
  C: [
    "#d83034",
    "#f9e858",
    "#ff9d3a",
    "#4ecb8d",
    "#bdd373",
    "#7e4794",
    "#c701ff",
    "#ff73b6",
    "#008dff",
    "#003a7d",
  ],
};

const DASHBOARD_COMPONENTS = {
  positionCountsIndex: 0,
  ptmCountsIndex: 1,
  annotationsIndex: 2,
  presenceMapIndex: 3,
  contactMapIndex: 3,
  aaCountsIndex: 4,
  positions: null,
  modifications: null,
  positionCountsSeries: null,
  ptmCountsSeries: null,
  aaCountsSeries: null,
  presenceMapSeries: null,
  contactMapSeries: null,
  annotationSeries: null,
  annotationsLabels: null,
};

const AMINO_ACIDS = [
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

const UPLOAD_PDB_HTML = `
    <p>You can provide a structure in .pdb format to continue.</p>
    <br>
    <input id="optional-pdb-input" type="file" data-role="file" data-mode="drop">`;

/**
 * Inizializes all client side elements of the PTMVision application.
 */
function init() {
  __url = API_PARAMETERS["URL"];
  $("#menu")[0].style.display = "flex";
  $("#panel")[0].style.display = "block";
  __proteinsOverviewTable = new Tabulator("#panel-table-tabulator", {
    selectable: 1,
    columns: [
      {
        title: "ID",
        field: "id",
        sorter: "string",
        width: "10%",
      },
      {
        title: "Name",
        field: "name",
        sorter: "string",
        width: "40%",
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
        visible: false,
      },
    ],
  });
  __proteinsOverviewTable.on(
    "rowSelectionChanged",
    function (data, rows, selected, deselected) {
      if (data.length > 0) {
        $("#dashboard-display-button")[0].disabled = false;
      } else {
        $("#dashboard-display-button")[0].disabled = true;
      }
    }
  );
  __overviewChart = echarts.init($("#panel-overview-chart")[0], {
    devicePixelRatio: 2,
    renderer: "canvas",
    width: "auto",
    height: "auto",
  });
  __overviewSizeObserver = new ResizeObserver((entries) => {
    __overviewChart.resize({
      width: entries[0].width,
      height: entries[0].height,
    });
  });
  __overviewSizeObserver.observe($("#panel-overview-chart")[0]);
  __structureView = $3Dmol.createViewer($("#panel-dashboard-structure"), {
    backgroundColor: "#FAFAFC",
    antialias: true,
    cartoonQuality: 6,
  });
  hideStructure();
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
    console.log(response.data);
    const D = new Date();
    downloadBlob(
      response.data,
      "ptmvis-" +
        [D.getFullYear(), D.getMonth() + 1, Date.now()].join("-") +
        ".zlib"
    );
  });
}

function startExampleSession() {
  displayNotification("Initializing example session.");
  axios
    .get(__url + "/example_session")
    .then((_) => {
      initializeOverviewTable(); // Init. table.
      initializeOverviewChart(); // Init.overview chart.
      togglePanel("panel-inputs");
    })
    .catch((error) => {
      console.error(error);
      handleError(error.message);
    })
    .finally(() => {
      removeNotification();
    });
}

async function restartSession() {
  displayNotification("Re-initializing session.");
  request = null;
  if ($("#session-input-form")[0].files.length == 0) {
    handleError("No session data was supplied.");
    $("body").css("cursor", "auto");
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
      initializeOverviewTable(); // Init. table.
      initializeOverviewChart(); // Init.overview chart.
      togglePanel("panel-inputs");
    })
    .catch((error) => {
      console.error(error);
      handleError(error.message);
    })
    .finally(() => {
      removeNotification();
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
 * Displays a toast element with a custom error message to the user.
 *
 * @param {String} msg The error message to display.
 */
function handleError(msg) {
  Metro.toast.create("Fatal Error: " + msg, null, 6000, "error-toast");
}

/**
 * Redirects the browser to the GitHub repository of this project.
 */
function redirectSource() {
  window.open(
    "https://github.com/Integrative-Transcriptomics/PTMVision",
    "_blank"
  );
}

/**
 * Redirects the browser to the about page of this project.
 */
function redirectAbout() {
  window.open(__url + "/about", "_blank");
}

function toggleElement(id) {
  if ($("#" + id).is(":visible")) {
    $("#" + id).hide();
  } else {
    $("#" + id).show();
  }
}

/**
 * Toggles the visibility of the PTMVision structure panel.
 */
function toggleStructure() {
  let id = $("#panel-dashboard-structure")[0].parentNode.parentNode.id;
  if ($("#" + id).is(":visible")) {
    hideStructure();
  } else {
    showStructure();
  }
}

function togglePanel(id) {
  $("#" + id).css("height", "4vh");
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

function showStructure() {
  let id = $("#panel-dashboard-structure")[0].parentNode.parentNode.id;
  $("#panel-dashboard-toggle-structure-button").html("Hide Protein Structure");
  $("#" + id).show();
}

function hideStructure() {
  let id = $("#panel-dashboard-structure")[0].parentNode.parentNode.id;
  $("#panel-dashboard-toggle-structure-button").html("Show Protein Structure");
  $("#" + id).hide();
}

function setTableFilters(_) {
  const id_values = Metro.getPlugin("#panel-table-filter-id", "select").val();
  const mod_values = Metro.getPlugin(
    "#panel-table-filter-modification",
    "select"
  ).val();
  __proteinsOverviewTable.clearFilter(true);
  const filters = [];
  for (let mod_value of mod_values) {
    filters.push({
      field: "modifications",
      type: "like",
      value: mod_value,
    });
  }
  for (let id_value of id_values) {
    let id_value_fields = id_value.split("$");
    filters.push({
      field: id_value_fields[0] == "id" ? "id" : "name",
      type: "like",
      value: id_value_fields[1],
    });
  }
  __proteinsOverviewTable.setFilter(filters);
}

/**
 * Sends the specified search enginge output data to the PTMVision backend and loads the results in the overview table.
 */
async function uploadData() {
  displayNotification("Transfer and process entered data.");
  request = {
    contentType: null,
    content: null,
  };
  if ($("#data-input-form")[0].files.length == 0) {
    handleError("No search engine output data was supplied.");
    $("body").css("cursor", "auto");
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
      initializeOverviewTable(); // Init. table.
      initializeOverviewChart(); // Init.overview chart.
      togglePanel("panel-inputs");
    })
    .catch((error) => {
      console.error(error);
      handleError(error.message);
    })
    .finally(() => {
      removeNotification();
    });
}

function initializeOverviewTable() {
  axios
    .get(__url + "/get_available_proteins")
    .then((response) => {
      __proteinsOverviewTable.setData(response.data);
      __modifications = new Set();
      identifiers = new Set();
      for (let entry of response.data) {
        entry.modifications.split("$").forEach((m) => __modifications.add(m));
        identifiers.add(entry.id + "$" + entry.name);
      }
      __modifications = [...__modifications];
      modifications_data_string = "";
      __modifications.forEach((m) => {
        modifications_data_string +=
          `<option value="` + m + `">` + m + `</option>`;
      });
      Metro.getPlugin("#panel-table-filter-modification", "select").data(
        modifications_data_string
      );
      identifiers = [...identifiers];
      identifiers_data_string = "";
      identifiers.forEach((i) => {
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
    })
    .catch((error) => {
      throw error;
    });
}

function initializeOverviewChart() {
  axios
    .get(__url + "/get_overview_data")
    .then((response) => {
      __overviewOption = getOverviewOption(response.data);
      __overviewChart.setOption(__overviewOption);
      $("#panel-overview h6 button").attr("disabled", false);
      window.scrollTo({
        top: $("#panel-overview").get(0).offsetTop + 40,
        behavior: "smooth",
      });
    })
    .catch((error) => {
      throw error;
    });
}

function initializeDashboard(cutoff_value, pdb_text_value) {
  displayNotification("Initializing dashboard.");
  if (cutoff_value == null) {
    cutoff_value = 4.69;
  }
  if (pdb_text_value == null) {
    pdb_text_value = null;
  }
  request = {
    uniprot_id: __proteinsOverviewTable.getSelectedData()[0].id,
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
      if (response.data.status == "Failed: No protein structure available.") {
        _requestPdb();
      } else {
        __data = response.data.content;
        __structureView = getStructureView(
          __data.structure,
          $("#panel-dashboard-structure")
        );
        // Initialize dashboard components.
        [__dashboardChart, __dashboardSizeObserver] = constructChartInstance(
          $("#panel-dashboard-chart")[0],
          {}
        );
        [
          DASHBOARD_COMPONENTS.positionCountsSeries,
          DASHBOARD_COMPONENTS.positions,
          DASHBOARD_COMPONENTS.ptmCountsSeries,
          DASHBOARD_COMPONENTS.modifications,
          DASHBOARD_COMPONENTS.aaCountsSeries,
        ] = getDashboardCountsComponents(
          DASHBOARD_COMPONENTS.positionCountsIndex,
          DASHBOARD_COMPONENTS.ptmCountsIndex,
          DASHBOARD_COMPONENTS.aaCountsIndex
        );
        [
          DASHBOARD_COMPONENTS.presenceMapSeries,
          DASHBOARD_COMPONENTS.annotationSeries,
          DASHBOARD_COMPONENTS.annotationsLabels,
        ] = getPresenceMapComponents(
          DASHBOARD_COMPONENTS.presenceMapIndex,
          DASHBOARD_COMPONENTS.annotationsIndex
        );
        DASHBOARD_COMPONENTS.contactMapSeries = getContactMapComponents(
          DASHBOARD_COMPONENTS.contactMapIndex
        );
        setDashboardPresenceMap();
        $("#panel-dashboard h6 button").attr("disabled", false);
        structureViewSetDefaultStyle();
        window.scrollTo({
          top: $("#panel-dashboard").get(0).offsetTop + 40,
          behavior: "smooth",
        });
      }
      togglePanel("panel-overview");
      togglePanel("panel-table");
    })
    .catch((error) => {
      console.error(error);
      handleError(error.message);
    })
    .finally(() => {
      removeNotification();
    });
}

function setDashboardPresenceMap() {
  $("#panel-dashboard").attr("map_content", "presence");
  hideStructure();
  const axisProperties = {
    nameGap: 35,
    nameLocation: "center",
    nameTextStyle: {
      fontWeight: "bold",
      fontSize: 11,
    },
    axisTick: {
      alignWithLabel: true,
    },
    axisLabel: {
      fontWeight: "lighter",
      fontSize: 10,
    },
  };
  const titleProperties = {
    textStyle: {
      fontSize: 12,
    },
    borderRadius: 4,
    backgroundColor: "#f2f3f2",
  };
  var option = {
    title: [
      {
        text: "PTM Counts per Position",
        ...titleProperties,
        top: "1%",
        left: "3%",
      },
      {
        text: "PTM Occurrence per Class and Amino Acid",
        ...titleProperties,
        top: "29%",
        left: "3%",
      },
      {
        text: "Modifications Presence Map and Position Annotation",
        ...titleProperties,
        top: "2%",
        left: "59%",
      },
    ],
    grid: [
      {
        id: "grid_position_counts",
        top: "6%",
        left: "4%",
        height: "20%",
        width: "46%",
        containLabel: false,
        zlevel: 0,
      },
      {
        id: "grid_ptms_counts",
        top: "34%",
        left: "4%",
        height: "20%",
        width: "46%",
        containLabel: false,
        zlevel: 0,
      },
      {
        id: "grid_annotations",
        top: "5%",
        right: "2%",
        height: "19%",
        width: "39%",
        containLabel: false,
        zlevel: 1,
      },
      {
        id: "grid_presence_map",
        top: "25%",
        right: "2%",
        height: "69%",
        width: "39%",
        containLabel: false,
        zlevel: 1,
      },
      {
        id: "grid_aa_counts",
        top: "55%",
        left: "4%",
        height: "34%",
        width: "46%",
        containLabel: false,
        zlevel: 0,
      },
    ],
    xAxis: [
      {
        data: [],
        name: "Protein Position",
        ...axisProperties,
        gridIndex: DASHBOARD_COMPONENTS.positionCountsIndex,
      },
      {
        data: [],
        axisLabel: {
          show: false,
        },
        axisTick: {
          show: false,
        },
        gridIndex: DASHBOARD_COMPONENTS.ptmCountsIndex,
      },
      {
        data: [],
        axisLabel: {
          show: false,
        },
        axisTick: {
          show: false,
        },
        gridIndex: DASHBOARD_COMPONENTS.annotationsIndex,
      },
      {
        data: [],
        name: "Protein Position",
        ...axisProperties,
        gridIndex: DASHBOARD_COMPONENTS.presenceMapIndex,
      },
      {
        data: [],
        name: "Modification",
        ...axisProperties,
        gridIndex: DASHBOARD_COMPONENTS.aaCountsIndex,
      },
    ],
    yAxis: [
      {
        type: "value",
        name: "No. PTMs",
        ...axisProperties,
        gridIndex: DASHBOARD_COMPONENTS.positionCountsIndex,
      },
      {
        type: "value",
        name: "Count",
        ...axisProperties,
        gridIndex: DASHBOARD_COMPONENTS.ptmCountsIndex,
      },
      {
        data: [],
        inverse: false,
        name: "Annotations",
        ...axisProperties,
        nameGap: 140,
        gridIndex: DASHBOARD_COMPONENTS.annotationsIndex,
      },
      {
        data: [],
        name: "Modification",
        ...axisProperties,
        nameGap: 140,
        interval: 0,
        gridIndex: DASHBOARD_COMPONENTS.presenceMapIndex,
      },
      {
        data: AMINO_ACIDS,
        name: "Amino Acid",
        ...axisProperties,
        gridIndex: DASHBOARD_COMPONENTS.aaCountsIndex,
      },
    ],
    dataZoom: [
      {
        type: "inside",
        xAxisIndex: [
          DASHBOARD_COMPONENTS.positionCountsIndex,
          DASHBOARD_COMPONENTS.annotationsIndex,
          DASHBOARD_COMPONENTS.presenceMapIndex,
        ],
        throttle: 0,
      },
      {
        type: "inside",
        xAxisIndex: [
          DASHBOARD_COMPONENTS.ptmCountsIndex,
          DASHBOARD_COMPONENTS.aaCountsIndex,
        ],
        yAxisIndex: [DASHBOARD_COMPONENTS.presenceMapIndex],
        throttle: 0,
      },
      {
        type: "inside",
        yAxisIndex: [DASHBOARD_COMPONENTS.aaCountsIndex],
        throttle: 0,
      },
    ],
    legend: [
      {
        id: "legend",
        bottom: "bottom",
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
        selector: true,
      },
    ],
    axisPointer: [
      {
        show: true,
        triggerTooltip: false,
        link: [
          {
            xAxisIndex: [
              DASHBOARD_COMPONENTS.presenceMapIndex,
              DASHBOARD_COMPONENTS.annotationsIndex,
            ],
          },
          {
            xAxisIndex: [
              DASHBOARD_COMPONENTS.ptmCountsIndex,
              DASHBOARD_COMPONENTS.aaCountsIndex,
            ],
          },
        ],
        label: { backgroundColor: "rgba(51, 51, 51, 0.7)" },
      },
    ],
    tooltip: [
      {
        trigger: "axis",
        backgroundColor: "rgba(51, 51, 51, 0.7)",
        borderColor: "transparent",
        textStyle: {
          fontWeight: "lighter",
          fontSize: 10,
          color: "#f0f5f5",
        },
      },
    ],
    visualMap: [
      {
        bottom: "5%",
        left: "10%",
        color: ["#111111", "#d1d1d1"],
        orient: "horizontal",
        itemWidth: 10,
        itemHeight: 100,
        precision: 0,
        textStyle: {
          fontSize: 10,
          fontWeight: "lighter",
        },
      },
    ],
  };
  option.xAxis[DASHBOARD_COMPONENTS.positionCountsIndex].data =
    DASHBOARD_COMPONENTS.positions;
  option.xAxis[DASHBOARD_COMPONENTS.presenceMapIndex].data =
    DASHBOARD_COMPONENTS.positions;
  option.xAxis[DASHBOARD_COMPONENTS.annotationsIndex].data =
    DASHBOARD_COMPONENTS.positions;
  option.xAxis[DASHBOARD_COMPONENTS.aaCountsIndex].data =
    DASHBOARD_COMPONENTS.modifications;
  option.xAxis[DASHBOARD_COMPONENTS.ptmCountsIndex].data =
    DASHBOARD_COMPONENTS.modifications;
  option.yAxis[DASHBOARD_COMPONENTS.presenceMapIndex].data =
    DASHBOARD_COMPONENTS.modifications;
  option.yAxis[DASHBOARD_COMPONENTS.annotationsIndex].data =
    DASHBOARD_COMPONENTS.annotationsLabels;
  option.series = [
    ...DASHBOARD_COMPONENTS.positionCountsSeries,
    ...DASHBOARD_COMPONENTS.ptmCountsSeries,
    ...DASHBOARD_COMPONENTS.presenceMapSeries,
    ...DASHBOARD_COMPONENTS.annotationSeries,
    ...DASHBOARD_COMPONENTS.aaCountsSeries,
  ];
  option.tooltip[0].formatter = (params, ticket, callback) => {
    let content = ``;
    let axisIndex = params[0].axisIndex;
    if (
      axisIndex == DASHBOARD_COMPONENTS.positionCountsIndex ||
      axisIndex == DASHBOARD_COMPONENTS.presenceMapIndex ||
      axisIndex == DASHBOARD_COMPONENTS.annotationsIndex
    ) {
      content +=
        "<b>Position " +
        parseInt(params[0].name) +
        "</b> (" +
        __data.sequence[params[0].name - 1] +
        ")</br>";
    } else if (
      axisIndex == DASHBOARD_COMPONENTS.ptmCountsIndex ||
      axisIndex == DASHBOARD_COMPONENTS.aaCountsIndex
    ) {
      content += "<b>Modification " + params[0].name + "</b></br>";
    }
    for (let entry of params) {
      if (
        entry.axisIndex == DASHBOARD_COMPONENTS.positionCountsIndex ||
        entry.axisIndex == DASHBOARD_COMPONENTS.ptmCountsIndex
      ) {
        content += entry.seriesName + " " + entry.data[1] + "</br>";
      } else if (entry.axisIndex == DASHBOARD_COMPONENTS.aaCountsIndex) {
        continue;
      } else if (
        entry.axisIndex == DASHBOARD_COMPONENTS.presenceMapIndex ||
        entry.axisIndex == DASHBOARD_COMPONENTS.annotationsIndex
      ) {
        if (entry.seriesId.startsWith("annotation@")) {
          content += "&#8226; " + entry.seriesName + "</br>";
          for (let _ of entry.value[3]) {
            content +=
              _[0] +
              " (Pos. " +
              _[1] +
              " to " +
              _[2] +
              ") <small>" +
              _[3] +
              "</small></br>";
          }
        }
      }
    }
    return content;
  };
  option.visualMap[0].seriesIndex = option.series.length - 1;
  option.visualMap[0].min = 1;
  option.visualMap[0].max = Math.max(
    ...option.series[option.series.length - 1].data.map((_) => _[2])
  );
  option.visualMap[0].text = [
    option.visualMap[0].max,
    "Modification Occurrence " + option.visualMap[0].min,
  ];
  option.legend[0].data = Array.from(
    new Set(
      option.series
        .filter(
          (_) => !_.id.startsWith("annotation") && !_.id.startsWith("aa_counts")
        )
        .map((_) => {
          return { name: _.name };
        })
    )
  );
  __dashboardChart.setOption(option, { notMerge: true });
  __dashboardChart.on("click", (params) => {
    if (
      params.seriesId.startsWith("position_counts") ||
      params.seriesId.startsWith("presence_map") ||
      params.seriesId.startsWith("annotation")
    ) {
      structureViewHighlightContacts(params.data[0] + 1, 4.69);
    }
  });
}

function setDashboardContactMap() {
  $("#panel-dashboard").attr("map_content", "contact");
  showStructure();
  const axisProperties = {
    nameGap: 35,
    nameLocation: "center",
    nameTextStyle: {
      fontWeight: "bold",
      fontSize: 11,
    },
    axisTick: {
      alignWithLabel: true,
    },
    axisLabel: {
      fontWeight: "lighter",
      fontSize: 10,
    },
  };
  const titleProperties = {
    textStyle: {
      fontSize: 12,
    },
    borderRadius: 4,
    backgroundColor: "#f2f3f2",
  };
  var option = {
    title: [
      {
        text: "PTM Counts per Position",
        ...titleProperties,
        top: "1%",
        left: "3%",
      },
      {
        text: "PTM Occurrence per Class and Amino Acid",
        ...titleProperties,
        top: "29%",
        left: "3%",
      },
      {
        text: "Contact Map and Position Annotation",
        ...titleProperties,
        top: "2%",
        left: "59%",
      },
      {
        id: "title_contact_detail",
        ...titleProperties,
        text: "Contact Modifications",
        top: "60%",
        left: "28%",
      },
    ],
    grid: [
      {
        id: "grid_position_counts",
        top: "6%",
        left: "4%",
        height: "20%",
        width: "46%",
        containLabel: false,
        zlevel: 0,
      },
      {
        id: "grid_ptms_counts",
        top: "34%",
        left: "4%",
        height: "20%",
        width: "46%",
        containLabel: false,
        zlevel: 0,
      },
      {
        id: "grid_annotations",
        top: "5%",
        right: "2%",
        height: "19%",
        width: "39%",
        containLabel: false,
        zlevel: 1,
      },
      {
        id: "grid_contact_map",
        top: "25%",
        right: "2%",
        height: "69%",
        width: "39%",
        containLabel: false,
        zlevel: 1,
      },
      {
        id: "grid_contact_detail",
        top: "63%",
        left: "28%",
        height: "33%",
        width: "22%",
        containLabel: false,
        show: true,
        backgroundColor: "#fdfefd",
        zlevel: 0,
      },
    ],
    xAxis: [
      {
        data: [],
        name: "Protein Position",
        ...axisProperties,
        gridIndex: DASHBOARD_COMPONENTS.positionCountsIndex,
      },
      {
        data: [],
        name: "Modification",
        ...axisProperties,
        gridIndex: DASHBOARD_COMPONENTS.ptmCountsIndex,
      },
      {
        data: [],
        axisLabel: {
          show: false,
        },
        axisTick: {
          show: false,
        },
        gridIndex: DASHBOARD_COMPONENTS.annotationsIndex,
      },
      {
        data: [],
        name: "Protein Position",
        ...axisProperties,
        gridIndex: DASHBOARD_COMPONENTS.contactMapIndex,
      },
    ],
    yAxis: [
      {
        type: "value",
        name: "No. PTMs",
        ...axisProperties,
        gridIndex: DASHBOARD_COMPONENTS.positionCountsIndex,
      },
      {
        type: "value",
        name: "Count",
        ...axisProperties,
        gridIndex: DASHBOARD_COMPONENTS.ptmCountsIndex,
      },
      {
        data: [],
        inverse: false,
        name: "Annotations",
        ...axisProperties,
        nameGap: 140,
        gridIndex: DASHBOARD_COMPONENTS.annotationsIndex,
      },
      {
        data: [],
        name: "Protein Position",
        ...axisProperties,
        nameGap: 40,
        gridIndex: DASHBOARD_COMPONENTS.contactMapIndex,
      },
    ],
    dataZoom: [
      {
        type: "inside",
        xAxisIndex: [
          DASHBOARD_COMPONENTS.positionCountsIndex,
          DASHBOARD_COMPONENTS.annotationsIndex,
          DASHBOARD_COMPONENTS.contactMapIndex,
        ],
        throttle: 0,
      },
      {
        type: "inside",
        yAxisIndex: [DASHBOARD_COMPONENTS.contactMapIndex],
        throttle: 0,
      },
      {
        type: "inside",
        xAxisIndex: [DASHBOARD_COMPONENTS.ptmCountsIndex],
        throttle: 0,
      },
    ],
    legend: [
      {
        id: "legend",
        bottom: "bottom",
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
        selector: true,
      },
    ],
    axisPointer: [
      {
        show: true,
        triggerTooltip: false,
        link: [
          {
            xAxisIndex: [
              DASHBOARD_COMPONENTS.contactMapIndex,
              DASHBOARD_COMPONENTS.annotationsIndex,
            ],
          },
        ],
        label: { backgroundColor: "rgba(51, 51, 51, 0.7)" },
      },
    ],
    tooltip: [
      {
        trigger: "axis",
        backgroundColor: "rgba(51, 51, 51, 0.7)",
        borderColor: "transparent",
        appendToBody: true,
        className: "panel-dashboard-tooltip",
        textStyle: {
          fontWeight: "lighter",
          fontSize: 10,
          color: "#f0f5f5",
        },
      },
    ],
  };
  option.xAxis[DASHBOARD_COMPONENTS.positionCountsIndex].data =
    DASHBOARD_COMPONENTS.positions;
  option.xAxis[DASHBOARD_COMPONENTS.ptmCountsIndex].data =
    DASHBOARD_COMPONENTS.modifications;
  option.xAxis[DASHBOARD_COMPONENTS.annotationsIndex].data =
    DASHBOARD_COMPONENTS.positions;
  option.xAxis[DASHBOARD_COMPONENTS.contactMapIndex].data =
    DASHBOARD_COMPONENTS.positions;
  option.yAxis[DASHBOARD_COMPONENTS.annotationsIndex].data =
    DASHBOARD_COMPONENTS.annotationsLabels;
  option.yAxis[DASHBOARD_COMPONENTS.contactMapIndex].data =
    DASHBOARD_COMPONENTS.positions;
  option.series = [
    ...DASHBOARD_COMPONENTS.positionCountsSeries,
    ...DASHBOARD_COMPONENTS.ptmCountsSeries,
    ...DASHBOARD_COMPONENTS.annotationSeries,
    ...DASHBOARD_COMPONENTS.contactMapSeries,
  ];
  option.tooltip[0].formatter = (params, ticket, callback) => {
    let content = ``;
    let axisIndex = params[0].axisIndex;
    if (
      axisIndex == DASHBOARD_COMPONENTS.positionCountsIndex ||
      axisIndex == DASHBOARD_COMPONENTS.contactMapIndex ||
      axisIndex == DASHBOARD_COMPONENTS.annotationsIndex
    ) {
      content +=
        "<b>Position " +
        parseInt(params[0].name) +
        "</b> (" +
        __data.sequence[params[0].name - 1] +
        ")</br>";
    } else if (axisIndex == DASHBOARD_COMPONENTS.ptmCountsIndex) {
      content += "<b>Modification " + params[0].name + "</b></br>";
    }
    for (let entry of params) {
      if (
        entry.axisIndex == DASHBOARD_COMPONENTS.positionCountsIndex ||
        entry.axisIndex == DASHBOARD_COMPONENTS.ptmCountsIndex
      ) {
        content += entry.seriesName + " " + entry.data[1] + "</br>";
      } else if (
        entry.axisIndex == DASHBOARD_COMPONENTS.contactMapIndex ||
        entry.axisIndex == DASHBOARD_COMPONENTS.annotationsIndex
      ) {
        if (entry.seriesId.startsWith("annotation@")) {
          content += "&#8226; " + entry.seriesName + "</br>";
          for (let _ of entry.value[3]) {
            content +=
              _[0] +
              " (Pos. " +
              _[1] +
              " to " +
              _[2] +
              ") <small>" +
              _[3] +
              "</small></br>";
          }
        }
      }
    }
    return content;
  };
  option.legend[0].data = Array.from(
    new Set(
      option.series
        .filter(
          (_) =>
            !_.id.startsWith("annotation") && !_.id.startsWith("contact_map")
        )
        .map((_) => {
          return { name: _.name };
        })
    )
  );
  __dashboardChart.setOption(option, { notMerge: true });
  __dashboardChart.on("click", (params) => {
    if (
      params.seriesId.startsWith("position_counts") ||
      params.seriesId.startsWith("contact_map") ||
      params.seriesId.startsWith("annotation")
    ) {
      structureViewHighlightContacts(params.data[0] + 1, 4.69);
    }
    if (
      params.seriesId.startsWith("contact_map") &&
      params.seriesId == "contact_map@upper"
    ) {
      let option = getContactDetail(params.data[0], params.data[1]);
      __dashboardChart.setOption({ title: [option[0]], series: [option[1]] });
    } else {
      __dashboardChart.setOption({
        title: [
          {
            id: "title_contact_detail",
            ...titleProperties,
            text: "Contact Modifications",
            top: "60%",
            left: "28%",
          },
        ],
        series: [
          {
            id: "detail@contact",
            data: [],
            links: [],
          },
        ],
      });
    }
  });
}

function switchDashboardMap() {
  if ($("#panel-dashboard").attr("map_content") == "presence") {
    $("#panel-dashboard-switch-map-button").html("Show Presence Map");
    setDashboardContactMap();
  } else if ($("#panel-dashboard").attr("map_content") == "contact") {
    $("#panel-dashboard-switch-map-button").html("Show Contact Map");
    setDashboardPresenceMap();
  }
}

function _requestPdb() {
  Swal.fire({
    title: "Failed to fetch Protein Structure!",
    html: UPLOAD_PDB_HTML,
    width: "1000%",
    position: "bottom",
    showCloseButton: false,
    showCancelButton: true,
    cancelButtonColor: "#dc5754",
    showConfirmButton: true,
    confirmButtonColor: "#62a8ac",
    confirmButtonText: "Continue",
    grow: true,
    heightAuto: true,
    color: "#333333",
    background: "#fafafcd9",
    backdrop: `
      rgba(239, 240, 248, 0.1)
      left top
      no-repeat
    `,
  }).then((result) => {
    if (result.isConfirmed) {
      initializeDashboard(null, readFile($("#optional-pdb-input")[0].files[0]));
    }
  });
}

function constructChartInstance(html_container, echart_option) {
  var chart = echarts.init(html_container, {
    devicePixelRatio: 4,
    renderer: "canvas",
    width: "auto",
    height: "auto",
  });
  var observer = new ResizeObserver((entries) => {
    chart.resize({
      width: entries[0].width,
      height: entries[0].height,
    });
  });
  observer.observe(html_container);
  chart.setOption(echart_option, true);
  return [chart, observer];
}

function getStructureView(pdb_text, html_container) {
  let view = $3Dmol.createViewer(html_container, {
    backgroundColor: "#FAFAFC",
    antialias: true,
    cartoonQuality: 6,
  });
  view.clear();
  view.addModel(pdb_text, "pdb");
  view.zoomTo();
  view.addSurface(
    "SAS",
    {
      color: "#d4d4d4",
      opacity: 0.4,
    },
    {}
  );
  structureViewSetDefaultStyle();
  view.render();
  return view;
}

function structureViewHighlightContacts(position, cutoff) {
  __structureView.removeAllLabels();
  __structureView.removeAllShapes();
  structureViewSetDefaultStyle();
  let selectedResidue = position;
  let selectedResidueCA = __structureView.getAtomsFromSel({
    resi: [selectedResidue],
    atom: "CA",
  })[0];
  __structureView.addStyle(
    {
      resi: [selectedResidue],
    },
    {
      stick: {
        color: "#dc5754",
        radius: 0.7,
      },
    }
  );
  __structureView.addResLabels(
    {
      resi: [selectedResidue],
    },
    {
      backgroundColor: "rgb(51, 51, 51)",
      backgroundOpacity: 0.7,
      fontColor: "#f0f5f5",
      fontSize: 11,
    }
  );
  if (__data.contacts.hasOwnProperty(selectedResidue - 1)) {
    for (let entry of __data.contacts[selectedResidue - 1]) {
      let contactResidue = entry[0] + 1;
      __structureView.addStyle(
        {
          resi: [contactResidue],
        },
        {
          stick: {
            color: "#dd9ac2",
            radius: 0.4,
          },
        }
      );
      __structureView.addResLabels(
        {
          resi: [contactResidue],
        },
        {
          backgroundColor: "rgb(51, 51, 51)",
          backgroundOpacity: 0.7,
          fontColor: "#f0f5f5",
          fontSize: 11,
        }
      );
      let contactEntryCA = __structureView.getAtomsFromSel({
        resi: [contactResidue],
        atom: "CA",
      })[0];
      __structureView.addLine({
        color: "#000000",
        hidden: false,
        dashed: false,
        start: {
          x: selectedResidueCA.x,
          y: selectedResidueCA.y,
          z: selectedResidueCA.z,
        },
        end: {
          x: contactEntryCA.x,
          y: contactEntryCA.y,
          z: contactEntryCA.z,
        },
      });
      __structureView.addLabel(entry[1].toFixed(2) + " Å", {
        backgroundColor: "rgb(51, 51, 51)",
        backgroundOpacity: 0.4,
        fontColor: "#f0f5f5",
        fontSize: 10,
        position: {
          x: (selectedResidueCA.x + contactEntryCA.x) / 2,
          y: (selectedResidueCA.y + contactEntryCA.y) / 2,
          z: (selectedResidueCA.z + contactEntryCA.z) / 2,
        },
      });
    }
  }
  __structureView.render();
}

function structureViewSetDefaultStyle(selection) {
  if (selection == undefined) {
    selection = {};
  }
  __structureView.setStyle(selection, {
    cartoon: {
      color: "#d4d4d4",
    },
  });
}

function getDashboardCountsComponents(
  axisIndexPositionCounts,
  axisIndexPTMCounts,
  axisIndexAACounts
) {
  let modifications = [];
  let aa_counts = {};
  let ptms_counts = {};
  let position_counts = {};
  let aa_counts_series = {
    id: "aa_counts",
    name: "amino_acid_counts",
    type: "heatmap",
    data: [],
    xAxisIndex: axisIndexAACounts,
    yAxisIndex: axisIndexAACounts,
  };
  let ptm_counts_series = {};
  let position_counts_series = {};
  for (const [position, position_data] of Object.entries(__data.positions)) {
    const aa = __data.sequence[position - 1];
    if (!aa_counts.hasOwnProperty(aa)) {
      aa_counts[aa] = {};
    }
    if (!position_counts.hasOwnProperty(position)) {
      position_counts[position] = {};
    }
    for (const modification of Object.values(position_data.modifications)) {
      const modification_name = modification.modification_unimod_name;
      const modification_classification =
        modification.modification_classification;
      modifications.push(modification_name);
      if (!aa_counts[aa].hasOwnProperty(modification_name)) {
        aa_counts[aa][modification_name] = 1;
      } else {
        aa_counts[aa][modification_name] += 1;
      }
      if (!ptms_counts.hasOwnProperty(modification_name)) {
        ptms_counts[modification_name] = {};
      }
      if (
        !ptms_counts[modification_name].hasOwnProperty(
          modification_classification
        )
      ) {
        ptms_counts[modification_name][modification_classification] = 1;
      } else {
        ptms_counts[modification_name][modification_classification] += 1;
      }
      if (
        !position_counts[position].hasOwnProperty(modification_classification)
      ) {
        position_counts[position][modification_classification] = 1;
      } else {
        position_counts[position][modification_classification] += 1;
      }
    }
  }
  __modifications = [...new Set(modifications)].sort();
  for (const [amino_acid, _] of Object.entries(aa_counts)) {
    for (const [modification_name, count] of Object.entries(_)) {
      aa_counts_series.data.push([
        __modifications.indexOf(modification_name),
        AMINO_ACIDS.indexOf(amino_acid),
        count,
      ]);
    }
  }
  for (const [modification_name, _] of Object.entries(ptms_counts)) {
    for (const [modification_class, count] of Object.entries(_)) {
      if (!ptm_counts_series.hasOwnProperty(modification_class)) {
        ptm_counts_series[modification_class] = {
          id: "ptm_counts@" + modification_class,
          name: modification_class,
          stack: "total2",
          type: "bar",
          emphasis: {
            focus: "series",
          },
          data: [],
          itemStyle: {
            color: getColor(modification_class),
          },
          xAxisIndex: axisIndexPTMCounts,
          yAxisIndex: axisIndexPTMCounts,
        };
      }
      ptm_counts_series[modification_class].data.push([
        __modifications.indexOf(modification_name),
        count,
      ]);
    }
  }
  for (const [position, _] of Object.entries(position_counts)) {
    for (const [modification_class, count] of Object.entries(_)) {
      if (!position_counts_series.hasOwnProperty(modification_class)) {
        position_counts_series[modification_class] = {
          id: "position_counts@" + modification_class,
          name: modification_class,
          type: "bar",
          stack: "total",
          emphasis: {
            focus: "series",
          },
          data: [],
          itemStyle: {
            color: getColor(modification_class),
          },
          xAxisIndex: axisIndexPositionCounts,
          yAxisIndex: axisIndexPositionCounts,
        };
      }
      position_counts_series[modification_class].data.push([
        position - 1,
        count,
      ]);
    }
  }
  return [
    Object.values(position_counts_series),
    [...Array(__data.sequence.length).keys()].map((x) => x + 1),
    Object.values(ptm_counts_series),
    __modifications,
    [aa_counts_series],
  ];
}

function getPresenceMapComponents(axisIndexPresenceMap, axisIndexAnnotations) {
  let presence_map_series = {};
  for (let [position, position_data] of Object.entries(__data.positions)) {
    for (let modification of Object.values(position_data.modifications)) {
      if (
        !presence_map_series.hasOwnProperty(
          modification.modification_classification
        )
      ) {
        presence_map_series[modification.modification_classification] = {
          id: "presence_map@" + modification.modification_classification,
          name: modification.modification_classification,
          type: "heatmap",
          data: [],
          itemStyle: {
            color: getColor(modification.modification_classification),
          },
          xAxisIndex: axisIndexPresenceMap,
          yAxisIndex: axisIndexPresenceMap,
        };
      }
      let i = parseInt(position) - 1;
      presence_map_series[modification.modification_classification].data.push([
        i,
        __modifications.indexOf(modification.modification_unimod_name),
        1,
        modification.modification_unimod_name,
        modification.modification_classification,
        __data.sequence[i],
      ]);
    }
  }
  let annotation_types = {
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
  let annotation_series = {
    "Molecule processing": {
      id: "annotation@molecule_processing",
      name: "Molecule processing",
      type: "heatmap",
      data: {},
      xAxisIndex: axisIndexAnnotations,
      yAxisIndex: axisIndexAnnotations,
      itemStyle: {
        color: "#262626",
      },
    },
    Region: {
      id: "annotation@region",
      name: "Region",
      type: "heatmap",
      data: {},
      xAxisIndex: axisIndexAnnotations,
      yAxisIndex: axisIndexAnnotations,
      itemStyle: {
        color: "#082a54",
      },
    },
    Site: {
      id: "annotation@site",
      name: "Site",
      type: "heatmap",
      data: {},
      xAxisIndex: axisIndexAnnotations,
      yAxisIndex: axisIndexAnnotations,
      itemStyle: {
        color: "#e02b35",
      },
    },
    "Amino acid modifications": {
      id: "annotation@amino_acid_modifications",
      name: "Amino acid modifications",
      type: "heatmap",
      data: {},
      xAxisIndex: axisIndexAnnotations,
      yAxisIndex: axisIndexAnnotations,
      itemStyle: {
        color: "#f0c571",
      },
    },
    "Natural variations": {
      id: "annotation@natural_variations",
      name: "Natural variations",
      type: "heatmap",
      data: {},
      xAxisIndex: axisIndexAnnotations,
      yAxisIndex: axisIndexAnnotations,
      itemStyle: {
        color: "#59a89c",
      },
    },
    "Experimental info": {
      id: "annotation@experimental_info",
      name: "Experimental info",
      type: "heatmap",
      data: {},
      xAxisIndex: axisIndexAnnotations,
      yAxisIndex: axisIndexAnnotations,
      itemStyle: {
        color: "#a559aa",
      },
    },
    "Secondary structure": {
      id: "annotation@secondary_structure",
      name: "Secondary structure",
      type: "heatmap",
      data: {},
      xAxisIndex: axisIndexAnnotations,
      yAxisIndex: axisIndexAnnotations,
    },
  };
  let annotation_secondary_structure_colors = {
    Helix: "#c31e23",
    Turn: "#c99b38",
    "Beta strand": "#3594cc",
  };
  let annotation_labels = [...Object.keys(annotation_series)];
  for (let annotation_entry of Object.values(__data.annotation.features)) {
    if (annotation_entry.type == "Chain") {
      continue;
    }
    for (
      let index = parseInt(annotation_entry.location.start.value) - 1;
      index < parseInt(annotation_entry.location.end.value);
      index++
    ) {
      let y_index = annotation_labels.indexOf(
        annotation_types[annotation_entry.type]
      );
      if (annotation_types[annotation_entry.type] == "Secondary structure") {
        annotation_series[annotation_types[annotation_entry.type]].data[
          index + "@" + y_index
        ] = {
          value: [
            index,
            y_index,
            1,
            [
              [
                annotation_entry.type,
                annotation_entry.location.start.value,
                annotation_entry.location.end.value,
                "",
              ],
            ],
          ],
          itemStyle: {
            color: annotation_secondary_structure_colors[annotation_entry.type],
          },
        };
      } else {
        if (
          annotation_series[
            annotation_types[annotation_entry.type]
          ].data.hasOwnProperty(index + "@" + y_index)
        ) {
          annotation_series[annotation_types[annotation_entry.type]].data[
            index + "@" + y_index
          ][3].push([
            annotation_entry.type,
            annotation_entry.location.start.value,
            annotation_entry.location.end.value,
            annotation_entry.description,
          ]);
        } else {
          annotation_series[annotation_types[annotation_entry.type]].data[
            index + "@" + y_index
          ] = [
            index,
            y_index,
            1,
            [
              [
                annotation_entry.type,
                annotation_entry.location.start.value,
                annotation_entry.location.end.value,
                annotation_entry.description,
              ],
            ],
          ];
        }
      }
    }
  }
  let annotation_series_list = [];
  for (const value of Object.values(annotation_series)) {
    value.data = [...Object.values(value.data)];
    annotation_series_list.push(value);
  }
  return [
    Object.values(presence_map_series),
    annotation_series_list,
    annotation_labels,
  ];
}

function getContactMapComponents(axisIndexContactMap) {
  let contact_map_series = [
    {
      id: "contact_map@upper",
      name: "contact_map_upper",
      type: "heatmap",
      data: [],
      xAxisIndex: axisIndexContactMap,
      yAxisIndex: axisIndexContactMap,
      itemStyle: {
        color: "#dc5754",
      },
    },
    {
      id: "contact_map@lower",
      type: "heatmap",
      name: "contact_map_lower",
      data: [],
      xAxisIndex: axisIndexContactMap,
      yAxisIndex: axisIndexContactMap,
      itemStyle: {
        color: "#333333",
      },
    },
  ];
  for (let [residue_index_i, contacts_data] of Object.entries(
    __data.contacts
  )) {
    let x = parseInt(residue_index_i);
    for (let contact_data of contacts_data) {
      let y = contact_data[0];
      if (x > y) {
        var xModifications = new Set();
        var yModifications = new Set();
        if (__data.positions.hasOwnProperty(x + 1)) {
          __data.positions[x + 1].modifications.forEach((_) =>
            xModifications.add(_.modification_unimod_name)
          );
        }
        if (__data.positions.hasOwnProperty(y + 1)) {
          __data.positions[y + 1].modifications.forEach((_) =>
            yModifications.add(_.modification_unimod_name)
          );
        }
        xModifications = Array.from(xModifications);
        yModifications = Array.from(yModifications);
        var modificationsCount = new Set([...xModifications, ...yModifications])
          .size;
        if (xModifications.length > 0 || yModifications.length > 0) {
          contact_map_series[0].data.push([
            x,
            y,
            contact_data[1],
            modificationsCount,
            xModifications,
            yModifications,
          ]);
        }
      } else {
        contact_map_series[1].data.push([x, y, contact_data[1]]);
      }
    }
  }
  return contact_map_series;
}

function getContactDetail(index_i, index_j) {
  let position_i = index_i + 1;
  let position_j = index_j + 1;
  let residue_i = __data.sequence[index_i];
  let residue_j = __data.sequence[index_j];
  let modifications = {};
  if (__data.positions.hasOwnProperty(position_i)) {
    for (let modification of __data.positions[position_i].modifications) {
      if (
        !modifications.hasOwnProperty(modification.modification_unimod_name)
      ) {
        modifications[modification.modification_unimod_name] = {};
      }
      modifications[modification.modification_unimod_name][position_i] =
        modification.modification_classification;
    }
  }
  if (__data.positions.hasOwnProperty(position_j)) {
    for (let modification of __data.positions[position_j].modifications) {
      if (
        !modifications.hasOwnProperty(modification.modification_unimod_name)
      ) {
        modifications[modification.modification_unimod_name] = {};
      }
      modifications[modification.modification_unimod_name][position_j] =
        modification.modification_classification;
    }
  }
  let series = {
    type: "sankey",
    id: "detail@contact",
    top: "63.5%",
    left: "28.5%",
    height: "31.5%",
    width: "21%",
    draggable: false,
    levels: [
      { depth: 0, label: { position: "right", fontSize: 10 } },
      { depth: 1, label: { position: "left", fontSize: 10 } },
    ],
    data: [],
    links: [],
  };
  for (let modification_name of Object.keys(modifications)) {
    let in_i = false;
    let in_j = false;
    if (modifications[modification_name].hasOwnProperty(position_i)) {
      series.data.push({
        name: position_i + " " + modification_name,
        value: 1,
        depth: 0,
        itemStyle: {
          color: getColor(modifications[modification_name][position_i]),
        },
      });
      in_i = true;
    }
    if (modifications[modification_name].hasOwnProperty(position_j)) {
      series.data.push({
        name: position_j + " " + modification_name,
        value: 1,
        depth: 1,
        itemStyle: {
          color: getColor(modifications[modification_name][position_j]),
        },
      });
      in_j = true;
    }
    if (in_i && in_j) {
      series.links.push({
        source: position_i + " " + modification_name,
        target: position_j + " " + modification_name,
        value: 1,
      });
    }
  }
  return [
    {
      id: "title_contact_detail",
      text:
        "Contact Modifications: " +
        "Positions " +
        position_i +
        " " +
        residue_i +
        " and " +
        position_j +
        " " +
        residue_j,
      textStyle: {
        fontSize: 12,
      },
      borderRadius: 4,
      backgroundColor: "#f2f3f2",
      top: "60%",
      left: "28%",
    },
    series,
  ];
}

function getColor(key) {
  if (!COLORS.M.hasOwnProperty(key)) {
    COLORS.M[key] = COLORS.i;
    COLORS.i += 1;
    if (COLORS.i > COLORS.C.length) {
      COLORS.i = 0;
    }
  }
  return COLORS.C[COLORS.M[key]];
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
