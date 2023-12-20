_URL = "http://127.0.0.1:5000/";
_PROTEINS_OVERVIEW_TABLE = null;
// _DASHBOARD_SIZE_OBSERVER = null;
_MODIFICATIONS_GRAPH = null;
_MODIFICATIONS_GRAPH_SIZE_OBSERVER = null;
_DASHBOARD_STRUCTURE_VIEW = null;
_DASHBOARD_MAP_CHART = null;
_DASHBOARD_MAP_OBSERVER = null;
_DASHBOARD_CONTACT_MAP_OPTION = null;
_DASHBOARD_PRESENCE_MAP_OPTION = null;

/**
 * Inizializes all client side elements of the PTMVision application.
 */
function init() {
  _URL = API_PARAMETERS["URL"];
  $("#menu")[0].style.display = "flex";
  $("#panel")[0].style.display = "block";
  _PROTEINS_OVERVIEW_TABLE = new Tabulator("#panel-overview-table-tabulator", {
    height: "73%",
    maxHeight: "73%",
    selectable: 1,
    columns: [
      {
        title: "ID",
        field: "id",
        sorter: "string",
        width: "14%",
      },
      {
        title: "Name",
        field: "name",
        sorter: "string",
        width: "14%",
      },
      {
        title: "# Positions with Modifications",
        field: "modified_positions",
        sorter: "number",
        width: "35%",
      },
      {
        title: "# Distinct Modifications",
        field: "unique_modifications",
        sorter: "number",
        width: "35%",
      },
      {
        title: "Modifications",
        field: "modifications",
        visible: false,
      },
    ],
  });
  _PROTEINS_OVERVIEW_TABLE.on(
    "rowSelectionChanged",
    function (data, rows, selected, deselected) {
      if (data.length > 0) {
        $("#dashboard-display-button")[0].disabled = false;
      } else {
        $("#dashboard-display-button")[0].disabled = true;
      }
    }
  );
  _MODIFICATIONS_GRAPH = echarts.init($("#panel-overview-graph-chart")[0], {
    devicePixelRatio: 2,
    renderer: "canvas",
    width: "auto",
    height: "auto",
  });
  _MODIFICATIONS_GRAPH_SIZE_OBSERVER = new ResizeObserver((entries) => {
    _MODIFICATIONS_GRAPH.resize({
      width: entries[0].width,
      height: entries[0].height,
    });
  });
  _MODIFICATIONS_GRAPH_SIZE_OBSERVER.observe(
    $("#panel-overview-graph-chart")[0]
  );
  _DASHBOARD_STRUCTURE_VIEW = $3Dmol.createViewer($("#dashboard-structure"), {
    backgroundColor: "#FAFAFC",
    antialias: true,
    cartoonQuality: 6,
  });
}

function _set_table_filters(_) {
  const id_values = Metro.getPlugin(
    "#panel-overview-table-filter-id",
    "select"
  ).val();
  const mod_values = Metro.getPlugin(
    "#panel-overview-table-filter-modification",
    "select"
  ).val();
  _PROTEINS_OVERVIEW_TABLE.clearFilter(true);
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
  _PROTEINS_OVERVIEW_TABLE.setFilter(filters);
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

/**
 * Downloads an example dataset to the client.
 */
function downloadExampleData() {
  axios
    .get(_URL + "/example_data", { responseType: "blob" })
    .then((response) => {
      downloadBlob(response.data, "example_data.csv");
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
 * Displays a {@link sweetalert2} container yielding legal notice information.
 */
function redirectHelp() {
  Swal.fire({
    title: "Legal Notice",
    html: _redirect_help_html,
    width: "100%",
    position: "bottom",
    showCloseButton: true,
    showCancelButton: false,
    showConfirmButton: false,
    grow: true,
    heightAuto: true,
    color: "#333333",
    background: "#fafafcd9",
    backdrop: `
      rgba(239, 240, 248, 0.1)
      left top
      no-repeat
    `,
  });
}

/**
 * Animates the progress indicator icon and displays a custom message.
 *
 * @param {String} hint Progress message to display.
 */
function toggleProgress(hint) {
  $("#process-indicator").first().toggleClass("fa-flip");
  if (hint) {
    $("#process-indicator").first().attr({ "data-hint-text": hint });
  } else {
    $("#process-indicator")
      .first()
      .attr({ "data-hint-text": "No Process Running" });
  }
}

/**
 * Sends the specified search enginge output data to the PTMVision backend and loads the results in the overview table.
 */
async function uploadData() {
  toggleProgress("Uploading Data");
  request = {
    contentType: null,
    content: null,
  };
  if ($("#data-input-form")[0].files.length == 0) {
    handleError("No search engine output data was supplied.");
    toggleProgress();
    return;
  }
  await readFile($("#data-input-form")[0].files[0]).then((response) => {
    request.content = response;
  });
  request.contentType = $("#data-type-form")[0].value;
  axios
    .post(
      _URL + "/process_search_engine_output",
      pako.deflate(JSON.stringify(request)),
      {
        headers: {
          "Content-Type": "application/octet-stream",
          "Content-Encoding": "zlib",
        },
      }
    )
    .then((_) => {
      Metro.toast.create(
        "Your data has been uploaded successfully.",
        null,
        5000
      );
      // Init. table.
      axios.get(_URL + "/get_available_proteins").then((response) => {
        _PROTEINS_OVERVIEW_TABLE.setData(response.data);
        modifications = new Set();
        identifiers = new Set();
        for (let entry of response.data) {
          entry.modifications.split("$").forEach((m) => modifications.add(m));
          identifiers.add(entry.id + "$" + entry.name);
        }
        modifications = [...modifications];
        modifications_data_string = "";
        modifications.forEach((m) => {
          modifications_data_string +=
            `<option value="` + m + `">` + m + `</option>`;
        });
        Metro.getPlugin(
          "#panel-overview-table-filter-modification",
          "select"
        ).data(modifications_data_string);
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
        Metro.getPlugin("#panel-overview-table-filter-id", "select").data(
          identifiers_data_string
        );
      });
      // Init. graph.
      axios.get(_URL + "/get_modifications_graph").then((response) => {
        const opt = response.data;
        opt["tooltip"]["formatter"] = (params, ticket, callback) => {
          if (params.dataType == "edge") {
            return (
              `<b><u>` +
              params.data.source +
              `</u> and ` +
              `<u>` +
              params.data.target +
              `</u></b>
            <br>` +
              params.data.value +
              ` joint occurrences.`
            );
          } else if (params.dataType == "node") {
            return (
              `<b>Modification <u>` +
              params.data.name +
              `</u></b>
            <br>` +
              params.data.count +
              ` occurrences. Occurs<br>in ` +
              params.data.value +
              `% of proteins at least once.`
            );
          }
        };
        opt["series"][0]["itemStyle"]["color"] = (params) => {
          let R = 0; // to 255
          let G = 190; // to 60
          let B = 210; // to 0
          let Rs = 2.55;
          let Gs = 1.3;
          let Bs = 2.1;
          let frequency = params.data.frequency;
          return (
            "rgb(" +
            (R + Rs * frequency) +
            "," +
            (G - Gs * frequency) +
            "," +
            (B - Bs * frequency) +
            ")"
          );
        };
        _MODIFICATIONS_GRAPH.setOption(opt);
      });
    })
    .catch((error) => {
      handleError(error.message);
    })
    .finally((_) => {
      toggleProgress();
    });
}

function getDashboard(cutoff_value, pdb_text_value) {
  if (cutoff_value == null) {
    cutoff_value = 4.69;
  }
  if (pdb_text_value == null) {
    pdb_text_value = null;
  }
  toggleProgress("Initialize Dashboard");
  request = {
    uniprot_id: _PROTEINS_OVERVIEW_TABLE.getSelectedData()[0].id,
    pdb_text: pdb_text_value,
    cutoff: cutoff_value,
  };
  axios
    .post(_URL + "/get_protein_data", pako.deflate(JSON.stringify(request)), {
      headers: {
        "Content-Type": "application/octet-stream",
        "Content-Encoding": "zlib",
      },
    })
    .then((response) => {
      if (response.data.status == "Failed: No protein structure available.") {
        _requestPdb();
      } else {
        // Initialize single components.
        _DASHBOARD_STRUCTURE_VIEW = _getStructureView(
          response.data.content.structure,
          $("#dashboard-structure")
        );
        _DASHBOARD_CONTACT_MAP_OPTION = _getContactMapOption(
          response.data.content
        );
        _DASHBOARD_PRESENCE_MAP_OPTION = _getPresenceMapOption(
          response.data.content
        );
        [_DASHBOARD_MAP_CHART, _DASHBOARD_MAP_OBSERVER] =
          _constructEChartInstance(
            $("#dashboard-map")[0],
            _DASHBOARD_CONTACT_MAP_OPTION
          );
      }
    })
    .catch((error) => {
      handleError(error.message);
      console.log(error);
    })
    .finally((_) => {
      toggleProgress();
    });
}

function _requestPdb() {
  Swal.fire({
    title: "Failed to fetch Protein Structure!",
    html: _upload_pdb_html,
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
      getDashboard(null, readFile($("#optional-pdb-input")[0].files[0]));
    }
  });
}

function _constructEChartInstance(html_container, echart_option) {
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
  chart.setOption(echart_option);
  return [chart, observer];
}

var _upload_pdb_html = `
    <p>You can provide a structure in .pdb format to continue.</p>
    <br>
    <input id="optional-pdb-input" type="file" data-role="file" data-mode="drop">`;

var _redirect_help_html = `
    <div>
      <div class="remark m-4">
        <h2 class="p-2 text-center">
          <small>3rd Party Libraries and Resources</small>
        </h2>
        <p class="text-left">
          <i class="fa-solid fa-circle-chevron-right"></i> We use
          <a href="https://sweetalert2.github.io/">sweetalert2.js</a> for popups.
        </p>
        <p class="text-left">
          <i class="fa-solid fa-circle-chevron-right"></i> The styling of our
          platform is based on <a href="https://metroui.org.ua/">Metro 4</a> and we
          use
          <a href="https://github.com/normanzb/resize-sensor#readme"
            >Element Sensor</a
          >
          for our layout. Many of the used icons are from
          <a href="https://fontawesome.com/">fontawesome</a>.
        </p>
        <p class="text-left">
          <i class="fa-solid fa-circle-chevron-right"></i>
          <a href="https://axios-http.com/docs/intro">Axios</a> is used for HTTP
          requests and communication with our server.
        </p>
        <p class="text-left">
          <i class="fa-solid fa-circle-chevron-right"></i> We utilize
          <a href="https://jquery.com/">jQuery.js</a> in our scripts and use
          <a href="https://github.com/foliojs/brotli.js">brotli.js</a> and
          <a href="https://github.com/nodeca/pako">pako.js</a> for data compression.
        </p>
        <p class="text-left">
          <i class="fa-solid fa-circle-chevron-right"></i> We utilize
          <a href="https://echarts.apache.org/en/index.html">Apache ECharts</a> for
          our visualization.
        </p>
        <p class="text-left">
          <i class="fa-solid fa-circle-chevron-right"></i> Our backend is based on
          <a href="https://flask.palletsprojects.com/en/2.2.x/">Flask</a> and
          <a href="https://flask-session.readthedocs.io/en/latest/">Flask-Session</a
          >.
        </p>
      </div>
      <div class="remark m-4">
        <h2 class="p-2 text-center">
          <small>Disclaimer</small>
        </h2>
        <p class="indent text-just">
          <b>Content liability</b> The contents of our pages were created with the
          utmost care. However, we cannot guarantee the accuracy, completeness and
          actuality of the content. As a service provider, we are responsible for
          our own content on these pages under the general laws according to § 7
          TMG. According to §§ 8 to 10 TMG, however, we are not obligated as a
          service provider to monitor transmitted or stored third-party information
          or to investigate circumstances that indicate illegal activity.
          Obligations to remove or block the use of information under the general
          laws remain unaffected. However, liability in this regard is only possible
          from the point in time at which a concrete infringement of the law becomes
          known. If we become aware of such infringements, we will remove this
          content immediately.
        </p>
        <p class="indent text-just">
          <b>Link liability</b> This page may contain links to external websites of
          third parties, on whose contents we have no influence. We cannot take any
          liability for these external contents. The respective provider or operator
          of the pages is always responsible for the content of the linked pages.
          The linked pages were checked for possible legal violations at the time of
          linking. Illegal contents were not recognizable at the time of linking.
          However, a permanent control of the contents of the linked pages is not
          reasonable without concrete evidence of a violation of the law. If we
          become aware of any infringements, we will remove such links immediately.
        </p>
        <p class="indent text-just">
          <b>Data privacy</b> The use of our website is possible without providing
          personal data and no processed data is stored. If personal data (such as
          name, address or e-mail addresses) is collected, this is on a voluntary
          basis. This data will not be passed on to third parties without your
          express consent. We point out that data transmission over the Internet
          (e.g. communication by e-mail) may yield security gaps. Complete
          protection of data against access by third parties is not possible. The
          use of contact data published within the framework of the impressum by
          third parties for the purpose of sending advertising and information
          material is hereby expressly prohibited. The operators of the pages
          expressly reserve the right to take legal action in the event of the
          sending of advertising information, such as spam e-mails.
        </p>
        <p class="indent text-just">
          <b>Copyright</b> This program is free software: you can redistribute it
          and/or modify it under the terms of the GNU General Public License as
          published by the Free Software Foundation, either version 3 of the
          License, or (at your option) any later version. This program is
          distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
          without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
          PARTICULAR PURPOSE. See the GNU General Public License for more details.
          (<a href="https://www.gnu.org/licenses/gpl-3.0.en.html" target="_blank">www.gnu.org/licenses/gpl-3.0.en</a>)
        </p>
      </div>
      <div class="remark m-4">
        <h2 class="p-2 text-center">
          <small>Cookie Usage</small>
        </h2>
        <p class="indent text-just">
          <b>We use a strictly necessary (http-only non-tracking) session cookie so that you can access your data!</b>
          (1) Your data will not be accessible to other people. (2) No personalized data is collected or shared with third parties and you will not be tracked by this cookie. (3) The cookie will expire after you close your tab.
        </p>
      </div>
      <div class="remark m-4">
        <h2 class="p-2 text-center">
          <small>Contact Information</small>
        </h2>
        <address>
          <b>Simon Hackl & Caroline Jachmann</b>
          <br />
          Sand 14<br />
          72076 Tübingen, Germany<br />
          <a href="mailto:simon.hackl@uni-tuebingen.de"
            >simon.hackl@uni-tuebingen.de</a
          > or <a href="mailto:caroline.jachmann@uni-tuebingen.de"
            >caroline.jachmann@uni-tuebingen.de</a
          >
        </address>
      </div>
    </div>`;

function _getStructureView(pdb_text, html_container) {
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
  view.setStyle(
    {},
    {
      cartoon: {
        color: "#d4d4d4",
      },
    }
  );
  view.render();
  return view;
}

function _getContactMapOption(_data) {
  let lower_data = [];
  let upper_data = [];
  for (let [residue_index_i, contacts_data] of Object.entries(_data.contacts)) {
    let x = parseInt(residue_index_i);
    for (let contact_data of contacts_data) {
      let y = contact_data[0];
      if (x > y) {
        upper_data.push([x, y, contact_data[1]]);
      } else {
        lower_data.push([x, y, contact_data[1]]);
      }
    }
  }
  return {
    title: {
      text: "Cβ Contact Map",
      textStyle: {
        fontSize: 12,
      },
    },
    legend: {
      bottom: 0,
      left: "center",
      orient: "horizontal",
      itemWidth: 10,
      itemHeight: 10,
      data: [],
      textStyle: {
        fontWeight: "lighter",
        fontSize: 10,
      },
      pageTextStyle: {
        fontWeight: "lighter",
        fontSize: 10,
      },
      pageIcons: {
        horizontal: [
          "M18.3 250.3c-3.1 3.1-3.1 8.2 0 11.3l216 216c3.1 3.1 8.2 3.1 11.3 0s3.1-8.2 0-11.3L35.3 256 245.7 45.7c3.1-3.1 3.1-8.2 0-11.3s-8.2-3.1-11.3 0l-216 216z",
          "M301.7 250.3c3.1 3.1 3.1 8.2 0 11.3l-216 216c-3.1 3.1-8.2 3.1-11.3 0s-3.1-8.2 0-11.3L284.7 256 74.3 45.7c-3.1-3.1-3.1-8.2 0-11.3s8.2-3.1 11.3 0l216 216z",
        ],
      },
      width: "100%",
      backgroundColor: "#f0f5f5",
    },
    dataZoom: [
      {
        type: "inside",
        xAxisIndex: [0, 1],
        yAxisIndex: 0,
        throttle: 0,
      },
    ],
    toolbox: {
      feature: {
        saveAsImage: {
          pixelRatio: 4,
        },
      },
    },
    grid: {
      top: "10%",
      left: "10%",
      height: "80%",
      width: "80%",
    },
    xAxis: [
      {
        data: [...Array(_data.sequence.length).keys()].map((v) => v + 1),
        name: "Protein Position",
        nameGap: 25,
        nameLocation: "center",
        nameTextStyle: {
          fontWeight: "bold",
          fontSize: 11,
        },
        axisLabel: {
          formatter: (p) => {
            return parseInt(p) + 1;
          },
          fontWeight: "lighter",
          fontSize: 10,
        },
      },
    ],
    yAxis: [
      {
        data: [...Array(_data.sequence.length).keys()].map((v) => v + 1),
        inverse: true,
        name: "Protein Position",
        nameGap: 25,
        nameLocation: "center",
        nameTextStyle: {
          fontWeight: "bold",
          fontSize: 11,
        },
        axisLabel: {
          fontWeight: "lighter",
          fontSize: 10,
        },
      },
    ],
    axisPointer: {
      show: true,
      label: {
        backgroundColor: "rgba(51, 51, 51, 0.7)",
        fontWeight: "lighter",
        fontSize: 10,
        color: "#f0f5f5",
      },
      link: {
        xAxisIndex: [0, 1],
      },
      triggerTooltip: false,
    },
    tooltip: {
      formatter: (params, ticket, callback) => {
        return (
          `Residues ` +
          params.data[0] +
          ` and ` +
          params.data[1] +
          ` in contact (` +
          parseFloat(params.data[2]).toFixed(2) +
          `Å).`
        );
      },
      backgroundColor: "rgba(51, 51, 51, 0.7)",
      borderColor: "transparent",
      textStyle: {
        fontWeight: "lighter",
        fontSize: 10,
        color: "#f0f5f5",
      },
    },
    series: [
      {
        type: "heatmap",
        name: "Contacts (Lower Triangle)",
        data: lower_data,
        itemStyle: {
          color: "#333333",
          borderColor: "#FAFAFC",
          borderWidth: 0.1,
          borderRadius: 5,
        },
      },
    ],
  };
}

function _getPresenceMapOption(_data) {
  let modifications = [];
  for (let position_data of Object.values(_data.positions)) {
    for (let modification of Object.values(position_data.modifications)) {
      modifications.push(modification.modification_unimod_name);
    }
  }
  modifications = [...new Set(modifications)].sort();
  let data_series = {};
  for (let [position, position_data] of Object.entries(_data.positions)) {
    for (let modification of Object.values(position_data.modifications)) {
      if (
        !data_series.hasOwnProperty(modification.modification_classification)
      ) {
        data_series[modification.modification_classification] = {
          type: "heatmap",
          name: modification.modification_classification,
          data: [],
          xAxisIndex: 0,
          yAxisIndex: 0,
        };
      }
      let residue_index = parseInt(position);
      data_series[modification.modification_classification].data.push([
        residue_index,
        modifications.indexOf(modification.modification_unimod_name),
        1,
        modification.modification_unimod_name,
        modification.modification_classification,
        _data.sequence[residue_index],
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
      type: "heatmap",
      data: {},
      xAxisIndex: 1,
      yAxisIndex: 1,
      name: "Molecule processing",
      id: "ann_molecule_processing",
      itemStyle: {
        color: "#211A1D",
      },
    },
    Region: {
      type: "heatmap",
      data: {},
      xAxisIndex: 1,
      yAxisIndex: 1,
      name: "Region",
      id: "ann_region",
      itemStyle: {
        color: "#F0B67F",
      },
    },
    Site: {
      type: "heatmap",
      data: {},
      xAxisIndex: 1,
      yAxisIndex: 1,
      name: "Site",
      id: "ann_site",
      itemStyle: {
        color: "#57B3DE",
      },
    },
    "Amino acid modifications": {
      type: "heatmap",
      data: {},
      xAxisIndex: 1,
      yAxisIndex: 1,
      name: "Amino acid modifications",

      id: "ann_amino_acid_modifications",
      itemStyle: {
        color: "#E70D56",
      },
    },
    "Natural variations": {
      type: "heatmap",
      data: {},
      xAxisIndex: 1,
      yAxisIndex: 1,
      name: "Natural variations",
      id: "ann_natural_variations",
      itemStyle: {
        color: "#54B679",
      },
    },
    "Experimental info": {
      type: "heatmap",
      data: {},
      xAxisIndex: 1,
      yAxisIndex: 1,
      name: "Experimental info",
      id: "ann_experimental_info",
      itemStyle: {
        color: "#927C85",
      },
    },
    "Secondary structure": {
      type: "heatmap",
      data: {},
      xAxisIndex: 1,
      yAxisIndex: 1,
      name: "Secondary structure",
      id: "ann_secondary_structure",
    },
  };
  let annotation_secondary_structure_colors = {
    Helix: "#E02222",
    Turn: "#ADD46E",
    "Beta strand": "#226CE0",
  };
  let annotation_labels = [...Object.keys(annotation_series)];
  for (let annotation_entry of Object.values(_data.annotation.features)) {
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
              annotation_entry.type,
              annotation_entry.location.start.value,
              annotation_entry.location.end.value,
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

  return {
    title: {
      text: "Modifications Presence Map",
      textStyle: {
        fontSize: 12,
      },
    },
    legend: {
      bottom: 0,
      left: "center",
      orient: "horizontal",
      itemWidth: 10,
      itemHeight: 10,
      data: Object.values(data_series).map((S) => {
        return { name: S.name, icon: "roundRect" };
      }),
      textStyle: {
        fontWeight: "lighter",
        fontSize: 10,
      },
      pageTextStyle: {
        fontWeight: "lighter",
        fontSize: 10,
      },
      pageIcons: {
        horizontal: [
          "M18.3 250.3c-3.1 3.1-3.1 8.2 0 11.3l216 216c3.1 3.1 8.2 3.1 11.3 0s3.1-8.2 0-11.3L35.3 256 245.7 45.7c3.1-3.1 3.1-8.2 0-11.3s-8.2-3.1-11.3 0l-216 216z",
          "M301.7 250.3c3.1 3.1 3.1 8.2 0 11.3l-216 216c-3.1 3.1-8.2 3.1-11.3 0s-3.1-8.2 0-11.3L284.7 256 74.3 45.7c-3.1-3.1-3.1-8.2 0-11.3s8.2-3.1 11.3 0l216 216z",
        ],
      },
      width: "100%",
      backgroundColor: "#f0f5f5",
    },
    dataZoom: [
      {
        type: "inside",
        xAxisIndex: [0, 1],
        yAxisIndex: 0,
        throttle: 0,
      },
    ],
    toolbox: {
      feature: {
        saveAsImage: {
          pixelRatio: 4,
        },
      },
    },
    grid: [
      {
        top: "22%",
        left: "18%",
        height: "66%",
        width: "80%",
      },
      {
        top: "5%",
        left: "18%",
        height: "15%",
        width: "80%",
        backgroundColor: "#769213",
      },
    ],
    xAxis: [
      {
        data: [...Array(_data.sequence.length).keys()],
        name: "Protein Position",
        nameGap: 25,
        nameLocation: "center",
        nameTextStyle: {
          fontWeight: "bold",
          fontSize: 11,
        },
        axisLabel: {
          formatter: (p) => {
            return parseInt(p) + 1;
          },
          fontWeight: "lighter",
          fontSize: 10,
        },
        gridIndex: 0,
      },
      {
        data: [...Array(_data.sequence.length).keys()],
        axisLabel: {
          show: false,
        },
        gridIndex: 1,
      },
    ],
    yAxis: [
      {
        data: modifications,
        inverse: true,
        name: "Modification",
        nameGap: 140,
        nameLocation: "center",
        nameTextStyle: {
          fontWeight: "bold",
          fontSize: 11,
        },
        axisLabel: {
          fontWeight: "lighter",
          fontSize: 10,
        },
        gridIndex: 0,
      },
      {
        data: annotation_labels,
        name: "Annotations",
        nameGap: 140,
        nameLocation: "center",
        nameTextStyle: {
          fontWeight: "bold",
          fontSize: 11,
        },
        axisLabel: {
          fontWeight: "lighter",
          fontSize: 9,
        },
        interval: 0,
        gridIndex: 1,
      },
    ],
    axisPointer: {
      show: true,
      label: {
        backgroundColor: "rgba(51, 51, 51, 0.7)",
        fontWeight: "lighter",
        fontSize: 10,
        color: "#f0f5f5",
      },
      link: {
        xAxisIndex: [0, 1],
      },
      triggerTooltip: false,
    },
    tooltip: {
      formatter: (params, ticket, callback) => {
        let content;
        if (params.seriesId.startsWith("ann_")) {
          if (params.seriesId == "ann_secondary_structure") {
            content =
              `Secondary structure <b>` +
              params.data.value[3][0] +
              `</b>, Positions ` +
              params.data.value[3][1] +
              ` to ` +
              params.data.value[3][2];
          } else {
            content = params.seriesName + ` annotation`;
            for (let entry of params.data[3]) {
              content +=
                `<hr><br>` +
                entry[0] +
                ` (` +
                (entry[3] == "" ? "no description" : entry[3]) +
                `), Positions ` +
                entry[1] +
                ` to ` +
                entry[2];
            }
          }
        } else {
          content =
            `<b>Position: ` +
            params.data[0] +
            `<br>Aminoacid: ` +
            params.data[5] +
            `</b><br>Modification: ` +
            params.data[3] +
            ` [` +
            params.data[4] +
            `]`;
        }
        return content;
      },
      backgroundColor: "rgba(51, 51, 51, 0.7)",
      borderColor: "transparent",
      textStyle: {
        fontWeight: "lighter",
        fontSize: 10,
        color: "#f0f5f5",
      },
    },
    color: [
      "#E54B4B",
      "#05668D",
      "#167911",
      "#654597",
      "#FFA630",
      "#3BA9C7",
      "#7DC95E",
      "#AB81CD",
      "#F4D35E",
      "#A0C6E8",
      "#BBDB9B",
      "#E2ADF2",
    ],
    series: Object.values(data_series).concat(annotation_series_list),
  };
}

function _generatePresenceMapOption() {}
