_URL = "http://127.0.0.1:5000/";
_ACTIVE_PANEL = "main-panel-1";
_PROTEINS_OVERVIEW_DATA = null;
_PROTEINS_OVERVIEW_TABLE = null;
_PROTEIN_MODIFICATIONS_DATA = null;
_DASHBOARD_SIZE_OBSERVER = null;
_DASHBOARD = null;

window.addEventListener("keyup", (event) => {
  if (event.key == "ArrowDown") {
    stepNext();
  } else if (event.key == "ArrowUp") {
    stepBack();
  }
});

/**
 * Inizializes all client side elements of the PTMVision application.
 */
function init() {
  _URL = API_PARAMETERS["URL"];
  $("#menu")[0].style.display = "flex";
  $("#main-panel-1")[0].style.display = "block";
  _PROTEINS_OVERVIEW_TABLE = new Tabulator("#main-panel-1-table", {
    height: "88%",
    maxHeight: "88%",
    selectable: 1,
    columns: [
      {
        title: "Identifier",
        field: "id",
        sorter: "string",
        width: 150,
      },
      {
        title: "Gene Name",
        field: "name",
        sorter: "string",
        width: 150,
      },
      {
        title: "No. Modified Positions",
        field: "modified_positions",
        sorter: "number",
        width: 250,
      },
      {
        title: "No. Unique Modifications",
        field: "unique_modifications",
        sorter: "number",
        width: 250,
      },
      {
        title: "Modifications",
        field: "modifications",
        visible: false,
      },
    ],
  });
  _DASHBOARD = echarts.init($("#main-panel-2-dashboard")[0], {
    devicePixelRatio: 2,
    renderer: "canvas",
    width: "auto",
    height: "auto",
  });
  _DASHBOARD_SIZE_OBSERVER = new ResizeObserver((entries) => {
    _DASHBOARD.resize({
      width: entries[0].width,
      height: entries[0].height,
    });
  });
  _DASHBOARD_SIZE_OBSERVER.observe($("#main-panel-2-dashboard")[0]);
}

function _set_table_filter(_, _, values) {
  _PROTEINS_OVERVIEW_TABLE.clearFilter(true);
  const filters = [];
  for (let v of values) {
    filters.push({
      field: "modifications",
      type: "like",
      value: v,
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
    html: `
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
          (<a href="http://_URL.gnu.org/licenses">http://_URL.gnu.org/licenses</a>)
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
    </div>`,
    width: "80vw",
    padding: "0.5em",
    position: "bottom",
    showCancelButton: false,
    grow: true,
    heightAuto: true,
    confirmButtonColor: "#62a8ac",
    confirmButtonText: "Ok",
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
 * Switches back to the first panel, if the second application panel is currently displayed.
 */
function stepBack() {
  if (_ACTIVE_PANEL === "main-panel-2") {
    $("#main-panel-2").first().slideToggle("medium");
    $("#main-panel-1").first().slideToggle("medium");
    $("#panel-back-button")[0].disabled = true;
    $("#panel-next-button")[0].disabled = false;
    $("#progress-icon-2")[0].classList.toggle("progress-active");
    $("#progress-icon-1")[0].classList.toggle("progress-active");
    _ACTIVE_PANEL = "main-panel-1";
  }
}

/**
 * Switches to the second panel, if the first application panel is currently displayed.
 */
function stepNext() {
  return;
  if (_ACTIVE_PANEL === "main-panel-1") {
    init_Dashboard();
    $("#main-panel-1").first().slideToggle("medium");
    $("#main-panel-2").first().slideToggle("medium");
    $("#panel-back-button")[0].disabled = false;
    $("#panel-next-button")[0].disabled = true;
    $("#progress-icon-1")[0].classList.toggle("progress-active");
    $("#progress-icon-2")[0].classList.toggle("progress-active");
    _ACTIVE_PANEL = "main-panel-2";
  }
}

/**
 * Animates the progress indicator icon and displays a custom message.
 *
 * @param {String} hint Progress message to display.
 */
function toggleProgress(hint) {
  $("#progress-indicator").first().toggleClass("fa-flip");
  if (hint) {
    $("#progress-indicator").first().attr({ "data-hint-text": hint });
  } else {
    $("#progress-indicator")
      .first()
      .attr({ "data-hint-text": "No Progress Running" });
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
    toggleProgress(null);
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
      axios.get(_URL + "/get_available_proteins").then((response) => {
        console.log(response);
        _PROTEINS_OVERVIEW_TABLE.setData(response.data);
        modifications = new Set();
        for (let entry of response.data) {
          entry.modifications.split("$").forEach((m) => modifications.add(m));
        }
        Metro.getPlugin(
          "#main-panel-1-table-filter",
          "taginput"
        ).setAutocompleteList([...modifications]);
      });
    })
    .catch((error) => {
      handleError(error.message);
    })
    .finally((_) => {
      toggleProgress(null);
      return;
    });
}

function init_Dashboard() {
  toggleProgress("Initialize _Dashboard");
  if (STATE__PROTEINS_OVERVIEW_DATA_CHANGED) {
    axios.get(_URL + "/get__proteins_OVERVIEW_DATA").then((response) => {
      let options = {};
      for (const [k, v] of Object.entries(response.data["proteins"])) {
        options[k] = "(" + String(v) + ") " + k;
      }
      Metro.getPlugin("#main-panel-2-protein-select", "select").data(options);
    });
  }
  toggleProgress();
}

async function process_DashboardRequest() {
  toggleProgress("Preparing _Dashboard");
  request = {
    uniprot_id: Metro.getPlugin("#main-panel-2-protein-select", "select").val(),
    opt_pdb_text: null,
    distance_cutoff: $("#main-panel-2-distance-cutoff")[0].value,
    /*contact_attr: Metro.getPlugin(
      "#main-panel-2-contact-attribute",
      "select"
    ).val(),
    annotation_obj: null,*/
  };
  axios
    .post(_URL + "/get_protein_data", pako.deflate(JSON.stringify(request)), {
      headers: {
        "Content-Type": "application/octet-stream",
        "Content-Encoding": "zlib",
      },
    })
    .then((response) => {
      if (response.data.status.startsWith("Failed")) {
        Swal.fire({
          title: "Unable to match UniProt ID " + request.uniprot_id,
          icon: "question",
          iconHtml:
            '<i class="fa-solid fa-circle-exclamation fa-shake fa-xl" style="color: #ff6663;"></i>',
          html: upload_pdb_html,
          width: "100vw",
          padding: "0.5em",
          position: "bottom",
          showCancelButton: true,
          grow: true,
          heightAuto: true,
          cancelButtonColor: "#d4d4d4",
          cancelButtonText: "Cancel",
          confirmButtonColor: "#dc5754",
          confirmButtonText: "Ok",
          color: "#333333",
          background: "#f0f5f5",
          backdrop: `
            rgba(161, 210, 206, 0.2)
            left top
            no-repeat
          `,
        }).then(async (result) => {
          if (result.isConfirmed) {
            await readFile($("#optional-pdb-input")[0].files[0])
              .then((response) => {
                request.opt_pdb_text = response;
                axios
                  .post(
                    _URL + "/get_protein_data",
                    pako.deflate(JSON.stringify(request)),
                    {
                      headers: {
                        "Content-Type": "application/octet-stream",
                        "Content-Encoding": "zlib",
                      },
                    }
                  )
                  .then((response) => {
                    if (response.data.status == "Ok") {
                      update_DashboardInformation(response.data);
                    } else {
                      throw Error(
                        "Failed to parse provided .pdb structure: " +
                          response.data.status
                      );
                    }
                  })
                  .catch((error) => {
                    handleError(error.message);
                  });
              })
              .catch((error) => {
                handleError(error.message);
              });
          } else {
            return;
          }
        });
      } else if (response.data.status == "Ok") {
        update_DashboardInformation(response.data);
      }
    })
    .catch((error) => {
      handleError(error.message);
    })
    .finally((_) => {
      toggleProgress(null);
      return;
    });
}

var upload_pdb_html = `
    <p>You can provide a structure in .pdb format to continue.</p>
    <br>
    <input id="optional-pdb-input" type="file" data-role="file" data-mode="drop">`;

function update_DashboardInformation(data) {
  toggleProgress("Update _Dashboard");
  Papa.parse(data.ptms, {
    skipEmptyLines: "greedy",
    complete: (results) => {
      headers = results.data[0];
      parsed_entries = results.data.slice(1);
      per_position_ptms = {};
      for (const entry of parsed_entries) {
        let position = entry[2];
        let ptm = entry[3];
        let ptm_acc = entry[4];
        let ptm_classification = entry[5];
        let ptm_annotation = entry[6];
        if (position in per_position_ptms) {
          per_position_ptms[position].add(
            ptm +
              "$" +
              ptm_acc +
              "$" +
              ptm_classification +
              "$" +
              ptm_annotation
          );
        } else {
          per_position_ptms[position] = new Set([
            ptm +
              "$" +
              ptm_acc +
              "$" +
              ptm_classification +
              "$" +
              ptm_annotation,
          ]); // Store in Set to remove duplicates.
        }
      }
      for (const [k, v] of Object.entries(per_position_ptms)) {
        per_position_ptms[k] = Array.from(v);
      }
      DATA_PER_POSITION_PTMS = per_position_ptms;
    },
  });
  DATA_CONTACTS = data.contacts;
  DATA_SEQUENCE = data.sequence;
  update_DashboardOption(
    DATA_PER_POSITION_PTMS,
    DATA_CONTACTS,
    DATA_SEQUENCE,
    _DASHBOARD
  );
  Metro.toast.create(
    "_Dashboard information was updated successfully.",
    null,
    5000
  );
  toggleProgress();
}
