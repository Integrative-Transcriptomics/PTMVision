WWW = "http://127.0.0.1:5000/";
ACTIVE_PANEL = "main-panel-1";
STATE_AVAILABLE_PROTEINS_CHANGED = false; // This should be stored on the server in some appropriate way for long term usage.
AVAILABLE_PROTEINS = null;
DATA_PER_POSITION_PTMS = null;
DATA_CONTACTS = null;
DATA_SEQUENCE = null;
DASHBOARD_SIZE_OBSERVER = null;
DASHBOARD = null;

window.addEventListener("keyup", (event) => {
  if (event.key == "ArrowDown") {
    stepNext();
  } else if (event.key == "ArrowUp") {
    stepBack();
  }
});

/**
 *
 */
function init() {
  WWW = API_PARAMETERS["URL"];
  $("#menu")[0].style.display = "flex";
  $("#main-panel-1")[0].style.display = "block";
  DASHBOARD = echarts.init($("#main-panel-2-dashboard")[0], {
    devicePixelRatio: 2,
    renderer: "canvas",
    width: "auto",
    height: "auto",
  });
  DASHBOARD_SIZE_OBSERVER = new ResizeObserver((entries) => {
    DASHBOARD.resize({
      width: entries[0].width,
      height: entries[0].height,
    });
  });
  DASHBOARD_SIZE_OBSERVER.observe($("#main-panel-2-dashboard")[0]);
}

/**
 *
 * @param {*} blob
 * @param {*} name
 */
function downloadBlob(blob, name) {
  var download_link = document.createElement("a");
  download_link.href = window.URL.createObjectURL(new Blob([blob]));
  download_link.download = name;
  download_link.click();
  download_link.remove();
}

/**
 *
 * @param {*} file
 * @returns
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
 *
 * @param {*} msg
 */
function handleError(msg) {
  Metro.toast.create("Fatal Error: " + msg, null, 6000, "error-toast");
}

/**
 *
 */
function redirectSource() {
  window.open(
    "https://github.com/Integrative-Transcriptomics/MUSIAL",
    "_blank"
  );
}

/**
 *
 */
function redirectHelp() {
  window.open("https://www.google.com/webhp", "_blank");
}

/**
 *
 * @returns
 */
function stepBack() {
  if (ACTIVE_PANEL === "main-panel-1") {
    return;
  } else if (ACTIVE_PANEL === "main-panel-2") {
    $("#main-panel-2").first().slideToggle("medium");
    $("#main-panel-1").first().slideToggle("medium");
    $("#panel-back-button")[0].disabled = true;
    $("#panel-next-button")[0].disabled = false;
    $("#progress-icon-2")[0].classList.toggle("progress-active");
    $("#progress-icon-1")[0].classList.toggle("progress-active");
    ACTIVE_PANEL = "main-panel-1";
  }
}

/**
 *
 * @returns
 */
function stepNext() {
  if (ACTIVE_PANEL === "main-panel-1") {
    initDashboard();
    $("#main-panel-1").first().slideToggle("medium");
    $("#main-panel-2").first().slideToggle("medium");
    $("#panel-back-button")[0].disabled = false;
    $("#panel-next-button")[0].disabled = true;
    $("#progress-icon-1")[0].classList.toggle("progress-active");
    $("#progress-icon-2")[0].classList.toggle("progress-active");
    ACTIVE_PANEL = "main-panel-2";
  } else if (ACTIVE_PANEL === "main-panel-2") {
    return;
  }
}

/**
 *
 * @param {*} hint
 */
function toggleProgress(hint) {
  $("#progress-indicator").first().toggleClass("fa-spin");
  $("#progress-indicator").first().toggleClass("fa-circle");
  $("#progress-indicator").first().toggleClass("fa-circle-dashed");
  if (hint) {
    $("#progress-indicator").first().attr({ "data-hint-text": hint });
  } else {
    $("#progress-indicator")
      .first()
      .attr({ "data-hint-text": "No Progress Running" });
  }
}

/**
 *
 */
function initDashboard() {
  toggleProgress("Initialize Dashboard");
  if (STATE_AVAILABLE_PROTEINS_CHANGED) {
    axios.get(WWW + "/get_available_proteins").then((response) => {
      let options = {};
      for (const [k, v] of Object.entries(response.data["proteins"])) {
        options[k] = "(" + String(v) + ") " + k;
      }
      Metro.getPlugin("#main-panel-2-protein-select", "select").data(options);
    });
    STATE_AVAILABLE_PROTEINS_CHANGED = false;
  }
  toggleProgress();
}

/**
 *
 */
async function uploadData() {
  toggleProgress("Uploading Data");
  request = {
    contentType: null,
    content: null,
  };
  await readFile($("#data-input-form")[0].files[0]).then((response) => {
    request.content = response;
  });
  request.contentType = $("#data-type-form")[0].value;
  axios
    .post(
      WWW + "/process_search_engine_output",
      pako.deflate(JSON.stringify(request)),
      {
        headers: {
          "Content-Type": "application/octet-stream",
          "Content-Encoding": "zlib",
        },
      }
    )
    .then((response) => {
      Metro.toast.create(
        "Your data has been uploaded successfully.\nYou can proceed to the dashboard.",
        null,
        5000
      );
      STATE_AVAILABLE_PROTEINS_CHANGED = true;
    })
    .catch((error) => {
      handleError(error.message);
    })
    .finally((_) => {
      toggleProgress(null);
      return;
    });
}

/**
 *
 */
async function processDashboardRequest() {
  toggleProgress("Preparing Dashboard");
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
    .post(WWW + "/get_protein_data", pako.deflate(JSON.stringify(request)), {
      headers: {
        "Content-Type": "application/octet-stream",
        "Content-Encoding": "zlib",
      },
    })
    .then((response) => {
      if (response.data.status.startsWith("Failed")) {
        console.log(response);
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
                    WWW + "/get_protein_data",
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
                      updateDashboardInformation(response.data);
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
        updateDashboardInformation(response.data);
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

/**
 *
 * @param {*} data
 */
function updateDashboardInformation(data) {
  toggleProgress("Update Dashboard");
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
        if (position in per_position_ptms) {
          per_position_ptms[position].add(
            ptm + "$" + ptm_acc + "$" + ptm_classification
          );
        } else {
          per_position_ptms[position] = new Set([
            ptm + "$" + ptm_acc + "$" + ptm_classification,
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
  updateDashboardOption(
    DATA_PER_POSITION_PTMS,
    DATA_CONTACTS,
    DATA_SEQUENCE,
    DASHBOARD
  );
  Metro.toast.create(
    "Dashboard information was updated successfully.",
    null,
    5000
  );
  toggleProgress();
}
