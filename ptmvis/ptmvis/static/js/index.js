WWW = "http://127.0.0.1:5000/";
ACTIVE_PANEL = "main-panel-1";
PTM_DATA = null;
CONTACT_DATA = null;
DASHBOARD_OPTION = {
  grid: [
    {
      id: "grid_primary_sequence_profile",
      show: true,
      left: 0,
      top: 0,
      right: "38%",
      bottom: "85%",
      backgroundColor: "#ECC8AF",
    },
    {
      id: "grid_primary_sequence_annotation",
      show: true,
      left: 0,
      top: "15%",
      right: "38%",
      bottom: "80%",
      backgroundColor: "#E7AD99",
    },
    {
      id: "grid_contact_map",
      show: true,
      left: 0,
      top: "20%",
      right: "38%",
      bottom: 0,
      backgroundColor: "#CE796B",
    },
    {
      id: "grid_composition_one",
      show: true,
      left: "61%",
      top: 0,
      right: "19%",
      bottom: "50%",
      backgroundColor: "#495867",
    },
    {
      id: "grid_composition_two",
      show: true,
      left: "61%",
      top: "50%",
      right: "19%",
      bottom: 0,
      backgroundColor: "#495867",
    },
    {
      id: "grid_residue_view",
      show: true,
      left: "80%",
      top: 0,
      right: 0,
      bottom: 0,
      backgroundColor: "#420039",
    },
  ],
};
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
  axios.get(WWW + "/get_available_proteins").then((response) => {
    let options = {};
    response.data.forEach((entry) => {
      options[entry] = entry;
    });
    Metro.getPlugin("#main-panel-2-protein-select", "select").data(options);
  });
  DASHBOARD = echarts.init($("#main-panel-2-dashboard")[0], {
    devicePixelRatio: 2,
    renderer: "svg",
    width: "auto",
    height: "auto",
  });
  DASHBOARD.setOption(DASHBOARD_OPTION);
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
      console.log(response);
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
                      console.log(response.data);
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
        console.log(response.data);
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
