EXAMPLE_PARAMETER = ""; // Exemplary parameter.
WWW = "http://127.0.0.1:5000/";
ACTIVE_PANEL = "main-panel-1";
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

window.addEventListener("wheel", (event) => {
  if (event.deltaY > 0) {
    stepNext();
  } else {
    stepBack();
  }
});

/**
 *
 */
function init() {
  EXAMPLE_PARAMETER = null;
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
    $("#main-panel-1").first().slideToggle("medium");
    $("#main-panel-2").first().slideToggle("medium");
    $("#panel-back-button")[0].disabled = false;
    $("#panel-next-button")[0].disabled = true;
    $("#progress-icon-1")[0].classList.toggle("progress-active");
    $("#progress-icon-2")[0].classList.toggle("progress-active");
    ACTIVE_PANEL = "main-panel-2";
    initDashboard();
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
  DASHBOARD = echarts.init($("#main-panel-2-dashboard")[0], {
    devicePixelRatio: 2,
    renderer: "svg",
    width: "auto",
    height: "auto",
  });
  DASHBOARD.setOption(DASHBOARD_OPTION);
}

async function uploadData() {
  request = {
    contentType: null,
    content: null,
  };
  await readFile($("#data-input-form")[0].files[0]).then((response) => {
    request.content = response;
  });
  request.contentType = $("#data-type-form")[0].value;
  toggleProgress("Processing Data");
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
      console.log(response);
    })
    .catch((error) => {
      console.log(error);
    })
    .finally((_) => {
      toggleProgress(null);
      return;
    });
}
