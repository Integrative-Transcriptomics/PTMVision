EXAMPLE_PARAMETER = ""; // Exemplary parameter.
WWW = "";
ACTIVE_PANEL = "main-panel-1";

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

function downloadBlob(blob, name) {
  var download_link = document.createElement("a");
  download_link.href = window.URL.createObjectURL(new Blob([blob]));
  download_link.download = name;
  download_link.click();
  download_link.remove();
}

function redirectSource() {
  window.open(
    "https://github.com/Integrative-Transcriptomics/MUSIAL",
    "_blank"
  );
}

function redirectHelp() {
  window.open("https://www.google.com/webhp", "_blank");
}

function stepBack() {
  if (ACTIVE_PANEL === "main-panel-1") {
    return;
  } else if (ACTIVE_PANEL === "main-panel-2") {
    $("#main-panel-2").first().slideToggle("medium");
    $("#main-panel-1").first().slideToggle("medium");
    $("#panel-back-button")[0].disabled = true;
    $("#panel-next-button")[0].disabled = false;
    $("#progress-no-2")[0].classList.toggle("progress-active");
    $("#progress-icon-2")[0].classList.toggle("progress-active");
    $("#progress-no-1")[0].classList.toggle("progress-active");
    $("#progress-icon-1")[0].classList.toggle("progress-active");
    ACTIVE_PANEL = "main-panel-1";
  } else if (ACTIVE_PANEL === "main-panel-3") {
    $("#main-panel-3").first().slideToggle("medium");
    $("#main-panel-2").first().slideToggle("medium");
    $("#panel-back-button")[0].disabled = false;
    $("#panel-next-button")[0].disabled = false;
    $("#progress-no-3")[0].classList.toggle("progress-active");
    $("#progress-icon-3")[0].classList.toggle("progress-active");
    $("#progress-no-2")[0].classList.toggle("progress-active");
    $("#progress-icon-2")[0].classList.toggle("progress-active");
    ACTIVE_PANEL = "main-panel-2";
  }
}

function stepNext() {
  if (ACTIVE_PANEL === "main-panel-1") {
    $("#main-panel-1").first().slideToggle("medium");
    $("#main-panel-2").first().slideToggle("medium");
    $("#panel-back-button")[0].disabled = false;
    $("#panel-next-button")[0].disabled = false;
    $("#progress-no-1")[0].classList.toggle("progress-active");
    $("#progress-icon-1")[0].classList.toggle("progress-active");
    $("#progress-no-2")[0].classList.toggle("progress-active");
    $("#progress-icon-2")[0].classList.toggle("progress-active");
    ACTIVE_PANEL = "main-panel-2";
  } else if (ACTIVE_PANEL === "main-panel-2") {
    $("#main-panel-2").first().slideToggle("medium");
    $("#main-panel-3").first().slideToggle("medium");
    $("#panel-back-button")[0].disabled = false;
    $("#panel-next-button")[0].disabled = true;
    $("#progress-no-2")[0].classList.toggle("progress-active");
    $("#progress-icon-2")[0].classList.toggle("progress-active");
    $("#progress-no-3")[0].classList.toggle("progress-active");
    $("#progress-icon-3")[0].classList.toggle("progress-active");
    ACTIVE_PANEL = "main-panel-3";
  } else if (ACTIVE_PANEL === "main-panel-3") {
    return;
  }
}

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
