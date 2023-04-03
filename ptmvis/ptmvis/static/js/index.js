EXAMPLE_PARAMETER = ""; // Exemplary parameter.
WWW = "";

/**
 *
 */
function init() {
  EXAMPLE_PARAMETER = null;
  WWW = API_PARAMETERS["URL"];
  $("#menu")[0].style.display = "flex";
  $(".main")[0].style.display = "block";
}

function downloadBlob(blob, name) {
  var download_link = document.createElement("a");
  download_link.href = window.URL.createObjectURL(new Blob([blob]));
  download_link.download = name;
  download_link.click();
  download_link.remove();
}
