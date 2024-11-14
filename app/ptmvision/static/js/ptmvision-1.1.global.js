/*!
 * PTMVision v1.1 (https://github.com/Integrative-Transcriptomics/PTMVision)
 * Copyright 2023-2024 by Caroline Jachmann and Simon Hackl
 * Licensed under GPL-3.0
 !*/

/**
 * Redirects the browser to the specified page of this project.
 *
 * This function is necessary as the single pages of the application do not know the base URL of the project.
 *
 * @param {String} pageName The name of the page to redirect to.
 * @param {String} target The target of the redirection; Should be one of '_blank' or '_self'.
 */
function redirectTo(pageName, target) {
  window.open(window.location.origin + "/" + pageName, target);
}

/**
 * Displays a notification to the user.
 *
 * @param {String} text The text to display in the notification.
 */
function displayNotification(text) {
  $("#menu").append(
    `<div class='notification'><i class="fa-duotone fa-spinner-third fa-spin fa-2xl"></i> ` +
      text +
      `</div>`
  );
}

/**
 * Removes all notifications from the user interface.
 */
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

/**
 * Scrolls to a specific DOM element.
 *
 * @param {String} id DOM element identifier to scroll to.
 */
function scroll_to(id) {
  let offset = $("#" + id).get(0).offsetTop - 40;
  window.scrollTo({
    top: offset,
    behavior: "smooth",
  });
}

/**
 * Downloads the current session data to the client.
 */
function downloadSessionData() {
  axios
    .get(window.location.origin + "/download_session")
    .then((response) => {
      const D = new Date();
      downloadBlob(
        response.data,
        "ptmvision-" +
          [D.getFullYear(), D.getMonth() + 1, Date.now()].join("-") +
          ".zlib"
      );
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

/**
 * Sends a request to download resources from the server.
 *
 * @param {String} name The full file name to request from the server.
 */
function downloadResource(name) {
  axios
    .get(window.location.origin + "/resource?name=" + name)
    .then((response) => {
      downloadBlob(response.data, name);
    });
}
