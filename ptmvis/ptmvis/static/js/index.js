SUCCESS_CODE = ""; // Code returned by the server for succ. requests.
ERROR_CODE = ""; // Code returned by the server for faulty requests.
APPLICATION_ISSUE_CODE = ""; // Code returned by the server, if session result data is faulty.
RESULT = ""; // Server session key to store results.
APPLICATION_ERROR_LOG = ""; // Server session key to store application error log.
SERVER_ERROR_LOG = ""; // Server session key to store server error log.
APPLICATION_RUN_LOG = ""; // Server session key to store application run log.
WWW = "";

/**
 *
 */
function init() {
  SUCCESS_CODE = API_PARAMETERS["SUCCESS_CODE"];
  ERROR_CODE = API_PARAMETERS["ERROR_CODE"];
  APPLICATION_ISSUE_CODE = API_PARAMETERS["APPLICATION_ISSUE_CODE"];
  RESULT = API_PARAMETERS["RESULT_KEY"];
  APPLICATION_ERROR_LOG = API_PARAMETERS["APPLICATION_ERROR_LOG_KEY"];
  SERVER_ERROR_LOG = API_PARAMETERS["SERVER_ERROR_LOG_KEY"];
  APPLICATION_RUN_LOG = API_PARAMETERS["APPLICATION_RUN_LOG_KEY"];
  WWW = API_PARAMETERS["URL"];
  checkForSession();
  $("#menu")[0].style.display = "flex";
  $(".main")[0].style.display = "block";
}

/**
 *
 * @param {*} response
 */
function handleResponseCode(response) {
  if (response.data == ERROR_CODE) {
    Swal.fire({
      iconHtml: `<i class="error-icon fa-solid fa-triangle-exclamation"></i>`,
      title: "Server Error",
      confirmButtonColor: "#6d81ad",
      color: "#747474",
      background: "#fafafcd9",
      backdrop: `
            rgba(96, 113, 150, 0.4)
            left top
            no-repeat
          `,
      html:
        `
        <div class="remark secondary text-left">
            A server side error occurred. Please check your input/session data. If you cannot solve your problem, feel free to open an issue <a href='https://github.com/Integrative-Transcriptomics/MUSIAL/issues' target='_blan'>here</a>.
        </div>
        <div class="remark secondary text-left"><span class="input-info-tag">LOG:</span><br>
            ` +
        response.data[SERVER_ERROR_LOG] +
        `
        </div>
      `,
    });
  } else if (response.data == APPLICATION_ISSUE_CODE) {
    axios
      .get(WWW + "/log")
      .then((response) => {
        Swal.fire({
          iconHtml: `<i class="error-icon fa-solid fa-bug"></i>`,
          title: "Request Error",
          confirmButtonColor: "#6d81ad",
          color: "#747474",
          background: "#fafafcd9",
          backdrop: `
            rgba(96, 113, 150, 0.4)
            left top
            no-repeat
          `,
          html:
            `
            <div class="remark secondary text-left">
                An application/request error occurred. Please check your input/session data. If you cannot solve your problem, feel free to open an issue <a href='https://github.com/Integrative-Transcriptomics/MUSIAL/issues' target='_blan'>here</a>.
            </div>
            <div class="remark secondary text-left"><span class="input-info-tag">LOG:</span><br>
            ` +
            response.data[APPLICATION_ERROR_LOG] +
            `
            </div>
          `,
        });
      })
      .catch((error) => {
        handleError(error);
      });
  }
}

/**
 *
 * @param {*} text
 */
function displayWarningPopup(text) {
  Swal.fire({
    icon: "warning",
    title: "Oops...",
    confirmButtonColor: "#6d81ad",
    text: text,
  });
}

function handleError(error) {
  // console.log( error );
  Swal.fire({
    iconHtml: `<i class="error-icon fa-solid fa-triangle-exclamation"></i>`,
    title: "Application Error",
    confirmButtonColor: "#6d81ad",
    color: "#747474",
    background: "#fafafcd9",
    backdrop: `
            rgba(96, 113, 150, 0.4)
            left top
            no-repeat
          `,
    html:
      `
        <div class="remark secondary text-left">
            An application error occurred. Please check your input/session data. If you cannot solve your problem, feel free to open an issue <a href='https://github.com/Integrative-Transcriptomics/MUSIAL/issues' target='_blan'>here</a>.
        </div>
        <div class="remark secondary text-left"><span class="input-info-tag">LOG:</span><br>
            ` +
      error.message +
      `
        </div>
      `,
  });
}

function downloadBlob(blob, name) {
  var download_link = document.createElement("a");
  download_link.href = window.URL.createObjectURL(new Blob([blob]));
  download_link.download = name;
  download_link.click();
  download_link.remove();
}

function checkForSession() {
  let CLR_NONE = "#eff0f8";
  let CLR_ACTIVE = "#39c093cc";
  let CLR_ISSUE = "#fe4848cc";
  axios
    .get(WWW + "/has_session")
    .then((response) => {
      if (response.data == SUCCESS_CODE) {
        $("#menu-active-session-indicator")[0].style.color = CLR_ACTIVE;
        $("#menu-link-results").attr("disabled", false);
      } else if (response.data == ERROR_CODE) {
        $("#menu-active-session-indicator")[0].style.color = CLR_NONE;
        $("#menu-link-results").attr("disabled", true);
      } else if (response.data == APPLICATION_ISSUE_CODE) {
        $("#menu-active-session-indicator")[0].style.color = CLR_ISSUE;
        $("#menu-link-results").attr("disabled", true);
      } else {
        $("#menu-active-session-indicator")[0].style.color = CLR_NONE;
        $("#menu-link-results").attr("disabled", true);
      }
    })
    .catch((error) => {
      handleError(error);
    });
}

function displayLoader(text) {
  Swal.fire({
    html:
      `
          <span class="tag-translucent">
          ` +
      text +
      `
          </span>
          <br>
          <div>
          <img src="` +
      LOADING_GIF +
      `" style="height: 200px; width: 200px">
          </div>
          `,
    width: "100%",
    padding: "3em",
    color: "#747474",
    background: "transparent",
    showConfirmButton: false,
    showCancelButton: false,
    allowOutsideClick: false,
    allowEscapeKey: false,
    backdrop: `
            rgba(239, 240, 248, 0.1)
            left top
            no-repeat
          `,
  });
}
