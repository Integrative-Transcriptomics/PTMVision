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
