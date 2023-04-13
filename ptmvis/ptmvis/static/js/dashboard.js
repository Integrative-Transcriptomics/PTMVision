var _styles_data_zoom = {
  backgroundColor: "transparent",
  fillerColor: "transparent",
  borderColor: "#33333320",
  brushStyle: {
    color: "#33333320",
  },
  handleStyle: {
    color: "#33333320",
    borderColor: "#33333320",
  },
  moveHandleStyle: {
    color: "#33333320",
    borderColor: "#33333320",
  },
  moveHandleSize: 3,
  emphasis: {
    moveHandleStyle: {
      color: "#33333380",
      borderColor: "#33333380",
    },
    handleStyle: {
      color: "#33333380",
      borderColor: "#33333380",
    },
  },
  showDataShadow: false,
  realtime: false,
};
var _style_sunburst_composition = {
  radius: ["0%", "5%"],
  emphasis: {
    focus: "ancestor",
  },
  levels: [
    {
      r0: "0%",
      r: "1%",
      itemStyle: {
        color: "#222222",
      },
    },
    {
      r0: "5%",
      r: "10%",
      label: {
        show: false,
      },
      emphasis: {
        label: {
          show: true,
          rotate: 0,
          position: "outside",
          distance: 30,
        },
      },
      itemStyle: {
        borderWidth: 3,
        borderRadius: 3,
      },
    },
    {
      r0: "10%",
      r: "30%",
      label: {
        position: "inside",
        align: "right",
        padding: 3,
        silent: false,
        formatter: (value, index) => {
          return value.name.substr(0, 8) + "\u2026";
        },
        fontSize: 8,
      },
      itemStyle: {
        borderWidth: 2,
        borderRadius: 4,
      },
    },
  ],
};
var _ptms = null;
var _ptms_accessors = null;
var _contacts = null;
var _sequence = null;
var _board = null;
var _board_option = {
  title: [
    {
      id: "composition_one_text",
      text: "Modification Composition Residue 10",
      left: "79.5%",
      bottom: "77.5%",
      textStyle: {
        fontWeight: "normal",
        fontSize: 14,
      },
      show: false,
    },
    {
      id: "composition_two_text",
      text: "Modification Composition Residue 10",
      left: "79.5%",
      bottom: "42%",
      textStyle: {
        fontWeight: "normal",
        fontSize: 14,
      },
      show: false,
    },
  ],
  grid: [
    {
      id: "grid_primary_profile",
      show: true,
      left: "30%",
      bottom: "85.8%",
      height: "10%",
      width: "40%",
      backgroundColor: "#fafafa",
    },
    {
      id: "grid_primary_annotation",
      show: true,
      left: "30%",
      bottom: "80.8%",
      height: "5%",
      width: "40%",
      backgroundColor: "#d4d4d4",
    },
    {
      id: "grid_secondary_profile",
      show: true,
      left: "72.82%",
      bottom: "10%",
      height: "70.8%",
      width: "5.65%",
      backgroundColor: "#fafafa",
    },
    {
      id: "grid_secondary_annotation",
      show: true,
      left: "70%",
      bottom: "10%",
      height: "70.8%",
      width: "2.82%",
      backgroundColor: "#d4d4d4",
    },
    {
      id: "grid_map",
      show: true,
      left: "30%",
      bottom: "10%",
      height: "70.8%",
      width: "40%",
      backgroundColor: "#fafafa",
    },
    {
      id: "grid_composition_one",
      show: true,
      left: "79.47%",
      bottom: "45.4%",
      height: "35.4%",
      width: "20%",
      backgroundColor: "transparent",
    },
    {
      id: "grid_composition_two",
      show: true,
      left: "79.47%",
      bottom: "10%",
      height: "35.4%",
      width: "20%",
      backgroundColor: "transparent",
    },
    {
      id: "grid_extra_view",
      show: true,
      left: "6%",
      bottom: "10%",
      height: "70.8%",
      width: "20%",
      backgroundColor: "#fafafa",
    },
  ],
  toolbox: {
    id: "toolbox",
    show: true,
    itemSize: 20,
    itemGap: 20,
    iconStyle: {
      color: "#a1d2ce",
      borderColor: "transparent",
      borderWidth: 0,
    },
    emphasis: {
      iconStyle: {
        color: "#50858b",
        textFill: "#333333",
      },
    },
    left: "left",
    top: "top",
    feature: {
      saveAsImage: {
        icon: "path://M192 112H64c-17.7 0-32 14.3-32 32v80H145.1c22.1-38.3 63.5-64 110.9-64s88.7 25.7 110.9 64H480V96c0-17.7-14.3-32-32-32H277.3c-6.9 0-13.7 2.2-19.2 6.4l-46.9 35.2c-5.5 4.2-12.3 6.4-19.2 6.4zM32 256V416c0 17.7 14.3 32 32 32H448c17.7 0 32-14.3 32-32V256H380c2.6 10.2 4 21 4 32c0 70.7-57.3 128-128 128s-128-57.3-128-128c0-11 1.4-21.8 4-32H32zM0 416V144c0-35.3 28.7-64 64-64H192l46.9-35.2C250 36.5 263.5 32 277.3 32H448c35.3 0 64 28.7 64 64V416c0 35.3-28.7 64-64 64H64c-35.3 0-64-28.7-64-64zM352 288a96 96 0 1 0 -192 0 96 96 0 1 0 192 0zM64 48c0-8.8 7.2-16 16-16h64c8.8 0 16 7.2 16 16s-7.2 16-16 16H80c-8.8 0-16-7.2-16-16z",
        // Font Awesome Pro 6.4.0 by @fontawesome - https://fontawesome.com License - https://fontawesome.com/license (Commercial License) Copyright 2023 Fonticons, Inc.
        pixelRatio: 3,
      },
      dataZoom: {
        icon: {
          zoom: "path://M208 32a176 176 0 1 1 0 352 176 176 0 1 1 0-352zm0 384c51.7 0 99-18.8 135.3-50L484.7 507.3c6.2 6.2 16.4 6.2 22.6 0s6.2-16.4 0-22.6L366 343.3c31.2-36.4 50-83.7 50-135.3C416 93.1 322.9 0 208 0S0 93.1 0 208S93.1 416 208 416zM192 304c0 8.8 7.2 16 16 16s16-7.2 16-16V224h80c8.8 0 16-7.2 16-16s-7.2-16-16-16H224V112c0-8.8-7.2-16-16-16s-16 7.2-16 16v80H112c-8.8 0-16 7.2-16 16s7.2 16 16 16h80v80z",
          // Font Awesome Pro 6.4.0 by @fontawesome - https://fontawesome.com License - https://fontawesome.com/license (Commercial License) Copyright 2023 Fonticons, Inc.
          back: "path://M94.7 360.2c-3.2-5-8.7-8.2-14.7-8.2c-12.3 0-20.3 12.8-13.7 23.2C106 438.2 176.1 480 256 480c123.7 0 224-100.3 224-224S379.7 32 256 32c-56.1 0-107.4 20.6-146.7 54.7L78.6 56c-5.1-5.1-12.1-8-19.3-8C44.2 48 32 60.2 32 75.3V176c0 8.8 7.2 16 16 16H148.7c15.1 0 27.3-12.2 27.3-27.3c0-7.2-2.9-14.2-8-19.3l-36-36C165.5 81.1 208.7 64 256 64c106 0 192 86 192 192s-86 192-192 192c-67.6 0-127.1-35-161.3-87.8zM64 86.6L137.4 160H64V86.6z",
          // Font Awesome Pro 6.4.0 by @fontawesome - https://fontawesome.com License - https://fontawesome.com/license (Commercial License) Copyright 2023 Fonticons, Inc.
        },
        xAxisIndex: [0],
        yAxisIndex: [0, 1],
      },
    },
  },
  tooltip: {
    trigger: "item",
    textStyle: {
      fontSize: 9,
    },
    formatter: function (params, ticket, callback) {
      if (params.seriesId == "contact_map") {
        return _formatTooltipContactMap(params, ticket, callback);
      } else if (params.seriesId == "modifications_map") {
        return _formatTooltipModificationsMap(params, ticket, callback);
      }
    },
  },
  axisPointer: {
    show: true,
    triggerOn: "mousemove",
    triggerTooltip: false,
    link: [
      {
        xAxisIndex: [0, 2],
      },
      {
        yAxisIndex: [0, 3],
      },
    ],
    label: {
      formatter: (params) => {
        if (
          (params.axisDimension == "x" && params.axisIndex == 2) ||
          (params.axisDimension == "y" && params.axisIndex == 3)
        ) {
          let ptms = _ptms[params.value];
          return typeof ptms === "undefined" ? "0" : String(ptms.length);
        } else {
          return params.value;
        }
      },
    },
  },
  xAxis: [
    {
      id: "grid_map_xaxis",
      type: "category",
      data: [],
      gridIndex: 4,
      name: "Residue",
      nameLocation: "middle",
      nameGap: 40,
      nameTextStyle: {
        fontWeight: "bold",
        fontSize: 14,
      },
      axisLabel: {
        show: true,
        fontSize: 11,
      },
      axisTick: {
        alignWithLabel: true,
      },
      splitLine: {
        show: true,
        interval: 0,
        lineStyle: {
          width: 0.2,
          color: "#e0e0e0",
        },
      },
      hideOverlap: true,
    },
    {
      id: "grid_extra_xaxis",
      type: "category",
      data: [],
      gridIndex: 7,
      name: "Modification",
      nameLocation: "middle",
      nameGap: 100,
      nameTextStyle: {
        fontWeight: "bold",
        fontSize: 14,
      },
      position: "top",
      axisLabel: {
        rotate: -50,
        interval: 0,
        show: true,
        margin: 25,
        fontSize: 9,
        formatter: (value, index) => {
          return value.substr(0, 10) + "\u2026";
        },
        hideOverlap: true,
      },
      axisTick: {
        length: 15,
        lineStyle: { width: 0.5 },
        alignWithLabel: true,
        interval: "auto",
      },
      splitLine: {
        show: true,
        interval: 0,
        lineStyle: {
          width: 0.5,
        },
      },
    },
    {
      id: "grid_primary_profile_xaxis",
      type: "category",
      data: [],
      gridIndex: 0,
      show: false,
      hideOverlap: true,
    },
    {
      id: "grid_secondary_profile_xaxis",
      type: "value",
      min: 0,
      gridIndex: 2,
      name: "No. Unique Mod.",
      nameLocation: "middle",
      nameGap: 40,
      minInterval: 1,
      nameTextStyle: {
        fontWeight: "bold",
        fontSize: 14,
      },
      axisLabel: {
        show: true,
        fontSize: 11,
      },
      axisTick: {
        alignWithLabel: true,
      },
      splitLine: {
        show: true,
        interval: 0,
        lineStyle: {
          width: 0.2,
          color: "#e0e0e0",
        },
      },
      hideOverlap: true,
    },
  ],
  yAxis: [
    {
      id: "grid_map_yaxis",
      type: "category",
      data: [],
      gridIndex: 4,
      name: "Residue",
      nameLocation: "middle",
      nameGap: 40,
      nameTextStyle: {
        fontWeight: "bold",
        fontSize: 14,
      },
      axisLabel: {
        show: true,
        fontSize: 11,
      },
      axisTick: {
        alignWithLabel: true,
      },
      splitLine: {
        show: true,
        interval: 0,
        lineStyle: {
          width: 0.2,
          color: "#e0e0e0",
        },
      },
      hideOverlap: true,
      inverse: true,
    },
    {
      id: "grid_extra_yaxis",
      type: "category",
      data: [],
      gridIndex: 7,
      show: false,
      hideOverlap: true,
      inverse: true,
    },
    {
      id: "grid_primary_profile_yaxis",
      type: "value",
      min: 0,
      gridIndex: 0,
      name: "No. Unique Mod.",
      nameLocation: "middle",
      nameGap: 40,
      minInterval: 1,
      nameTextStyle: {
        fontWeight: "bold",
        fontSize: 14,
      },
      axisLabel: {
        show: true,
        fontSize: 11,
      },
      axisTick: {
        alignWithLabel: true,
      },
      splitLine: {
        show: true,
        lineStyle: {
          width: 0.2,
          color: "#e0e0e0",
        },
      },
      hideOverlap: true,
    },
    {
      id: "grid_secondary_profile_yaxis",
      type: "category",
      data: [],
      gridIndex: 2,
      show: false,
      hideOverlap: true,
      inverse: true,
    },
  ],
  dataZoom: [
    {
      id: "main_xzoom",
      type: "slider",
      xAxisIndex: [0, 2],
      left: "30%",
      bottom: "80.8%",
      height: "5%",
      width: "40%",
      ..._styles_data_zoom,
    },
    {
      id: "extra_xzoom",
      type: "slider",
      xAxisIndex: [1],
      left: "6%",
      bottom: "80.8%",
      height: "2%",
      width: "20%",
      ..._styles_data_zoom,
    },
    {
      id: "main_yzoom",
      type: "slider",
      yAxisIndex: [0, 1, 3],
      left: "70%",
      bottom: "10%",
      height: "70.8%",
      width: "2.82%",
      ..._styles_data_zoom,
    },
  ],
  visualMap: [
    {
      id: "map_visual",
      type: "continuous",
      show: true,
      inRange: {
        color: [
          "#333333",
          "#333333",
          "#333333",
          "#333333",
          "#9DCEFF",
          "#F67C5E",
          "#FF4112",
          "#FF4112",
        ],
      },
      outOfRange: {
        color: "#333333",
      },
      range: [0, 1],
      precision: 2,
      min: -1,
      max: 1,
      seriesIndex: [0],
    },
  ],
  series: [
    {
      id: "contact_map",
      name: "Contact Map",
      type: "heatmap",
      data: [],
      xAxisIndex: 0,
      yAxisIndex: 0,
      itemStyle: {
        borderWidth: 0.5,
        borderColor: "#fafafa",
        borderCap: "round",
        borderJoin: "round",
        borderRadius: 3,
      },
    },
    {
      id: "modifications_map",
      name: "Modifications Map",
      type: "heatmap",
      data: [],
      xAxisIndex: 1,
      yAxisIndex: 1,
      itemStyle: {
        color: "#333333",
        borderWidth: 0.5,
        borderColor: "#fafafa",
        borderCap: "round",
        borderJoin: "round",
        borderRadius: 3,
      },
    },
    {
      id: "primary_profile",
      name: "Primary Profile",
      type: "bar",
      data: [],
      xAxisIndex: 2,
      yAxisIndex: 2,
      itemStyle: {
        color: "#333333",
      },
    },
    {
      id: "secondary_profile",
      name: "Secondary Profile",
      type: "bar",
      data: [],
      xAxisIndex: 3,
      yAxisIndex: 3,
      itemStyle: {
        color: "#333333",
      },
    },
    {
      id: "sunburst_composition_one",
      type: "sunburst",
      data: [],
      center: ["90%", "38%"], //["89.47", "55.4%"],
      ..._style_sunburst_composition,
    },
    {
      id: "sunburst_composition_two",
      type: "sunburst",
      data: [],
      center: ["90%", "73%"], //["89.47", "55.4%"],
      ..._style_sunburst_composition,
    },
  ],
};

function updateDashboardOption(ptms, contacts, sequence, board) {
  _ptms = ptms;
  _ptms_accessors = {};
  _contacts = contacts;
  _sequence = sequence;
  _board = board;
  // Reset composition overviews.
  _board_option.title.filter((t) => t.id == "composition_one_text")[0].text =
    "";
  _board_option.title.filter((t) => t.id == "composition_two_text")[0].text =
    "";
  _board_option.series.filter(
    (s) => s.id == "sunburst_composition_one"
  )[0].data = [];
  _board_option.series.filter(
    (s) => s.id == "sunburst_composition_two"
  )[0].data = [];
  // Add click event listener.
  _board.on("click", function (params) {
    if (params.seriesId == "modifications_map") {
      _handleModificationsMapClick(params);
    } else if (params.seriesId == "contact_map") {
      _handleContactMapClick(params);
    }
  });
  // Contact map.
  _board_option.xAxis.filter((ax) => ax.id == "grid_map_xaxis")[0].data =
    Array.from({ length: _sequence.length }, (_, i) => i + 1);
  _board_option.yAxis.filter((ax) => ax.id == "grid_map_yaxis")[0].data =
    Array.from({ length: _sequence.length }, (_, i) => i + 1);
  let contact_map_data = [];
  for (const [i, js] of Object.entries(_contacts)) {
    for (const j of js) {
      let posi = parseInt(i);
      let posj = parseInt(j);
      if (posi < posj) {
        contact_map_data.push([posi, posj, -1]);
      } else {
        contact_map_data.push([posi, posj, _applyCommonPTMSScore(i, j)]);
      }
    }
  }
  _board_option.series.filter((srs) => srs.id == "contact_map")[0].data =
    contact_map_data;
  // Mod. map.
  let per_modification_sites = {};
  for (let [modifications, sites] of Object.entries(_.invertBy(_ptms))) {
    for (let modification of modifications.split(",")) {
      let content = modification.split("$"),
        name = content[0],
        accessor = content[1];
      _ptms_accessors[name] = accessor;
      if (name in per_modification_sites) {
        per_modification_sites[name].push(...sites);
      } else {
        per_modification_sites[name] = sites;
      }
    }
  }
  let modifications = Array.from(Object.keys(per_modification_sites));
  _board_option.xAxis.filter((ax) => ax.id == "grid_extra_xaxis")[0].data =
    modifications;
  _board_option.yAxis.filter((ax) => ax.id == "grid_extra_yaxis")[0].data =
    Array.from({ length: _sequence.length }, (_, i) => i + 1);
  let modification_map_data = [];
  for (const [i, js] of Object.entries(per_modification_sites)) {
    for (const j of js) {
      modification_map_data.push([
        modifications.indexOf(i),
        parseInt(j) - 1,
        1,
      ]);
    }
  }
  _board_option.series.filter((srs) => srs.id == "modifications_map")[0].data =
    modification_map_data;
  // Primary Profile Simple.
  let per_site_modifications = [];
  for (let i = 1; i <= _sequence.length; i++) {
    if (i in _ptms) {
      per_site_modifications.push(_ptms[i].length);
    } else {
      per_site_modifications.push(0);
    }
  }
  _board_option.xAxis.filter(
    (ax) => ax.id == "grid_primary_profile_xaxis"
  )[0].data = Array.from({ length: _sequence.length }, (_, i) => i + 1);
  _board_option.series.filter((srs) => srs.id == "primary_profile")[0].data =
    per_site_modifications;
  // Secondary Profile Simple.
  _board_option.yAxis.filter(
    (ax) => ax.id == "grid_secondary_profile_yaxis"
  )[0].data = Array.from({ length: _sequence.length }, (_, i) => i + 1);
  _board_option.series.filter((srs) => srs.id == "secondary_profile")[0].data =
    per_site_modifications;
  // Update the board.
  _apply_option();
}

function _apply_option() {
  _board.showLoading();
  _board.setOption(_board_option, {
    notMerge: false,
    lazyUpdate: true,
    silent: true,
  });
  _board.hideLoading();
}

function _applyCommonPTMSScore(i, j) {
  let ptm_i = new Set(_ptms[parseInt(i) + 1]);
  let ptm_j = new Set(_ptms[parseInt(j) + 1]);
  if (ptm_i.size > 0 && ptm_j.size > 0) {
    let union = new Set([...ptm_i, ...ptm_j]);
    let intersect = new Set([...ptm_i].filter((e) => ptm_j.has(e)));
    return intersect.size / union.size;
  }
}

function _formatTooltipContactMap(params, ticket, callback) {
  let resx = parseInt(params.data[0]);
  let resy = parseInt(params.data[1]);
  let mode = parseFloat(params.data[2]);
  let suffix = "";
  if (mode >= 0) {
    suffix =
      `<h6><small>Jaccard Similarity:&nbsp;` +
      mode.toFixed(2) +
      `</small></h6>`;
    /*`&nbsp;<span class='donut'>` +
      mode +
      `/1</span></span>` + */
  }
  /*setTimeout(() => {
    document.querySelectorAll(".donut").forEach((e) =>
      peity(e, "donut", {
        fill: ["#6d81ad", "#fe4848cc"],
        innerRadius: 4,
        radius: 9,
      })
    );
  }, 0);*/
  return (
    `
      <h6><small>` +
    (resx + 1) +
    `&nbsp;(` +
    _sequence[resx] +
    `)&nbsp;<i class="fa-thin fa-arrows-left-right-to-line"></i>&nbsp;` +
    (resy + 1) +
    `&nbsp;(` +
    _sequence[resy] +
    `)</small></h6>` +
    suffix +
    `</br><span style="border-radius: 4px; background-color: #f1f1f1; padding: 1px;">Click for Details</span>`
  );
}

function _formatTooltipModificationsMap(params, ticket, callback) {
  let modification = params.name;
  let residue_index = parseInt(params.data[1]);
  return (
    `
      <h6><small>` +
    (residue_index + 1) +
    `&nbsp;(` +
    _sequence[residue_index] +
    `)&nbsp;<i class="fa-light fa-map-pin fa-rotate-90"></i>&nbsp;` +
    modification +
    `</br><span style="border-radius: 4px; background-color: #f1f1f1; padding: 1px;">Click to access UNIMOD</span>`
  );
}

function _handleModificationsMapClick(params) {
  window.open(
    "https://www.unimod.org/modifications_view.php?editid1=" +
      _ptms_accessors[params.name],
    "_blank"
  );
}

function _handleContactMapClick(params) {
  _setCompositionContent(params.data[0], params.data[1]);
}

function _setCompositionContent(resx_index, resy_index) {
  // Adjust composition titles.
  _board_option.title.filter((t) => t.id == "composition_one_text")[0].text =
    "PTMs of residue " +
    (resx_index + 1).toString() +
    ", " +
    _sequence[resx_index];
  _board_option.title.filter(
    (t) => t.id == "composition_one_text"
  )[0].show = true;
  _board_option.title.filter((t) => t.id == "composition_two_text")[0].text =
    "PTMs of residue " +
    (resy_index + 1).toString() +
    ", " +
    _sequence[resy_index];
  _board_option.title.filter(
    (t) => t.id == "composition_two_text"
  )[0].show = true;
  // Adjust sunburst composition view data.
  let ptms_resx = new Set(_ptms[resx_index + 1]);
  let ptms_resy = new Set(_ptms[resy_index + 1]);
  let predicate_common_ptm = (name) =>
    new Set(
      [...ptms_resx].filter((e) => ptms_resy.has(e)).map((e) => e.split("$")[0])
    ).has(name);
  let fill_to_composition = function (series, item_name, item_class) {
    if (series.length == 0) {
      series.push({
        name: item_class,
        value: 1,
        itemStyle: {
          color: item_class.toHex(),
        },
        children: [
          {
            name: item_name,
            value: 1,
            itemStyle: {
              color: predicate_common_ptm(item_name) ? "#FF4112" : "#9DCEFF",
            },
          },
        ],
      });
    } else {
      let matched_children = series.filter((entry) => entry.name == item_class);
      if (matched_children.length > 0) {
        matched_children[0].value += 1;
        matched_children[0].children.push({
          name: item_name,
          value: 1,
          itemStyle: {
            color: predicate_common_ptm(item_name) ? "#FF4112" : "#9DCEFF",
          },
        });
      } else {
        series.push({
          name: item_class,
          value: 1,
          itemStyle: {
            color: item_class.toHex(),
          },
          children: [
            {
              name: item_name,
              value: 1,
              itemStyle: {
                color: predicate_common_ptm(item_name) ? "#FF4112" : "#9DCEFF",
              },
            },
          ],
        });
      }
    }
  };
  // Fill first composition (top).
  let sunburst_composition_one_data = [];
  for (let ptm of ptms_resx) {
    let ptm_name = String(ptm.split("$")[0]);
    let ptm_class = String(ptm.split("$")[2]);
    fill_to_composition(sunburst_composition_one_data, ptm_name, ptm_class);
  }
  _board_option.series.filter(
    (s) => s.id == "sunburst_composition_one"
  )[0].data = sunburst_composition_one_data;
  // Fill second composition (bottom).
  let sunburst_composition_two_data = [];
  for (let ptm of ptms_resy) {
    let ptm_name = String(ptm.split("$")[0]);
    let ptm_class = String(ptm.split("$")[2]);
    fill_to_composition(sunburst_composition_two_data, ptm_name, ptm_class);
  }
  _board_option.series.filter(
    (s) => s.id == "sunburst_composition_two"
  )[0].data = sunburst_composition_two_data;
  _apply_option();
}

String.prototype.toHex = function () {
  var hash = 0;
  if (this.length === 0) return hash;
  for (var i = 0; i < this.length; i++) {
    hash = this.charCodeAt(i) + ((hash << 5) - hash);
    hash = hash & hash;
  }
  var color = "#";
  for (var i = 0; i < 3; i++) {
    var value = (hash >> (i * 8)) & 255;
    color += ("00" + value.toString(16)).substr(-2);
  }
  return color;
};
