function getOverviewOption(modifications, coOccurrenceData, classCounts) {
  var option = {
    title: [
      {
        text:
          "Total Modifications " +
          Object.values(classCounts).reduce((a, b) => a + b, 0) +
          " | Distinct Modifications " +
          Object.keys(modifications).length +
          " | Modification Classes " +
          Object.keys(classCounts).length,
        top: "top",
        left: "left",
        textStyle: {
          fontSize: 14,
          fontWeight: "lighter",
        },
      },
    ],
    grid: [
      {
        top: 70,
        bottom: "12%",
        left: "10%",
        height: "auto",
        width: "25%",
        show: true,
      },
      {
        top: 70,
        bottom: "12%",
        left: "37%",
        height: "auto",
        width: "14%",
        show: true,
      },
      {
        top: 70,
        bottom: "12%",
        left: "53%",
        height: "auto",
        width: "14%",
        show: true,
      },
      {
        top: 70,
        bottom: "30%",
        right: "5%",
        height: "auto",
        width: "22%",
        show: true,
      },
    ],
    xAxis: [
      {
        gridIndex: 0,
        type: "category",
        name: "Modification",
        nameLocation: "center",
        nameTextStyle: { fontWeight: "lighter", fontSize: 11 },
        nameGap: 30,
        data: Object.keys(modifications),
        axisTick: {
          alignWithLabel: true,
          interval: 0,
          length: 3,
        },
        axisLabel: {
          show: true,
          formatter: (i) => {
            if (modifications[i]["display_name"].length > 5) {
              return modifications[i]["display_name"].substring(0, 5) + "...";
            } else {
              return modifications[i]["display_name"];
            }
          },
          fontWeight: "lighter",
          fontSize: 9,
          rotate: 45,
        },
        axisPointer: {
          show: true,
          label: {
            show: false,
          },
          triggerEmphasis: false,
          triggerTooltip: false,
        },
      },
      {
        gridIndex: 1,
        type: "value",
        name: "Mass Shift [Da]",
        nameLocation: "center",
        nameTextStyle: { fontWeight: "lighter", fontSize: 11 },
        nameGap: 30,
        axisLabel: {
          show: true,
          interval: 0,
          fontWeight: "lighter",
          fontSize: 11,
          rotate: 45,
        },
      },
      {
        gridIndex: 2,
        type: "value",
        name: "Count",
        nameLocation: "center",
        nameTextStyle: { fontWeight: "lighter", fontSize: 11 },
        nameGap: 30,
        axisLabel: {
          show: true,
          interval: 0,
          fontWeight: "lighter",
          fontSize: 11,
          rotate: 45,
        },
      },
      {
        gridIndex: 3,
        type: "category",
        name: "Modification Class",
        nameLocation: "center",
        nameTextStyle: { fontSize: 11 },
        nameGap: 90,
        data: Object.keys(classCounts),
        axisTick: {
          show: false,
        },
        axisLabel: {
          show: true,
          interval: 0,
          fontWeight: "lighter",
          fontSize: 11,
          rotate: 45,
        },

        axisPointer: {
          show: true,
          label: {
            show: false,
          },
          triggerEmphasis: false,
        },
      },
    ],
    yAxis: [
      {
        gridIndex: 0,
        type: "category",
        name: "Modification",
        nameLocation: "center",
        nameTextStyle: { fontSize: 11 },
        nameGap: 140,
        data: Object.keys(modifications),
        inverse: true,
        axisTick: {
          alignWithLabel: true,
          interval: 0,
          length: 1,
        },
        axisLabel: {
          show: true,
          formatter: (i) => {
            if (modifications[i]["display_name"].length > 20) {
              return modifications[i]["display_name"].substring(0, 20) + "...";
            } else {
              return modifications[i]["display_name"];
            }
          },
          fontWeight: "lighter",
          fontSize: 11,
        },
        axisPointer: {
          show: true,
          label: {
            show: false,
          },
          triggerEmphasis: false,
          triggerTooltip: false,
        },
      },
      {
        gridIndex: 1,
        type: "category",
        data: Object.keys(modifications),
        show: false,
        inverse: true,
        axisPointer: {
          show: true,
          label: {
            show: false,
          },
          triggerEmphasis: false,
        },
      },
      {
        gridIndex: 2,
        type: "category",
        data: Object.keys(modifications),
        show: false,
        inverse: true,
        axisPointer: {
          show: true,
          label: {
            show: false,
          },
          triggerEmphasis: false,
        },
      },
      {
        gridIndex: 3,
        type: "value",
        name: "Count",
        nameLocation: "center",
        nameTextStyle: { fontSize: 11 },
        nameGap: 40,
        axisLabel: {
          show: true,
          fontWeight: "lighter",
          fontSize: 11,
        },
      },
    ],
    tooltip: {
      backgroundColor: "#fbfbfbe6",
      borderColor: "#fbfbfb",
      textStyle: {
        color: "#111111",
        fontSize: 13,
      },
      formatter: (params) => {
        if (Array.isArray(params)) params = params[0];
        if (params.seriesIndex == 0)
          return (
            "Modification <code>" +
            modifications[params.data[1]].display_name +
            "</code> (" +
            modifications[params.data[1]].mass_shift.toFixed(2) +
            " Da) and <code>" +
            modifications[params.data[0]].display_name +
            "</code> (" +
            modifications[params.data[0]].mass_shift.toFixed(2) +
            " Da) have <code>" +
            params.data[2] +
            "</code> co-occurrences."
          );
        if (params.seriesIndex == 1)
          return (
            "Modification <code>" +
            modifications[params.name].display_name +
            "</code> assigned mass shift is <code>" +
            params.data +
            "</code> Da."
          );
        if (params.seriesIndex == 2)
          return (
            "Modification <code>" +
            modifications[params.name].display_name +
            "</code> counted <code>" +
            params.data +
            "</code> times."
          );
        if (params.seriesIndex == 3)
          return (
            "Modification class <code>" +
            params.name +
            "</code> assigned <code>" +
            params.data +
            "</code> times."
          );
      },
    },
    visualMap: [
      {
        type: "continuous",
        seriesIndex: [0],
        inRange: {
          color: [
            "#dddddd",
            "#cccccc",
            "#888888",
            "#666666",
            "#444444",
            "#222222",
            "#000000",
          ],
        },
        outOfRange: {
          color: ["#444444"],
        },
        min: 1,
        max: Math.max(...coOccurrenceData.map((_) => _[2])),
        orient: "horizontal",
        top: 50,
        left: "10%",
        itemHeight: 120,
        itemWidth: 12,
        text: [
          Math.max(...coOccurrenceData.map((_) => _[2])),
          "No. co-occurrences 1",
        ],
        textStyle: { fontWeight: "lighter", fontSize: 11 },
      },
    ],
    dataZoom: [
      {
        type: "inside",
        yAxisIndex: [0, 1, 2],
        brushSelect: false,
        throttle: 0,
      },
      {
        type: "inside",
        xAxisIndex: [0],
        throttle: 0,
      },
    ],
    series: [
      {
        type: "heatmap",
        xAxisIndex: 0,
        yAxisIndex: 0,
        data: coOccurrenceData,
        cursor: "default",
        large: true,
        animation: false,
        itemStyle: {
          borderWidth: 0.1,
          borderRadius: 2,
          borderColor: "#fbfbfb",
        },
      },
      {
        type: "bar",
        xAxisIndex: 1,
        yAxisIndex: 1,
        itemStyle: {
          color: "#111111",
        },
        barWidth: "50%",
        data: modifications.map((m) => m["mass_shift"]),
        cursor: "default",
      },
      {
        type: "bar",
        xAxisIndex: 2,
        yAxisIndex: 2,
        itemStyle: {
          color: "#111111",
        },
        barWidth: "50%",
        data: modifications.map((m) => m["count"]),
        cursor: "default",
      },
      {
        type: "bar",
        xAxisIndex: 3,
        yAxisIndex: 3,
        itemStyle: {
          color: "#111111",
        },
        barWidth: "50%",
        data: Object.values(classCounts),
        cursor: "default",
      },
    ],
  };
  return option;
}

function _getOverviewOption(_data) {
  // Filter data for top L modifications (wrt. frequency) and adjust layout settings.
  // (i) Init. list of node objects.
  nodes_data = Object.entries(_data[0]);
  nodes_data.sort((a, b) => {
    return b[1][0] - a[1][0];
  });
  nodes_values = Object.values(_data[0]).map((_) => {
    return _[0];
  });
  nodes_value_min = Math.min(...nodes_values);
  nodes_value_max = Math.max(...nodes_values);
  nodes = [];
  nodes_keys = new Set();
  L = 45;
  for (let entry of nodes_data.slice(0, L - 1)) {
    nodes.push({
      name: entry[0],
      value: entry[1][0],
      frequency: entry[1][1],
      symbolSize: 10,
    });
    nodes_keys.add(entry[0]);
  }
  // (ii) Init. list of link objects.
  links_data = Object.entries(_data[1]);
  links_value_max = Math.max(...Object.values(_data[1]));
  links = [];
  for (let entry of links_data) {
    let [s, t] = entry[0].split("@");
    if (!(nodes_keys.has(s) && nodes_keys.has(t))) {
      continue;
    }
    let v = entry[1];
    links.push({
      source: s,
      target: t,
      value: v,
      lineStyle: {
        width: v / links_value_max,
      },
    });
  }
  return {
    animation: false,
    tooltip: {
      backgroundColor: "rgba(51, 51, 51, 0.7)",
      borderColor: "transparent",
      textStyle: {
        fontWeight: "lighter",
        fontSize: 11,
        color: "#f0f5f5",
      },
      formatter: (params, ticket, callback) => {
        if (params.dataType == "edge") {
          return (
            `<b><u>` +
            params.data.source +
            `</u> and ` +
            `<u>` +
            params.data.target +
            `</u></b> have ` +
            params.data.value +
            ` co-occurrences.`
          );
        } else if (params.dataType == "node") {
          return (
            `<b>Modification <u>` +
            params.data.name +
            `</u></b> has ` +
            params.data.value +
            ` occurrences.<br>Occurs at least once in ` +
            params.data.frequency * 100 +
            `% of proteins in the dataset.`
          );
        }
      },
    },
    visualMap: {
      bottom: "bottom",
      left: "left",
      color: ["#111111", "#d1d1d1"],
      orient: "horizontal",
      precision: 0,
      itemWidth: 10,
      text: [nodes_value_max, "Modification Occurrence " + nodes_value_min],
    },
    series: [
      {
        name: "Modifications Graph",
        emphasis: {
          focus: "adjacency",
          label: { show: false },
          itemStyle: { shadowColor: "#dc5754", shadowBlur: 10 },
        },
        label: {
          show: true,
          color: "#333333",
          fontWeight: "lighter",
          fontSize: 11,
          backgroundColor: "rgba(240, 245, 245, 0.6)",
          borderRadius: 4,
        },
        edgeLabel: { show: false },
        type: "graph",
        layout: "circular",
        circular: {
          rotateLabel: true,
        },
        data: nodes,
        links: links,
        roam: true,
        lineStyle: { color: "#333333", curveness: 0.2 },
        itemStyle: {},
      },
    ],
  };
}
