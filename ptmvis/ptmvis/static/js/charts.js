function getOverviewOption(_data) {
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
