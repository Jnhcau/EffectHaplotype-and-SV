library(networkD3)
library(tidyverse)
library(readxl)
library(htmlwidgets)
library(webshot)
library(cols4all) 
library(readxl)
library(ggalluvial)

df <- read_xlsx("桑吉图.xlsx")
sample_sizes <- c("6070" = 30, "8090" = 85, "0010" = 53, "pubus" = 74, "expvp" = 89)


df_gene <- df %>%
  filter(genename == "knr6") %>%
  mutate(group = as.character(group),
         percent = value / sample_sizes[group])

nodes <- data.frame(name = unique(c(df_gene$type, df_gene$group)),
                    stringsAsFactors = FALSE)
nodes$ID <- 0:(nrow(nodes) - 1)

links <- data.frame(
  source = match(df_gene$type, nodes$name) - 1,
  target = match(df_gene$group, nodes$name) - 1,
  value  = df_gene$percent
)

fixed_colors <- c(
  "Hap1" = "#D55E00", "Hap2" = "#0072B2", "Hap3" = "#009E73",
  "Hap4" = "#CC79A7", "Hap5" = "#E69F00", "6070" = "#56B4E9",
  "8090" = "#E42320", "0010" = "#008B8B", "pubus" = "#6A8EC9",
  "expvp" = "#652884"
)
node_names <- nodes$name

node_colors <- fixed_colors[node_names]
my_colors <- sprintf(
  'd3.scaleOrdinal().domain(%s).range(%s)',
  jsonlite::toJSON(node_names, auto_unbox = TRUE),
  jsonlite::toJSON(unname(node_colors), auto_unbox = TRUE)
)

# 创建桑基图
sankey_plot <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  fontSize = 14,
  fontFamily = "Times New Roman",
  nodeWidth = 20,
  nodePadding = 10,
  height = 400,
  width = 400,
  colourScale = JS(my_colors),
  iterations = 0
)
p1 <- htmlwidgets::onRender(sankey_plot, '
  function(el, x) {
    var svg = d3.select(el).select("svg");

    function createValidID(name) {
      if (!name) return "unknown";
      return name.replace(/[^a-zA-Z0-9-]/g, "_");
    }

    // 清除 NA 节点
    var nodesToRemove = [];
    svg.selectAll(".node").each(function(d) {
      if (!d.name || d.name === "NA") nodesToRemove.push(d.name);
    });

    nodesToRemove.forEach(function(nodeName) {
      svg.selectAll(".link").filter(function(link) {
        return link.target.name === nodeName;
      }).remove();

      svg.selectAll(".node").filter(function(node) {
        return node.name === nodeName;
      }).remove();
    });

    // 添加渐变定义并应用到 link 上
    svg.selectAll(".link").each(function(d) {
      var gradientID = "gradient-" + createValidID(d.source.name) + "-" + createValidID(d.target.name);
      var gradient = svg.append("defs")
        .append("linearGradient")
        .attr("id", gradientID)
        .attr("gradientUnits", "userSpaceOnUse")
        .attr("x1", d.source.x + d.source.dx / 2)
        .attr("y1", d.source.y + d.source.dy / 2)
        .attr("x2", d.target.x + d.target.dx / 2)
        .attr("y2", d.target.y + d.target.dy / 2);

      var sourceColor = d3.select(el).selectAll(".node")
        .filter(function(node) { return node.name === d.source.name; })
        .select("rect").style("fill");

      var targetColor = d3.select(el).selectAll(".node")
        .filter(function(node) { return node.name === d.target.name; })
        .select("rect").style("fill");

      gradient.append("stop")
        .attr("offset", "0%")
        .attr("stop-color", sourceColor);

      gradient.append("stop")
        .attr("offset", "100%")
        .attr("stop-color", targetColor);

      d3.select(this).style("stroke", "url(#" + gradientID + ")");
    });

    // 给节点加黑色边框
    svg.selectAll(".node rect")
      .style("stroke", "black")
      .style("stroke-width", "1.5px");
  }
')

# 显示图形
p1

library(htmlwidgets)
saveWidget(p1, "knr6.html", selfcontained = TRUE)

library(pagedown)
htmlwidgets::saveWidget(p1, "knr6.html", selfcontained = TRUE)

# 用 Chrome 打印为 PDF
pagedown::chrome_print("D:\\毕业\\CAU毕业\\发表文章\\gchap_SV\\FIG\\FIG3\\knr6.html", output = "D:\\毕业\\CAU毕业\\发表文章\\gchap_SV\\FIG\\FIG3\\knr6.pdf")


