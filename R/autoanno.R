
#' Automatic Annotation of Cell Clusters
#'
#' This function automatically annotates cell clusters based on the sum of
#' avg.exp.scaled values for each feature.groups in a DotPlot.
#'
#' @param sce_object A Seurat object containing the single-cell data.
#' @param marker_file Path to the marker genes file (e.g., CSV or TXT).
#' @param cols Color scheme for DotPlot (default is RdYlBu).
#' @return A data frame containing cluster IDs and their corresponding annotations.
#' @export
#' @examples
#' # Example usage:
#' # annotations_df <- autoanno(sce.all.filt, "my_markers.txt")
autoanno <- function(sce_object, marker_file, cols = "RdYlBu") {
  library(Seurat)
  library(data.table)
  library(ggplot2)
  # 读取标记基因文件
  a <- read.table(marker_file, sep = ",")
  gt <- split(a[,2], a[,1])
  # 创建DotPlot
  p <- DotPlot(sce_object, features = gt, cols = cols)
  # 提取DotPlot数据
  dot_data <- as.data.table(p$data)
  # 按照簇 (id) 和细胞类型 (feature.groups) 进行汇总
  summed_data <- dot_data[, list(avg_exp_scaled_sum = sum(avg.exp.scaled)), by = .(id, feature.groups)]
  # 找到每个簇中累加值最大的feature.groups
  annotations <- summed_data[, .SD[which.max(avg_exp_scaled_sum)], by = id]
  # 将注释结果转换为数据框
  annotations_df <- data.frame(Cluster = annotations$id, Annotation = annotations$feature.groups)
  # 返回注释数据框
  return(annotations_df)
}




