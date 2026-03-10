#' Forest Style Dot Plot
#'
#' Creates a forest-plot style visualization for event fractions with confidence intervals.
#'
#' @param data A dataframe containing the required columns:
#'   \itemize{
#'     \item \code{event_name}: Categorical (Y-axis)
#'     \item \code{fraction_value}: Numeric 0-1 (X-axis)
#'     \item \code{ci_min}: Numeric (Lower error bound)
#'     \item \code{ci_max}: Numeric (Upper error bound)
#'     \item \code{count_size}: Numeric (Point size scaling)
#'     \item \code{is_significant}: Logical (Coloring)
#'   }
#' @param title Optional string for plot title.
#' @param xlab Optional string for X-axis label.
#'
#' @return A ggplot object.
#' @export
plot_forest_dot <- function(data, title = "Event Fraction Forest Plot", xlab = "Fraction Value") {
  
  # Ensure data is ordered by the input row order (bottom to top usually for plots)
  # We convert event_name to factor with specific levels to preserve order
  data$event_name <- factor(data$event_name, levels = unique(data$event_name))
  
  # Prepare grid line data
  n_events <- length(levels(data$event_name))
  grid_data <- data.frame(y_pos = seq(0.5, n_events + 0.5, 1))
  
  p <- ggplot(data, aes(x = fraction_value, y = event_name, 
                        text = paste("Event:", event_name, "<br>Fraction:", round(fraction_value,3), "<br>Count:", count_size))) +
    # 1. Base Layer: Reference Line at 0.5 and Connecting Segments
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
    geom_line(aes(group = event_name), color = "grey60", linetype = "dotted", size = 0.5) +
    
    # 2. Geometries: Error Bars and Points with Dodging
    geom_errorbarh(aes(xmin = ci_min, xmax = ci_max, color = Group), 
                   height = 0.2, position = position_dodge(width = 0.6)) +
    geom_point(aes(size = count_size, color = Group), 
                   position = position_dodge(width = 0.6)) +
    
    # 3. Aesthetics
    scale_color_manual(values = c("Case" = "#F8766D", "Control" = "#00BFC4"), 
                       name = "Group") +
    scale_size_continuous(name = "Read Count") +
    scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
    coord_cartesian(xlim = c(-0.02, 1.02), clip = "off") +
    
    # 4. Styling
    theme_bw() +
    theme(
      panel.border = element_rect(colour = "grey", fill = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      axis.text.y = element_text(size = 10),
      plot.margin = margin(t=10, r=100, b=10, l=10, unit = "pt")
    ) +
    
    # 5. Separation Lines (Restricted to 0-1 range)
    geom_segment(data = grid_data, 
                 aes(y = y_pos, yend = y_pos, x = 0, xend = 1), 
                 inherit.aes = FALSE,
                 color = "grey90", size = 0.5) +
               
    # 6. Text Annotations (DeltaPSI and FDR)
    geom_text(data = subset(data, Group == "Case"),
              aes(x = 1.06, y = event_name, label = sprintf("dPSI: %.2f  FDR: %.1e", dpsi, fdr)),
              hjust = 0, size = 3, color = "black", inherit.aes = FALSE) +

    # 7. Labels
    labs(title = title, x = xlab, y = NULL)
    
  return(p)
}

#' Information Tooltip Helper
#'
#' @param content The text to display in the tooltip.
#' @param placement Where to place the tooltip (top, bottom, left, right).
#'
#' @return A shiny tag for an info icon with tooltip.
#' @export
info_icon <- function(content, placement = "right") {
  tags$span(
    style = "cursor: help; color: grey; margin-left: 5px;",
    `data-toggle` = "tooltip",
    `data-placement` = placement,
    `title` = content,
    icon("info-circle")
  )
}
