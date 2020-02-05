
#############################################
#' @title Overwrite facet_wrap stripes layout to facet_grid stripes layout at ggplot2 object.
#' 
#' @description 
#' Function overwritting facet_wrap stripes layout by facet_grid stripes layout at ggplot2 object.
#' This is to overcome the lack of full scale = free functionality at the current implementation of facet_grid().
#' facet_wrap() has this functionality in place already.
#'
#' @param p - ggplot2 object to modify.
#' @param facet_vars_rows - variable name(s) in data which should be used for defining faceting groups on the rows,
# by default this is the first variable from facets argument supplied to facet_wrap() function at object p.
#' @param panel_resize_factor Resize (height) upper panel by panel_resize_factor.
#'
#' @return gtable class object
#'
#' @export
#' @examples
#' 
#' df <- data.frame(facet = c(rep("A.Top", 10), rep("B.Bottom", 5)),
#'                  series = c(rep("rNa", 5), rep("coma", 5), rep("takeall_pct", 5)), 
#'                  sample_fraction = rep(c(0.15, 0.25, 0.32, 0.39, 0.41), 3),
#'                  value = c(1, 2, 1, 3, 3, 250, 270, 290, 304, 308, 40, 15, 18, 30, 90), 
#'                  stringsAsFactors = FALSE)
#' df <- rbind(cbind(col = "col1", df), cbind(col = "col2", df))#, cbind(col = "col3", df))
#' df[21:25, 5] <- df[21:25, 5] * 10 
#' 
#' p <- ggplot(data = df, mapping = aes(x = sample_fraction, y = value)) +
#'   geom_line(data = subset(df, facet == "A.Top"), mapping = aes(linetype = series)) +
#'   geom_point(data = subset(df, facet == "A.Top"), mapping = aes(shape = series), size = 2) +
#'   geom_bar(data = subset(df, facet == "B.Bottom"), mapping = aes(y = value), stat = "identity") +
#'   facet_wrap(~ facet + col, nrow = 2, scales = "free_y")
#'   
#' plot(p)
#' 
#' gt <- to_facet_grid_layout(p)
#' grid::grid.draw(gt)
#' 
to_facet_grid_layout <- function(p, facet_vars_rows = NULL, panel_resize_factor = 5) {
  
  # get faceting variable names from plot object p
  facet_vars <- names(p$facet$params$facets)
  
  if(is.null(facet_vars_rows))
    facet_vars_rows <- facet_vars[1] # set of variables defining faceting groups on the rows
  
  facet_vars_cols <- setdiff(facet_vars, facet_vars_rows)
  if(length(facet_vars_cols) == 0)
    facet_vars_cols <- "."
  
  # get theme strip settings from plot object p
  p_theme_strip_idx <- grep('strip.', names(p$theme))
  p_theme_strip <- p$theme[p_theme_strip_idx]
  if(length(p_theme_strip) != 0)
    attributes(p_theme_strip) <- c(list(names = names(p_theme_strip)), attributes(p$theme)[-1])
  
  # remove theme y-label and strips from plot object p, bind modified object to p1 variable
  p1 <- p + ylab(NULL)
  p1$theme[p_theme_strip_idx] <- NULL
  p1 <- p1 + theme(strip.text = element_blank())
  
  # create new plot object based on facet_grid 
  p2 <- ggplot(data = p1$data, mapping = p1$mapping) +
    facet_grid(rows = as.formula(paste(facet_vars_rows, "~", facet_vars_cols)),
               scales = "free_y", switch = "y") + 
    p_theme_strip + # preserve plot p theme strip settings on new grid
    theme(strip.placement = "outside") # except strip.placement, if exists
  
  # preserve labeller settings from plot object p
  p2$facet$params$labeller <- p$facet$params$labeller # labeller mapping
  
  #### gtable object modification

  # get gtable objects from original plot p1 and new plot p2
  gt1 <- ggplot_gtable(ggplot_build(p1))
  gt2 <- ggplot_gtable(ggplot_build(p2))
  
  # add top stripes to plot p1 from plot p2
  gt1 <- gtable::gtable_add_rows(gt1, heights = unit(1, 'cm'), pos = 4)
  gt1 <- gtable::gtable_add_grob(gt1, grobs = gt2$grobs[grep('strip-t', gt2$layout$name)], 
                                 t = 5, l = gt1$layout[grep('strip-t.+1$', gt1$layout$name),]$l)

  # add left stripes to plot p1 from plot p2
  gt1 <- gtable::gtable_add_cols(gt1, widths = unit(1, 'cm'), pos = 2)
  panel_1_1_l <- gt1$layout[grep('panel-1-1$', gt1$layout$name),]$l
  gt1 <- gtable::gtable_add_grob(gt1, gt2$grobs[grep('strip-l', gt2$layout$name)], 
                                 t = gt1$layout[grepl('panel-.*$', gt1$layout$name) & gt1$layout$l == panel_1_1_l, ]$t, l = 3)
  
  # resize top/bottom panels
  #gtable::gtable_show_layout(gt1)
  panel1_id <- gt1$layout[grep('panel-1-1$', gt1$layout$name),]$t
  gt1$heights[panel1_id] <- panel_resize_factor * gt1$heights[panel1_id]
  
  #grid::grid.draw(gt1)
  return(gt1)
  
}

