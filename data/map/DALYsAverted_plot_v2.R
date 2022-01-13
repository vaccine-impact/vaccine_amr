# Create map of GAS burden averted (pharyngitis, impetigo, invasive disease,
# cellulitis, rheumatic heart disease by GAS vaccination

# load libraries
library (ggplot2)
library (sf)
library (rnaturalearth)
library (rnaturalearthdata)
library (rgeos)
library (data.table)
library (ggforce)
library (ggpubr)
library (stats)

# clear workspace
rm (list = ls ())

# ------------------------------------------------------------------------------
# create map of GAS burden averted (pharyngitis, impetigo, invasive disease,
# cellulitis, rheumatic heart disease by GAS vaccination
# ------------------------------------------------------------------------------
create_map <- function (daly_file, 
                        map_file) {
  
  # read results for dalys averted
  daly_table <- fread (file = daly_file)
  
  # estimate dalys averted per 1000
  daly_table [, dalysAvertedPer1000 := (RHD_DALYs_averted + 
                                          Invasive_infection_DALYs_averted +
                                          Impetigo_DALYs_averted + 
                                          Pharyngitis_DALYs_averted +  
                                          Cellulitis_DALYs_averted) / numVaxx * 1000]
  
  # map tutorial
  # https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html
  world <- ne_countries (scale       = "medium", 
                         returnclass = "sf")
  setDT (world)
  setkey (world, iso_a3)
  
  # combine tables to add geometry
  dt <- merge (x    = daly_table, 
               y    = world, 
               by.x = "Code", 
               by.y = "iso_a3", 
               all  = T )
  
  # generate map of dalys averted per 1000 vaccinated individuals
  dalys_plot <- ggplot (data = dt) +
    geom_sf (aes (fill = dalysAvertedPer1000, geometry = geometry)) + 
    scale_fill_viridis_c (option = "viridis", direction = -1, na.value = "grey90") +
    # scale_fill_viridis_c (option = "plasma", direction = -1, na.value = "grey90") +
    ggtitle ("DALYS averted per 1000 fully vaccinated individuals") + 
    theme (legend.title = element_blank()) + 
    theme (axis.text.x = element_blank(), axis.ticks = element_blank()) + 
    theme (axis.text.y = element_blank(), axis.ticks = element_blank()) + 
    theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme (plot.title = element_text(size = 10))
  
  # save plot
  ggsave (filename = map_file, 
          width = 6, 
          height = 3, 
          dpi = 600)
  
  return ()
  
} # end of function -- create_map
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# create map of daly results
create_map (daly_file = "DALYsAvertedbyCountry_allconditions.csv", 
            map_file  = "dalysAvertedPer1000.png")
# ------------------------------------------------------------------------------

