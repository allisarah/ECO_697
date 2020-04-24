# ---- packages ----
{
  require(here)
  require(sf)
  require(maps)
  require(GISTools)
  require(ggplot2)
  require(cowplot)
  require(dismo)
  require(spData)
  
  source(here("data", "environment_vars.R"))
}

# ---- convex_hull  ----
{
  # Create some random points
  rbeta_inv = function(x1, x2, n) return(c(rbeta(n, x1, x2), rbeta(n, x2, x1)))
  # Create some random points  
  n = 10
  x = rnorm(n)
  y = rnorm(n)
  coord.tp = cbind(x, y)
  
  pts1 = st_multipoint(coord.tp)
  poly1 = st_convex_hull(pts1)
  
#buffer using st_buffer  
  poly_buff = st_buffer(poly1, 1)
  plot(poly_buff)
  plot(poly1, add = T)
  
  g_p = function(pts)
    return (geom_point(
      aes(x = X, y = Y), 
      data = data.frame(st_coordinates(pts)), 
      pch = 21, fill = "steelblue", cex = 4))
  g_h = function(hull)
    return (geom_polygon(
      aes(x = X, y = Y), 
      data = data.frame(st_coordinates(hull)), 
      fill = fill_hull, col = 1))
  
  n_pts1 = 7
  n_pts2 = 14
  set.seed(12343)
  
  t1 = theme(axis.title = element_blank())  
  gg = function() ggplot() + t1
  
  # CSR
  coords1 = matrix(runif(n_pts1 * 2), ncol = 2)
  
  # Funky
  coords2 = cbind(rbeta_inv(4, 1, n_pts2), rbeta_inv(10, 2, n_pts2))
  
  fill_hull = gray(0, 0.05)
  
  # Make sf objects
  pts1 = st_multipoint(coords1)
  pts2 = st_multipoint(coords2)
  
  # build convex hulls
  hull1 = st_convex_hull(pts1)  
  hull2 = st_convex_hull(pts2)  
  
  gg() + g_p(pts1)
  gg() + g_p(pts2)
  
  gg() + g_p(pts1) + g_h(hull1)
  gg() + g_p(pts2) + g_h(hull2)
  
  
  
  gg_hulls = plot_grid(
    gg() + g_p(pts1),
    gg() + g_p(pts2),
    gg() + g_p(pts1) + g_h(hull1),
    gg() + g_p(pts2) + g_h(hull2),
    nrow = 2)
}





# ---- voronoi ----
{
  n_pts = 15
  set.seed(12343)
  coords = data.frame(matrix(runif(n_pts * 2), ncol = 2))
  names(coords) = c("x", "y")
  polys_v = st_as_sf(voronoi(coords))
  polys_v$x = coords[, 1]
  polys_v$y = coords[, 2]
  
  gg_random_points = ggplot(polys_v) + 
    geom_sf(fill = "transparent", color = "transparent") + 
    geom_sf(fill = "transparent", data = st_union(polys_v)) + 
    geom_point(aes(x, y), data = polys_v) + 
    theme_minimal() +
    theme(axis.title = element_blank()) +
    ggtitle("CSR Points")
  
  gg_random_points_voronoi = ggplot(polys_v) + 
    geom_sf(fill = "transparent") + 
    geom_point(aes(x, y), data = polys_v) + 
    theme_minimal() +
    theme(axis.title = element_blank()) +
    ggtitle("CSR Points: Voronoi Diagram")
}


# ---- gis_operations ----
{
  us_crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
#base R
  plot(us_states)
 
    us_states = st_transform(us_states, us_crs)
    
    #ggplot
    ggplot(us_states) + geom_sf() 
    
  l_48_int = st_intersection(st_as_sf(us_states))
  l_48_union = st_union(st_as_sf(us_states))
  
  ggplot(l_48_union) + geom_sf() 
  us_states_buff_100k = lwgeom::st_make_valid(st_buffer(us_states, 1e5))
  
  gg_us_states = 
    ggplot(st_as_sf(us_states)) + 
    geom_sf(fill = "transparent") +
    ggtitle("State Polygons")
  
  gg_states_int = 
    ggplot(l_48_int) + 
    geom_sf(fill = "transparent") +
    ggtitle("Intersection")
  
  gg_states_union = 
    ggplot(l_48_union) + 
    geom_sf(fill = "transparent") +
    ggtitle("Union")
  
  gg_states_union
  
  gg_l_48_buf_100k = 
    ggplot(st_buffer(l_48_union, 1e5)) + 
    geom_sf(fill = "transparent")+
    geom_sf(data = l_48_union, col = "steelblue", fill = gray(0, 0.2)) +
    geom_sf(data = st_buffer(l_48_union, 1e6), col = "transparent", fill = "transparent") + 
    ggtitle("100 Kilometer Buffer")
  
  gg_l_48_buf_100k
  
  gg_l_48_buf_1000k = 
    ggplot(st_buffer(l_48_union, 1e6)) + 
    geom_sf(fill = "transparent") +
    geom_sf(data = l_48_union, col = "steelblue", fill = gray(0, 0.2)) +
    ggtitle("1000 Kilometer Buffer")
  
  gg_l_48_buf_1000k
  
  gg_states_buf_100k = 
    ggplot(st_buffer(us_states, 1e5)) + 
    geom_sf(alpha = 0.2, fill = "steelblue") +
    ggtitle("100 Kilometer Buffer")
  
}


#individual states buffered
states_buffered = st_buffer(us_states, 1e5)
plot(st_geometry(states_buffered))
#plot a subset of the polyf=gons eg - states 1 to 10
plot(states_buffered[1:10,])

plot(st_geometry(st_intersection(states_buffered[1:28,])))

#exclude problem polygons
plot(st_geometry(st_intersection(states_buffered[-29,])))

##fix problem polygons
require(lwgeom) ## st_make_valid works sometimes, but  didnt work here, can fix polygons in arcmap
plot(st_intersection(st_make_valid(states_buffered)))



# ---- us_voronoi ----
{
  us_crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  
  us_states = st_transform(us_states, us_crs)
  us_centroids = st_centroid(us_states)
  us_centroids
  plot(us_centroids, pch = 16)
  plot(st_geometry(us_centroids, pch = 16)) # st_geometry strips off all attributes except for coordinates
  
  us_voronoi = st_as_sf(voronoi(st_coordinates(us_centroids)))
  st_crs(us_voronoi) = st_crs(us_centroids)
  
  us_voronoi_2 = st_intersection(us_voronoi, st_union(st_as_sf(us_states)))
  
  gg_us_states_cen =
    ggplot(st_as_sf(us_states)) + 
    geom_sf(fill = "transparent") +
    geom_point(
      aes(x = X, y = Y), 
      data = data.frame(st_coordinates(us_centroids))) +
    ggtitle("Centroids of States")
  
  gg_us_vor =
    ggplot(us_voronoi_2) + 
    geom_sf(data = st_as_sf(us_states), col = "transparent", fill = "transparent")+
    geom_sf(fill = "transparent") +
    ggtitle("Voronoi Diagram: Centroids of States")
  
  gg_us_vor_cen =
    ggplot(us_voronoi_2) + 
    geom_sf(data = st_as_sf(us_states), col = "transparent", fill = "transparent")+
    geom_sf(fill = "transparent") +
    geom_point(
      aes(x = X, y = Y), 
      data = data.frame(st_coordinates(us_centroids))) +
    ggtitle("Voronoi Diagram: Centroids of States")
  
  print(gg_us_vor)
  print(gg_us_vor_cen)
}





# ----save_figures ----
{
  fig_width = 6
  fig_height = 3.5
  
  save_pdf = function(filename, grob, width = 1, height = 1)
  {
    pdf(file.path(tmp_img_dir, filename), width = width * fig_width, height = height * fig_height)
    print(grob)
    dev.off()
  }
  
  
  save_pdf("polygons_us_states.pdf", gg_us_states)
  save_pdf("polygons_us_states_union.pdf", gg_states_union)
  save_pdf("polygons_us_states_intersection.pdf", gg_states_int)
  
  save_pdf("polygons_lower_48_buffer_1000k.pdf", gg_l_48_buf_1000k)
  save_pdf("polygons_lower_48_buffer_100k.pdf", gg_l_48_buf_100k)
  save_pdf("polygons_us_states_buffer_100k.pdf", gg_states_buf_100k)
  
  save_pdf("polygons_us_states_centroids.pdf", gg_us_states_cen)
  save_pdf("polygons_us_states_voronoi.pdf", gg_us_vor_cen)
  
  save_pdf("voronoi_csr.pdf", plot_grid(gg_random_points, gg_random_points_voronoi, nrow = 1))
  
  save_pdf("convex_hulls.pdf", gg_hulls, 1.5, 1.5 * fig_width / fig_height)
}




library(mapproj)
library(ggmap)
library(DeducerSpatial)

require(maps)
require(ggmap)

map("world", "Ghana")
map.cities(country = "Ghana", capitals = 1)
map("region", "Greater Accra")

gama.map <- qmap(gama_utm)
get_map()
?register_google
has_google_account()
google_key()

register_google(key, account_type, client, signature, second_limit,
                day_limit, write = FALSE)
