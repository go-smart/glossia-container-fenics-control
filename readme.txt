mesher_cgal --centre_radius 0.1  --output_prefix ire-cuboid-4  \                                                                                                                        ⏎
--zonefield 0.01 --farfield 0.01 --granularity 0.001 --nearfield 0.01  --output_gmsh --centre -0.0108 0.0174 -0.004 --zones tumour1.vtp:7 needle1.vtp:1:0.0003 needle2.vtp:2:0.0003 needle3.vtp:3:0.0003 needle4.vtp:4:0.0003 needle5.vtp:5:0.0003 needle6.vtp:6:0.0003 exterior.vtp:11 --vessels vessel1.vtp:8 vessel2.vtp:9 vessel3.vtp:10 --organ cuboid.vtp:0 --tissueid 0 --mark_zone_boundaries --output_vtk --bounding_radius 1
# WAS go-smart-mesher_cgal --centre_radius 0.0001  --output_prefix ire \                                                                                                                                 ⏎
#       --zonefield 0.001 --farfield 0.01 --granularity .0010 \
#       --output_gmsh --centre -0.0108 0.0174 -0.004 --zones tumour1.vtp:7 \
#       needle1.vtp:1 needle2.vtp:2 needle3.vtp:3 needle4.vtp:4 \
#       needle5.vtp:5 needle6.vtp:6 --vessels vessel1.vtp:8 vessel2.vtp:9 \
#       vessel3.vtp:10 --extent exterior.vtp:11 --tissueid 0 --mark_zone_boundaries
