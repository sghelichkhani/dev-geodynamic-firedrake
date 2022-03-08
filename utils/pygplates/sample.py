import matplotlib.pyplot as plt

import cartopy
import cartopy.crs as ccrs
import pygplates
import numpy as np

lons = np.linspace(-180, 180, 361)
lats = np.linspace(-90, 90, 181)
lons_x, lats_x = np.meshgrid(lons, lats)
lons_v, lats_v = lons_x.reshape(np.prod(np.shape(lons_x))), lats_x.reshape(np.prod(np.shape(lats_x)))

all_points = pygplates.MultiPointOnSphere(np.column_stack((lats_v, lons_v)))

rot_fi = '/Applications/GPlates-2.2.0/SampleData/FeatureCollections/Rotations/Matthews_etal_GPC_2016_410-0Ma_GK07.rot'
#dynpoly_fi = '/Applications/GPlates-2.2.0/SampleData/FeatureCollections/DynamicPolygons/Matthews_etal_GPC_2016_MesozoicCenozoic_PlateTopologies.gpmlz'
dynpoly_fi = '/Applications/GPlates-2.2.0/SampleData/FeatureCollections/DynamicPolygons/Matthews_etal_GPC_2016_MesozoicCenozoic_PlateTopologies.gpmlz'

# Load one or more rotation files into a rotation model.
rotation_model = pygplates.RotationModel(rot_fi)

# Load the topological plate polygon features.
topology_features = pygplates.FeatureCollection(dynpoly_fi)

# Load the features that contain the geometries we will calculate velocities at.
# These can be generated in GPlates via the menu 'Features > Generate Velocity Domain Points'.
#velocity_domain_features = pygplates.FeatureCollection('lat_lon_velocity_domain_9_18.gpml')
velocity_domain_features = pygplates.FeatureCollection(\
		pygplates.Feature.create_tectonic_section(feature_type=pygplates.FeatureType.gpml_unclassified_feature, geometry=all_points))
#velocity_domain_features = pygplates.FeatureCollection(pygplates.Feature.create_flowline(seed_geometry=all_points, times=[0, 5]))

# Calculate velocities using a delta time interval of 1My.
delta_time = 1

# Our geological times will be from 0Ma to 'num_time_steps' Ma (inclusive) in 1 My intervals.
num_time_steps = 140

# 'time' = 0, 1, 2, ... , 140
time = 100

print('Time: %d' % time)

# All domain points and associated (magnitude, azimuth, inclination) velocities for the current time.
all_domain_points = []
all_velocities = []
all_ids = [] 

# Partition our velocity domain features into our topological plate polygons at the current 'time'.
partitioned_domain_features = pygplates.partition_into_plates(
    topology_features,
    rotation_model,
    velocity_domain_features,
    reconstruction_time = time)

for partitioned_domain_feature in partitioned_domain_features:

    # We need the newly assigned plate ID to get the equivalent stage rotation of that tectonic plate.
    partitioning_plate_id = partitioned_domain_feature.get_reconstruction_plate_id()

    # Get the stage rotation of partitioning plate from 'time + delta_time' to 'time'.
    equivalent_stage_rotation = rotation_model.get_rotation(time, partitioning_plate_id, time + delta_time)

    # A velocity domain feature usually has a single geometry but we'll assume it can be any number.
    # Iterate over them all.
    for partitioned_domain_geometry in partitioned_domain_feature.get_geometries():

        partitioned_domain_points = partitioned_domain_geometry.get_points()
        all_ids.extend([partitioning_plate_id for i in partitioned_domain_points])
        all_domain_points.extend([i.to_lat_lon() for i in partitioned_domain_points])

        # Calculate velocities at the velocity domain geometry points.
        # This is from 'time + delta_time' to 'time' on the partitioning plate.
        partitioned_domain_velocity_vectors = pygplates.calculate_velocities(
            partitioned_domain_points,
            equivalent_stage_rotation,
            delta_time)

        # Convert global 3D velocity vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
        partitioned_domain_velocities = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(
                partitioned_domain_points,
                partitioned_domain_velocity_vectors)

        # Append results for the current geometry to the final results.
        all_velocities.extend([i.to_xyz() for i in partitioned_domain_velocities])
pnts = np.array(all_domain_points) 
vel = np.array(all_velocities) 
plt.close(1)
fig = plt.figure(num=1)
ax = fig.add_subplot(111)
ax.scatter(pnts[:,1], pnts[:,0], c=all_ids, s=0.1)
#ax.scatter(lons_v, lats_v, c=all_ids, s=0.1)
#ax.quiver(pnts[:,1], pnts[:,0], vel[:,1], vel[:,0], headlength=10, width=0.001)
#ax.quiver(pnts[:,1], pnts[:,0], vel[:,1], vel[:,0], headlength=10, width=0.001)
fig.show()
