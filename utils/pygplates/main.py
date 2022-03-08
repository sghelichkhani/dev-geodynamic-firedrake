import pygplates 
import numpy as np

print('Imported pyGPlates version: %s' % pygplates.Version.get_imported_version())

gpmlz_fi    = '/Applications/GPlates-2.2.0/SampleData/FeatureCollections/Hotspots/Hotspots_Compilation_Whittaker_etal.gpmlz'
rot_fi      = '/Applications/GPlates-2.2.0/SampleData/FeatureCollections/Rotations/Matthews_etal_GPC_2016_410-0Ma_GK07.rot'
spoly_fi    = '/Applications/GPlates-2.2.0/SampleData/FeatureCollections/StaticPolygons/Muller_etal_AREPS_2016_StaticPolygons.gpmlz'
out_fi      = 'test.xy'

# the age 
age = 10

# generating time
lons = np.linspace(-180, 180, 361)
lats = np.linspace(-90, 90, 181)
lons_x, lats_x = np.meshgrid(lons, lats)
lons_v, lats_v = lons_x.reshape(np.prod(np.shape(lons_x))), lats_x.reshape(np.prod(np.shape(lats_x)))

# loading sample feature collection
features = pygplates.FeatureCollection(gpmlz_fi)
# rotation model given in rot_fi
rotation_model = pygplates.RotationModel(rot_fi)

# generating sample points a sphere
all_points = pygplates.MultiPointOnSphere(np.column_stack((lats_v, lons_v)))
a_point = pygplates.PointOnSphere((lats_v[0], lons_v[0]))

# Reconstructing features
ReconstructedFeatures = []
pygplates.reconstruct(features, rotation_model, ReconstructedFeatures, reconstruction_time=age)

assigned_point_features = pygplates.partition_into_plates(
        spoly_fi,
        rotation_model,
        features,#point_features,
        properties_to_copy = [
                    pygplates.PartitionProperty.reconstruction_plate_id,
                    pygplates.PartitionProperty.valid_time_period])

#equivalent_stage_rotation = rotation_model.get_rotation(10, moving_plate_id=101, from_time=10, fixed_plate_id=101)
# Get the rotation from 11Ma to 10Ma, and the feature's reconstruction plate ID.
equivalent_stage_rotation = rotation_model.get_rotation(10, 1)


# Calculate a velocity for each reconstructed point over the 1My time interval.
velocities = pygplates.calculate_velocities(
        all_points.get_points(),
        equivalent_stage_rotation,
        1,
        pygplates.VelocityUnits.cms_per_yr)



