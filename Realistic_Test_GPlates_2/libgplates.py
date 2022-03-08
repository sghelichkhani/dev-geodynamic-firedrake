import itertools, os.path, pygplates, numpy

# get the naming right 
data_root = '/Applications/GPlates-2.2.0/SampleData/FeatureCollections/'

topology_filenames = []
topology_filenames.append(os.path.join(data_root,'DynamicPolygons/Matthews_etal_GPC_2016_MesozoicCenozoic_PlateTopologies.gpmlz'))
topology_filenames.append(os.path.join(data_root,'DynamicPolygons/Matthews_etal_GPC_2016_Paleozoic_PlateTopologies.gpmlz'))

rotation_filenames = []
rotation_filenames.append(os.path.join(data_root, 'Rotations/Matthews_etal_GPC_2016_410-0Ma_GK07.rot'))

# factor to scale plate velocities to RMS velocity of model (RMS_Model / RMS_Earth)
# this is specially used in low-Rayleigh-number simulations 
plate_scaling_factor    = 1.0 
# Factor to non-dimensionalise gplates velocities: d/kappa
# non-dimensionalise the retrieved velocity from pygplates
velocity_non_dim_factor = 2890e3/1.0e-6 
# Factor to dimensionalise model time: d^2/kappa
# dimensionalising time to secibds 
time_dim_factor         = 2890e3**2/1.0e-6 
# 1 Myr in seconds (scaled)
gplates_stage_time      = 31536000000000.0 * (1./plate_scaling_factor)    
# what is the geologic zero, where do we want to start from
geologic_zero           = 250
start_time              = 0.0

# Make sure rotation & topology files exist
if False in [os.path.isfile(fi) for fi in topology_filenames+rotation_filenames]:
    raise ValueError("Unable to find:", list(itertools.compress(topology_filenames+rotation_filenames,\
                                    [[not os.path.isfile(fi) for fi in topology_filenames+rotation_filenames]])))

class pygplate_runner(object):
    def __init__(self, rotation_filenames, topology_filenames, delta_time=1.):
        """
        A class for getting surface velocities from pygplates
            rotation_filenames:
            topology_filenames: 
            delta_time:
        """
        
        # Rotation model(s)
        self.rotation_model = pygplates.RotationModel(rotation_filenames)
         
        # Topological plate polygon feature(s).
        self.topology_features = []
        for fname in topology_filenames:
            for f in pygplates.FeatureCollection(fname):
                self.topology_features.append(f)

        # time window for velocity interpolations
        self.delta_time = delta_timei

        self.model_time = None

    def set_time(self, t):
        
        # model_time: dimensionalised model time assuming that t starts from zero
        model_time = t * time_dim_factor

        geologic_time = geologic_zero - model_time/gplates_state_time 
        print( "GPlates: Stage time (s):",self.stage_time,"of:",gplates_stage_time, flush=True)

        # Opens the GPlates data file appropriate for the current simulation time (stage):
        total_stage_time = gplates_stage_time * (1. / plate_scaling_factor)
        new_stage        = int(numpy.ceil(total_stages - (self.current_time / gplates_stage_time) ))

        if (self.stage == new_stage):
          return
        self.stage = new_stage

        # If appropriate, read GPlates file:
        if(self.stage != total_stages):
          print( "GPlates *** Plate Stage Has Changed *** ",flush=True )
        print( "GPlates: Current dimensional time (s):",self.current_time,flush=True )
        print( "GPlates *** Reading from file: velocity_%s.00Ma.nc *** " % str(self.stage), flush=True )


    def compute_velocities(self, coords, reconstruction_time):
        
        # make sure the coordinates don't change in between (in case we switch to 
        # a temporally varying mesh
        if not self.coords:
            self.coords = coords
        else:
            reset_GPML_velocity_feature = False in self.coords == coords

        # Make sure velocity_domain_features are initiated before passing things on
        if not self.velocity_domain_features or reset_GPML_velocity_feature:
            self.velocity_domain_features = self.__make_GPML_velocity_feature(coords)
            self.coords = coords
        
        # calculate velocities here 
        all_velocities = self.__calc_velocities(velocity_domain_features=self.velocity_domain_features,
                                      topology_features=self.topology_features,
                                      rotation_model=self.rotation_model,
                                      time=reconstruction_time,
                                      delta_time=self.delta_time)
        self.boundary_velocity =  numpy.array([i.to_xyz() for i in all_velocities])

    def __make_GPML_velocity_feature(self, coords):
        """ function to make a velocity mesh nodes at an arbitrary set of points defined in 
             coords[# of points, 3] = x, y, z"""
   
        # Add points to a multipoint geometry
        multi_point = pygplates.MultiPointOnSphere(\
                  [pygplates.PointOnSphere(x=coords[i,0], y=coords[i,1], z=coords[i,2], normalise=True) for i in range(numpy.shape(coords)[0])])
    
        # Create a feature containing the multipoint feature, and defined as MeshNode type
        meshnode_feature = pygplates.Feature(pygplates.FeatureType.create_from_qualified_string('gpml:MeshNode'))
        meshnode_feature.set_geometry(multi_point)
        meshnode_feature.set_name('Velocity Mesh Nodes from pygplates')
    
        output_feature_collection = pygplates.FeatureCollection(meshnode_feature)
    
        # NB: at this point, the feature could be written to a file using
        # output_feature_collection.write('myfilename.gpmlz')
    
        # for use within the notebook, the velocity domain feature is returned from the function
        return output_feature_collection
    
    
    def __calc_velocities(self, velocity_domain_features, topology_features, rotation_model, time, delta_time):
        # All domain points and associated (magnitude, azimuth, inclination) velocities for the current time.
        all_domain_points = []
        all_velocities = []
    
        # Partition our velocity domain features into our topological plate polygons at the current 'time'.
        plate_partitioner = pygplates.PlatePartitioner(topology_features, rotation_model, time)
    
        for velocity_domain_feature in velocity_domain_features:
    
            # A velocity domain feature usually has a single geometry but we'll assume it can be any number.
            # Iterate over them all.
            for velocity_domain_geometry in velocity_domain_feature.get_geometries():
    
                for velocity_domain_point in velocity_domain_geometry.get_points():
    
                    all_domain_points.append(velocity_domain_point)
    
                    partitioning_plate = plate_partitioner.partition_point(velocity_domain_point)
                    if partitioning_plate:
    
                        # We need the newly assigned plate ID to get the equivalent stage rotation of that tectonic plate.
                        partitioning_plate_id = partitioning_plate.get_feature().get_reconstruction_plate_id()
    
                        # Get the stage rotation of partitioning plate from 'time + delta_time' to 'time'.
                        equivalent_stage_rotation = rotation_model.get_rotation(time, partitioning_plate_id, time + delta_time)
    
                        # Calculate velocity at the velocity domain point.
                        # This is from 'time + delta_time' to 'time' on the partitioning plate.
                        velocity_vectors = pygplates.calculate_velocities(
                            [velocity_domain_point],
                            equivalent_stage_rotation,
                            delta_time)
                        
                        # add it to the list
                        all_velocities.extend(velocity_vectors)
                    else:
                        print("Found error")
                        all_velocities.extend([pygplates.Vector3D(0,0,0)])
    
        return all_velocities


rec_model = pygplate_runner(rotation_filenames=rotation_filenames, topology_filenames=topology_filenames)

