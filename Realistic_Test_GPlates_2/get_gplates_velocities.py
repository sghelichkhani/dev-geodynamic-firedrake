from scipy.io.netcdf import netcdf_file
import scipy.interpolate, numpy

# Constants:
gplates_dir             = "gplates_data/nc_files_muller_2019/"
plate_scaling_factor    = 1.0 # Factor to scale plate velocities to RMS velocity of model (RMS_Model / RMS_Earth)
velocity_non_dim_factor = 2890e3/1.0e-6 # Factor to non-dimensionalise gplates velocities: d/kappa
time_dim_factor         = 2890e3**2/1.0e-6 # Factor to dimensionalise model time: d^2/kappa
deg2rad                 = numpy.pi/180.
gplates_stage_time      = 31536000000000.0 * (1./plate_scaling_factor)    # 1 Myr in seconds (scaled)
total_stages            = 250
start_time              = 0.0

# Coordinate transformation function:
def xyz2spherical(X):
  # Returns r, phi, theta for a x,y,z coordinate. Note that theta is a weird latitude that starts at 0
  # at the south pole and is pi at the north pole. Phi is a weird longitude, which starts at 0 on the date line.
  r = numpy.sqrt(X[:,0]**2+X[:,1]**2+X[:,2]**2)
  phi = numpy.arctan2(X[:,1], X[:,0]) + numpy.pi # Longitude
  theta = numpy.pi-numpy.arccos(X[:,2]/r) # Latitude
  return r, phi, theta

class Gplates_Interpolator(object):

  def __init__(self):
    self.current_time = None
    self.stage_time   = None
    self.stage        = None
    self.lats         = None
    self.lons         = None
    self.coords       = None
    self.set_time(start_time)
  
  def set_time(self, t):
    if (self.current_time == t * time_dim_factor):
      return
    self.current_time = t * time_dim_factor
    self.stage_time   = self.current_time % gplates_stage_time
    print ( "GPlates: Stage time (s):",self.stage_time,"of:",gplates_stage_time, flush=True)

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

    nc = netcdf_file(gplates_dir+'velocity_%s.00Ma.nc' % str(self.stage), 'r')
  
    # As lats and lons do not change between files, only store these once:
    if self.lats is None:
      self.lats = (nc.variables['lat'][1:-1] + 90.)*deg2rad # do not include poles
    if self.lons is None:
      self.lons = (nc.variables['lon'][:-1] + 180.)*deg2rad # do not include 0 meridian twice
    if self.coords is None:
      self.coords = numpy.array([(lat,lon) for lat in self.lats for lon in self.lons])

    # Set up interpolators:
    self.intpx  = scipy.interpolate.RectBivariateSpline(self.lats, self.lons, nc.variables['velocity'][0,:-1,1:-1].T,kx=1,ky=1)
    self.intpy  = scipy.interpolate.RectBivariateSpline(self.lats, self.lons, nc.variables['velocity'][1,:-1,1:-1].T,kx=1,ky=1)
    self.intpz  = scipy.interpolate.RectBivariateSpline(self.lats, self.lons, nc.variables['velocity'][2,:-1,1:-1].T,kx=1,ky=1)

  # Nodal velocity vector:
  def get_velocities(self,X):
    r, phi, theta = xyz2spherical(X)   

    return numpy.column_stack((self.intpx(theta,phi, grid=False) * velocity_non_dim_factor * plate_scaling_factor,
                    self.intpy(theta,phi, grid=False) * velocity_non_dim_factor * plate_scaling_factor, 
                    self.intpz(theta,phi, grid=False) * velocity_non_dim_factor * plate_scaling_factor))

# Set up interpolation class:
gplates_interpolator=Gplates_Interpolator()




