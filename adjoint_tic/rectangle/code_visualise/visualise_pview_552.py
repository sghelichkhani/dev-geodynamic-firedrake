# trace generated using paraview version 5.5.2

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
import sys

#my_pvd_file = '/Users/sghelichkhani/Data/cosgrove/Workplace/dev-geodynamic-firedrake/adjoint_tic/rectangle/cases/SteepestDescent_Backtracking_StrongWolf/visual/opt_file.pvd'
my_pvd_file = sys.argv[1]

if not '.pvd' == my_pvd_file[-4:]:
    raise ValueError('Input Error: Given argument should be a .pvd file')

# create a new 'PVD Reader'
opt_filepvd = PVDReader(FileName=my_pvd_file)
opt_filepvd.PointArrays = ['InitTemperature', 'FinTemperature', 'Velocity', 'Pressure', 'Gradient']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1804, 1123]

# show data in view
opt_filepvdDisplay = Show(opt_filepvd, renderView1)

# get color transfer function/color map for 'InitTemperature'
initTemperatureLUT = GetColorTransferFunction('InitTemperature')
initTemperatureLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.25016039264684925, 0.865003, 0.865003, 0.865003, 0.5003207852936985, 0.705882, 0.0156863, 0.14902]
initTemperatureLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'InitTemperature'
initTemperaturePWF = GetOpacityTransferFunction('InitTemperature')
initTemperaturePWF.Points = [0.0, 0.0, 0.5, 0.0, 0.5003207852936985, 1.0, 0.5, 0.0]
initTemperaturePWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
opt_filepvdDisplay.Representation = 'Surface'
opt_filepvdDisplay.ColorArrayName = ['POINTS', 'InitTemperature']
opt_filepvdDisplay.LookupTable = initTemperatureLUT
opt_filepvdDisplay.OSPRayScaleArray = 'InitTemperature'
opt_filepvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
opt_filepvdDisplay.SelectOrientationVectors = 'Velocity'
opt_filepvdDisplay.ScaleFactor = 0.1
opt_filepvdDisplay.SelectScaleArray = 'InitTemperature'
opt_filepvdDisplay.GlyphType = 'Arrow'
opt_filepvdDisplay.GlyphTableIndexArray = 'InitTemperature'
opt_filepvdDisplay.GaussianRadius = 0.005
opt_filepvdDisplay.SetScaleArray = ['POINTS', 'InitTemperature']
opt_filepvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
opt_filepvdDisplay.OpacityArray = ['POINTS', 'InitTemperature']
opt_filepvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
opt_filepvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
opt_filepvdDisplay.SelectionCellLabelFontFile = ''
opt_filepvdDisplay.SelectionPointLabelFontFile = ''
opt_filepvdDisplay.PolarAxes = 'PolarAxesRepresentation'
opt_filepvdDisplay.ScalarOpacityFunction = initTemperaturePWF
opt_filepvdDisplay.ScalarOpacityUnitDistance = 0.05210007309586914

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
opt_filepvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.5003207852936985, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
opt_filepvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.5003207852936985, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
opt_filepvdDisplay.DataAxesGrid.XTitleFontFile = ''
opt_filepvdDisplay.DataAxesGrid.YTitleFontFile = ''
opt_filepvdDisplay.DataAxesGrid.ZTitleFontFile = ''
opt_filepvdDisplay.DataAxesGrid.XLabelFontFile = ''
opt_filepvdDisplay.DataAxesGrid.YLabelFontFile = ''
opt_filepvdDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
opt_filepvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
opt_filepvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
opt_filepvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
opt_filepvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
opt_filepvdDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on initTemperatureLUT
initTemperatureLUT.NumberOfTableValues = 5

# Rescale transfer function
initTemperatureLUT.RescaleTransferFunction(0.0, 0.5)

# Rescale transfer function
initTemperaturePWF.RescaleTransferFunction(0.0, 0.5)

# get color legend/bar for initTemperatureLUT in view renderView1
initTemperatureLUTColorBar = GetScalarBar(initTemperatureLUT, renderView1)
initTemperatureLUTColorBar.Title = 'InitTemperature'
initTemperatureLUTColorBar.ComponentTitle = ''
initTemperatureLUTColorBar.TitleFontFile = ''
initTemperatureLUTColorBar.LabelFontFile = ''

# change scalar bar placement
initTemperatureLUTColorBar.Orientation = 'Horizontal'
initTemperatureLUTColorBar.WindowLocation = 'AnyLocation'
initTemperatureLUTColorBar.Position = [0.20, 0.50]
initTemperatureLUTColorBar.ScalarBarLength = 0.60

# Properties modified on initTemperatureLUTColorBar
initTemperatureLUTColorBar.AutoOrient = 0
initTemperatureLUTColorBar.TitleFontFamily = 'Times'
initTemperatureLUTColorBar.LabelFontFamily = 'Times'

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0


# show color bar/color legend
opt_filepvdDisplay.SetScalarBarVisibility(renderView1, False)


# Rescale transfer function
opt_filepvdDisplay.ScaleTransferFunction.RescaleTransferFunction(0.0937071353487, 0.761981275189)

# Rescale transfer function
opt_filepvdDisplay.OpacityTransferFunction.RescaleTransferFunction(0.0937071353487, 0.761981275189)


# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 2.7320508075688776]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.501

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# update the view to ensure updated data information
renderView1.Update()


#for i in range(int(animationScene1.EndTime)):
#    # save screenshot
#    png_name = my_pvd_file.replace('.pvd', '_{}.png'.format(i))
#    SaveScreenshot(png_name, renderView1, ImageResolution=[2000, 2000], TransparentBackground=1)
#    print('png file generate:', png_name)
#    animationScene1.GoToNext()
for i in [0, int(animationScene1.EndTime)]:
    # save screenshot
    png_name = my_pvd_file.replace('.pvd', '_{}.png'.format(i))
    SaveScreenshot(png_name, renderView1, ImageResolution=[2000, 2000], TransparentBackground=1)
    print('png file generate:', png_name)
    animationScene1.GoToLast()

#
#Hide(opt_filepvd)
## show color bar/color legend
#opt_filepvdDisplay.SetScalarBarVisibility(renderView1, True)
## update the view to ensure updated data information
#renderView1.Update()
#png_name = my_pvd_file.replace('.pvd', '_colorbar.png')
#SaveScreenshot(png_name, renderView1, ImageResolution=[2000, 2000], TransparentBackground=1)
#
