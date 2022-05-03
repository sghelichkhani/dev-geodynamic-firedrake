# trace generated using paraview version 5.9.0-RC1
import os

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
output_path = './PNGs/00_FWD_refmodel'

# create a new 'PVD Reader'
temperaturepvd = PVDReader(registrationName='temperature.pvd', FileName='/Users/sghelichkhani/Data/kerguelen/Workplace/dev-geodynamic-firedrake/adjoint_tic/rectangle/var_mu/outputs/00_FWD_refmodel/temperature.pvd')
temperaturepvd.CellArrays = []
temperaturepvd.PointArrays = ['Temperature', 'Viscosity']
temperaturepvd.ColumnArrays = []


# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# show data in view
temperaturepvdDisplay = Show(temperaturepvd, renderView1, 'UnstructuredGridRepresentation')


#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# set scalar coloring
ColorBy(temperaturepvdDisplay, ('POINTS', 'Viscosity'))

# rescale color and/or opacity maps used to include current data range
temperaturepvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
temperaturepvdDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Viscosity'
viscosityLUT = GetColorTransferFunction('Viscosity')
viscosityLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
viscosityLUT.InterpretValuesAsCategories = 0
viscosityLUT.AnnotationsInitialized = 0
viscosityLUT.ShowCategoricalColorsinDataRangeOnly = 0
viscosityLUT.RescaleOnVisibilityChange = 0
viscosityLUT.EnableOpacityMapping = 0
viscosityLUT.RGBPoints = [0.31622661914007094, 0.231373, 0.298039, 0.752941, 0.6507762498950012, 0.865003, 0.865003, 0.865003, 0.9853258806499312, 0.705882, 0.0156863, 0.14902]
viscosityLUT.UseLogScale = 0
viscosityLUT.UseOpacityControlPointsFreehandDrawing = 0
viscosityLUT.ShowDataHistogram = 0
viscosityLUT.AutomaticDataHistogramComputation = 0
viscosityLUT.DataHistogramNumberOfBins = 10
viscosityLUT.ColorSpace = 'Diverging'
viscosityLUT.UseBelowRangeColor = 0
viscosityLUT.BelowRangeColor = [0.0, 0.0, 0.0]
viscosityLUT.UseAboveRangeColor = 0
viscosityLUT.AboveRangeColor = [0.5, 0.5, 0.5]
viscosityLUT.NanColor = [1.0, 1.0, 0.0]
viscosityLUT.NanOpacity = 1.0
viscosityLUT.Discretize = 1
viscosityLUT.NumberOfTableValues = 256
viscosityLUT.ScalarRangeInitialized = 1.0
viscosityLUT.HSVWrap = 0
viscosityLUT.VectorComponent = 0
viscosityLUT.VectorMode = 'Magnitude'
viscosityLUT.AllowDuplicateScalars = 1
viscosityLUT.Annotations = []
viscosityLUT.ActiveAnnotatedValues = []
viscosityLUT.IndexedColors = []
viscosityLUT.IndexedOpacities = []

# get opacity transfer function/opacity map for 'Viscosity'
viscosityPWF = GetOpacityTransferFunction('Viscosity')
viscosityPWF.Points = [0.31622661914007094, 0.0, 0.5, 0.0, 0.9853258806499312, 1.0, 0.5, 0.0]
viscosityPWF.AllowDuplicateScalars = 1
viscosityPWF.UseLogScale = 0
viscosityPWF.ScalarRangeInitialized = 1

# invert the transfer function
viscosityLUT.InvertTransferFunction()

# convert to log space
viscosityLUT.MapControlPointsToLogSpace()

# Properties modified on viscosityLUT
viscosityLUT.UseLogScale = 1

# invert the transfer function
viscosityLUT.InvertTransferFunction()

# invert the transfer function
viscosityLUT.InvertTransferFunction()

# Properties modified on viscosityLUT
viscosityLUT.NumberOfTableValues = 10

animationScene1.Play()

animationScene1.GoToFirst()

# get color legend/bar for viscosityLUT in view renderView1
viscosityLUTColorBar = GetScalarBar(viscosityLUT, renderView1)
viscosityLUTColorBar.AutoOrient = 1
viscosityLUTColorBar.Orientation = 'Vertical'
viscosityLUTColorBar.WindowLocation = 'LowerRightCorner'
viscosityLUTColorBar.Position = [0.89, 0.02]
viscosityLUTColorBar.Title = 'Viscosity'
viscosityLUTColorBar.ComponentTitle = ''
viscosityLUTColorBar.TitleJustification = 'Centered'
viscosityLUTColorBar.HorizontalTitle = 0
viscosityLUTColorBar.TitleOpacity = 1.0
viscosityLUTColorBar.TitleFontFamily = 'Arial'
viscosityLUTColorBar.TitleFontFile = ''
viscosityLUTColorBar.TitleBold = 0
viscosityLUTColorBar.TitleItalic = 0
viscosityLUTColorBar.TitleShadow = 0
viscosityLUTColorBar.TitleFontSize = 16
viscosityLUTColorBar.LabelOpacity = 1.0
viscosityLUTColorBar.LabelFontFamily = 'Arial'
viscosityLUTColorBar.LabelFontFile = ''
viscosityLUTColorBar.LabelBold = 0
viscosityLUTColorBar.LabelItalic = 0
viscosityLUTColorBar.LabelShadow = 0
viscosityLUTColorBar.LabelFontSize = 16
viscosityLUTColorBar.AutomaticLabelFormat = 1
viscosityLUTColorBar.LabelFormat = '%-#6.3g'
viscosityLUTColorBar.DrawTickMarks = 1
viscosityLUTColorBar.DrawTickLabels = 1
viscosityLUTColorBar.UseCustomLabels = 0
viscosityLUTColorBar.CustomLabels = []
viscosityLUTColorBar.AddRangeLabels = 1
viscosityLUTColorBar.RangeLabelFormat = '%-#6.1e'
viscosityLUTColorBar.DrawAnnotations = 1
viscosityLUTColorBar.AddRangeAnnotations = 0
viscosityLUTColorBar.AutomaticAnnotations = 0
viscosityLUTColorBar.DrawNanAnnotation = 0
viscosityLUTColorBar.NanAnnotation = 'NaN'
viscosityLUTColorBar.TextPosition = 'Ticks right/top, annotations left/bottom'
viscosityLUTColorBar.ReverseLegend = 0
viscosityLUTColorBar.ScalarBarThickness = 16
viscosityLUTColorBar.ScalarBarLength = 0.33

# change scalar bar placement
viscosityLUTColorBar.Orientation = 'Horizontal'
viscosityLUTColorBar.WindowLocation = 'AnyLocation'
viscosityLUTColorBar.Position = [0.3247619047619048, 0.872504472271914]
viscosityLUTColorBar.ScalarBarLength = 0.3299999999999996

# Properties modified on viscosityLUTColorBar
viscosityLUTColorBar.AutoOrient = 0
viscosityLUTColorBar.TitleFontFamily = 'Times'

# Properties modified on viscosityLUTColorBar
viscosityLUTColorBar.ScalarBarThickness = 15

# Properties modified on viscosityLUTColorBar
viscosityLUTColorBar.ScalarBarLength = 0.6

# Rescale transfer function
viscosityLUT.RescaleTransferFunction(0.3, 1.0)

# Rescale transfer function
viscosityPWF.RescaleTransferFunction(0.3, 1.0)

# Properties modified on viscosityLUT
viscosityLUT.NumberOfTableValues = 7

# Properties modified on viscosityLUTColorBar
viscosityLUTColorBar.TitleFontSize = 8 
viscosityLUTColorBar.LabelFontSize = 8

# Properties modified on viscosityLUTColorBar
viscosityLUTColorBar.LabelFontFamily = 'Times'

# Properties modified on viscosityLUTColorBar
viscosityLUTColorBar.Position = [0.2, 0.9]

# get layout
layout1 = GetLayout()

LoadPalette(paletteName='WhiteBackground')


# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5,  8000.0]
renderView1.CameraFocalPoint = [0.5, 0.55, 0.0]
renderView1.CameraParallelScale = 0.6


# update the view to ensure updated data information
renderView1.Update()

os.makedirs(output_path)

for i in range(int(animationScene1.EndTime+1)):
    # save screenshot
    print("Saving figure %5.5i out of %5.5i" %(i, animationScene1.EndTime))
    SaveScreenshot(os.path.join(output_path, str('fig_%2.2i.png' %i)), renderView1, ImageResolution=[1800, 2000],
        FontScaling='Scale fonts proportionally',
        OverrideColorPalette='',
        StereoMode='No change',
        TransparentBackground=0, 
        # PNG options
        CompressionLevel='5')
    animationScene1.GoToNext()

## layout/tab size in pixels
#layout1.SetSize(1806, 1118)
#
## current camera placement for renderView1
#renderView1.InteractionMode = '2D'
#renderView1.CameraPosition = [0.5, 0.5, 10000.0]
#renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
#renderView1.CameraParallelScale = 0.7071067811865476
#
## save screenshot
#SaveScreenshot('/Users/sghelichkhani/Desktop/test.png', renderView1, ImageResolution=[2000, 2000],
#    FontScaling='Scale fonts proportionally',
#    OverrideColorPalette='',
#    StereoMode='No change',
#    TransparentBackground=0, 
#    # PNG options
#    CompressionLevel='5')
