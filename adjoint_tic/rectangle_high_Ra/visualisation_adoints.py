# trace generated using paraview version 5.9.0-RC1
import os

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


case = 'BFGS_linesearh_1'

output_path = './PNGs/'+case+'/'
# create a new 'PVD Reader'
opt_filepvd = PVDReader(registrationName='opt_file.pvd', FileName=\
                        '/home/sia/Workplace/dev-geodynamic-firedrake/adjoint_tic/rectangle_high_Ra/cases/'+case+'/visual/opt_file.pvd')

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0


# show data in view
opt_filepvdDisplay = Show(opt_filepvd, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
opt_filepvdDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()



# set scalar coloring
ColorBy(opt_filepvdDisplay, ('POINTS', 'InitTemperature'))


# rescale color and/or opacity maps used to include current data range
opt_filepvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
opt_filepvdDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'InitTemperature'
initTemperatureLUT = GetColorTransferFunction('InitTemperature')
initTemperatureLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.25016039264684925, 0.865003, 0.865003, 0.865003, 0.5003207852936985, 0.705882, 0.0156863, 0.14902]
initTemperatureLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'InitTemperature'
initTemperaturePWF = GetOpacityTransferFunction('InitTemperature')
initTemperaturePWF.Points = [0.0, 0.0, 0.5, 0.0, 0.5003207852936985, 1.0, 0.5, 0.0]
initTemperaturePWF.ScalarRangeInitialized = 1

# get color transfer function/color map for 'Viscosity'
initTemperatureLUT = GetColorTransferFunction('InitTemperature')
initTemperatureLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
initTemperatureLUT.InterpretValuesAsCategories = 0
initTemperatureLUT.AnnotationsInitialized = 0
initTemperatureLUT.ShowCategoricalColorsinDataRangeOnly = 0
initTemperatureLUT.RescaleOnVisibilityChange = 0
initTemperatureLUT.EnableOpacityMapping = 0
initTemperatureLUT.RGBPoints = [0.31622661914007094, 0.231373, 0.298039, 0.752941, 0.6507762498950012, 0.865003, 0.865003, 0.865003, 0.9853258806499312, 0.705882, 0.0156863, 0.14902]
initTemperatureLUT.UseLogScale = 0
initTemperatureLUT.UseOpacityControlPointsFreehandDrawing = 0
initTemperatureLUT.ShowDataHistogram = 0
initTemperatureLUT.AutomaticDataHistogramComputation = 0
initTemperatureLUT.DataHistogramNumberOfBins = 10
initTemperatureLUT.ColorSpace = 'Diverging'
initTemperatureLUT.UseBelowRangeColor = 0
initTemperatureLUT.BelowRangeColor = [0.0, 0.0, 0.0]
initTemperatureLUT.UseAboveRangeColor = 0
initTemperatureLUT.AboveRangeColor = [0.5, 0.5, 0.5]
initTemperatureLUT.NanColor = [1.0, 1.0, 0.0]
initTemperatureLUT.NanOpacity = 1.0
initTemperatureLUT.Discretize = 1
initTemperatureLUT.NumberOfTableValues = 256
initTemperatureLUT.ScalarRangeInitialized = 1.0
initTemperatureLUT.HSVWrap = 0
initTemperatureLUT.VectorComponent = 0
initTemperatureLUT.VectorMode = 'Magnitude'
initTemperatureLUT.AllowDuplicateScalars = 1
initTemperatureLUT.Annotations = []
initTemperatureLUT.ActiveAnnotatedValues = []
initTemperatureLUT.IndexedColors = []
initTemperatureLUT.IndexedOpacities = []

# get opacity transfer function/opacity map for 'Viscosity'
initTemperaturePWF = GetOpacityTransferFunction('Temperature')
initTemperaturePWF.Points = [0.31622661914007094, 0.0, 0.5, 0.0, 0.9853258806499312, 1.0, 0.5, 0.0]
initTemperaturePWF.AllowDuplicateScalars = 1
initTemperaturePWF.UseLogScale = 0
initTemperaturePWF.ScalarRangeInitialized = 1

# convert to log space
initTemperatureLUT.MapControlPointsToLogSpace()

# Properties modified on initTemperatureLUT
initTemperatureLUT.UseLogScale = 0

# invert the transfer function
initTemperatureLUT.InvertTransferFunction()

# invert the transfer function
initTemperatureLUT.InvertTransferFunction()

# Properties modified on initTemperatureLUT
initTemperatureLUT.NumberOfTableValues = 10

animationScene1.Play()

animationScene1.GoToFirst()

# get color legend/bar for initTemperatureLUT in view renderView1
initTemperatureLUTColorBar = GetScalarBar(initTemperatureLUT, renderView1)
initTemperatureLUTColorBar.AutoOrient = 1
initTemperatureLUTColorBar.Orientation = 'Vertical'
initTemperatureLUTColorBar.WindowLocation = 'LowerRightCorner'
initTemperatureLUTColorBar.Position = [0.89, 0.02]
initTemperatureLUTColorBar.Title = 'Temperature'
initTemperatureLUTColorBar.ComponentTitle = ''
initTemperatureLUTColorBar.TitleJustification = 'Centered'
initTemperatureLUTColorBar.HorizontalTitle = 0
initTemperatureLUTColorBar.TitleOpacity = 1.0
initTemperatureLUTColorBar.TitleFontFamily = 'Arial'
initTemperatureLUTColorBar.TitleFontFile = ''
initTemperatureLUTColorBar.TitleBold = 0
initTemperatureLUTColorBar.TitleItalic = 0
initTemperatureLUTColorBar.TitleShadow = 0
initTemperatureLUTColorBar.TitleFontSize = 16
initTemperatureLUTColorBar.LabelOpacity = 1.0
initTemperatureLUTColorBar.LabelFontFamily = 'Arial'
initTemperatureLUTColorBar.LabelFontFile = ''
initTemperatureLUTColorBar.LabelBold = 0
initTemperatureLUTColorBar.LabelItalic = 0
initTemperatureLUTColorBar.LabelShadow = 0
initTemperatureLUTColorBar.LabelFontSize = 16
initTemperatureLUTColorBar.AutomaticLabelFormat = 1
initTemperatureLUTColorBar.LabelFormat = '%-#6.3g'
initTemperatureLUTColorBar.DrawTickMarks = 1
initTemperatureLUTColorBar.DrawTickLabels = 1
initTemperatureLUTColorBar.UseCustomLabels = 0
initTemperatureLUTColorBar.CustomLabels = []
initTemperatureLUTColorBar.AddRangeLabels = 1
initTemperatureLUTColorBar.RangeLabelFormat = '%-#6.1e'
initTemperatureLUTColorBar.DrawAnnotations = 1
initTemperatureLUTColorBar.AddRangeAnnotations = 0
initTemperatureLUTColorBar.AutomaticAnnotations = 0
initTemperatureLUTColorBar.DrawNanAnnotation = 0
initTemperatureLUTColorBar.NanAnnotation = 'NaN'
initTemperatureLUTColorBar.TextPosition = 'Ticks right/top, annotations left/bottom'
initTemperatureLUTColorBar.ReverseLegend = 0
initTemperatureLUTColorBar.ScalarBarThickness = 16
initTemperatureLUTColorBar.ScalarBarLength = 0.33

# change scalar bar placement
initTemperatureLUTColorBar.Orientation = 'Horizontal'
initTemperatureLUTColorBar.WindowLocation = 'AnyLocation'
initTemperatureLUTColorBar.Position = [0.3247619047619048, 0.872504472271914]
initTemperatureLUTColorBar.ScalarBarLength = 0.3299999999999996

# Properties modified on initTemperatureLUTColorBar
initTemperatureLUTColorBar.AutoOrient = 0
initTemperatureLUTColorBar.TitleFontFamily = 'Times'

# Properties modified on initTemperatureLUTColorBar
initTemperatureLUTColorBar.ScalarBarThickness = 15

# Properties modified on initTemperatureLUTColorBar
initTemperatureLUTColorBar.ScalarBarLength = 0.6

# Rescale transfer function
initTemperatureLUT.RescaleTransferFunction(0.35, 0.55)

# Rescale transfer function
initTemperaturePWF.RescaleTransferFunction(0.35, 0.55)

# Properties modified on initTemperatureLUT
initTemperatureLUT.NumberOfTableValues = 20

# Properties modified on initTemperatureLUTColorBar
initTemperatureLUTColorBar.TitleFontSize = 8 
initTemperatureLUTColorBar.LabelFontSize = 8

# Properties modified on initTemperatureLUTColorBar
initTemperatureLUTColorBar.LabelFontFamily = 'Times'

# Properties modified on initTemperatureLUTColorBar
initTemperatureLUTColorBar.Position = [0.2, 0.9]

# get layout
layout1 = GetLayout()

LoadPalette(paletteName='WhiteBackground')

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5,  8000.0]
renderView1.CameraFocalPoint = [0.5, 0.55, 0.0]
renderView1.CameraParallelScale = 0.6

if not os.path.isdir(output_path):
    os.makedirs(output_path)

for i in range(int(animationScene1.EndTime+1)):
    print("Saving intial figure %5.5i out of %5.5i" %(i, animationScene1.EndTime))
    SaveScreenshot(os.path.join(output_path, str('fig_init_%2.2i.png' %i)), renderView1, ImageResolution=[1800, 2000],
        FontScaling='Scale fonts proportionally',
        OverrideColorPalette='',
        StereoMode='No change',
        TransparentBackground=0, 
        # PNG options
        CompressionLevel='5')
    animationScene1.GoToNext()


animationScene1.GoToFirst()

# set scalar coloring
ColorBy(opt_filepvdDisplay, ('POINTS', 'FinTemperature'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(initTemperatureLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
opt_filepvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
opt_filepvdDisplay.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'FinTemperature'
finTemperaturePWF = GetOpacityTransferFunction('FinTemperature')
finTemperaturePWF.Points = [-0.06588906619449987, 0.0, 0.5, 0.0, 0.506530234960336, 1.0, 0.5, 0.0]
finTemperaturePWF.ScalarRangeInitialized = 1

# get color transfer function/color map for 'Viscosity'
finTemperatureLUT = GetColorTransferFunction('FinTemperature')
finTemperatureLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
finTemperatureLUT.InterpretValuesAsCategories = 0
finTemperatureLUT.AnnotationsInitialized = 0
finTemperatureLUT.ShowCategoricalColorsinDataRangeOnly = 0
finTemperatureLUT.RescaleOnVisibilityChange = 0
finTemperatureLUT.EnableOpacityMapping = 0
finTemperatureLUT.RGBPoints = [0.31622661914007094, 0.231373, 0.298039, 0.752941, 0.6507762498950012, 0.865003, 0.865003, 0.865003, 0.9853258806499312, 0.705882, 0.0156863, 0.14902]
finTemperatureLUT.UseLogScale = 0
finTemperatureLUT.UseOpacityControlPointsFreehandDrawing = 0
finTemperatureLUT.ShowDataHistogram = 0
finTemperatureLUT.AutomaticDataHistogramComputation = 0
finTemperatureLUT.DataHistogramNumberOfBins = 10
finTemperatureLUT.ColorSpace = 'Diverging'
finTemperatureLUT.UseBelowRangeColor = 0
finTemperatureLUT.BelowRangeColor = [0.0, 0.0, 0.0]
finTemperatureLUT.UseAboveRangeColor = 0
finTemperatureLUT.AboveRangeColor = [0.5, 0.5, 0.5]
finTemperatureLUT.NanColor = [1.0, 1.0, 0.0]
finTemperatureLUT.NanOpacity = 1.0
finTemperatureLUT.Discretize = 1
finTemperatureLUT.NumberOfTableValues = 256
finTemperatureLUT.ScalarRangeInitialized = 1.0
finTemperatureLUT.HSVWrap = 0
finTemperatureLUT.VectorComponent = 0
finTemperatureLUT.VectorMode = 'Magnitude'
finTemperatureLUT.AllowDuplicateScalars = 1
finTemperatureLUT.Annotations = []
finTemperatureLUT.ActiveAnnotatedValues = []
finTemperatureLUT.IndexedColors = []
finTemperatureLUT.IndexedOpacities = []

# get opacity transfer function/opacity map for 'Viscosity'
initTemperaturePWF = GetOpacityTransferFunction('Temperature')
initTemperaturePWF.Points = [0.31622661914007094, 0.0, 0.5, 0.0, 0.9853258806499312, 1.0, 0.5, 0.0]
initTemperaturePWF.AllowDuplicateScalars = 1
initTemperaturePWF.UseLogScale = 0
initTemperaturePWF.ScalarRangeInitialized = 1

# convert to log space
finTemperatureLUT.MapControlPointsToLogSpace()

# Properties modified on finTemperatureLUT
finTemperatureLUT.UseLogScale = 0

# invert the transfer function
finTemperatureLUT.InvertTransferFunction()

# invert the transfer function
finTemperatureLUT.InvertTransferFunction()

# Properties modified on finTemperatureLUT
finTemperatureLUT.NumberOfTableValues = 10



# get color legend/bar for initTemperatureLUT in view renderView1
finTemperatureLUTColorBar = GetScalarBar(finTemperatureLUT, renderView1)
finTemperatureLUTColorBar.AutoOrient = 0
finTemperatureLUTColorBar.Orientation = 'Horizontal'
finTemperatureLUTColorBar.WindowLocation = 'AnyLocation'
finTemperatureLUTColorBar.Position = [0.89, 0.02]
finTemperatureLUTColorBar.Title = 'Temperature'
finTemperatureLUTColorBar.ComponentTitle = ''
finTemperatureLUTColorBar.TitleJustification = 'Centered'
finTemperatureLUTColorBar.HorizontalTitle = 0
finTemperatureLUTColorBar.TitleOpacity = 1.0
finTemperatureLUTColorBar.TitleFontFamily = 'Arial'
finTemperatureLUTColorBar.TitleFontFile = ''
finTemperatureLUTColorBar.TitleBold = 0
finTemperatureLUTColorBar.TitleItalic = 0
finTemperatureLUTColorBar.TitleShadow = 0
finTemperatureLUTColorBar.TitleFontSize = 16
finTemperatureLUTColorBar.LabelOpacity = 1.0
finTemperatureLUTColorBar.LabelFontFamily = 'Arial'
finTemperatureLUTColorBar.LabelFontFile = ''
finTemperatureLUTColorBar.LabelBold = 0
finTemperatureLUTColorBar.LabelItalic = 0
finTemperatureLUTColorBar.LabelShadow = 0
finTemperatureLUTColorBar.LabelFontSize = 16
finTemperatureLUTColorBar.AutomaticLabelFormat = 1
finTemperatureLUTColorBar.LabelFormat = '%-#6.3g'
finTemperatureLUTColorBar.DrawTickMarks = 1
finTemperatureLUTColorBar.DrawTickLabels = 1
finTemperatureLUTColorBar.UseCustomLabels = 0
finTemperatureLUTColorBar.CustomLabels = []
finTemperatureLUTColorBar.AddRangeLabels = 1
finTemperatureLUTColorBar.RangeLabelFormat = '%-#6.1e'
finTemperatureLUTColorBar.DrawAnnotations = 1
finTemperatureLUTColorBar.AddRangeAnnotations = 0
finTemperatureLUTColorBar.AutomaticAnnotations = 0
finTemperatureLUTColorBar.DrawNanAnnotation = 0
finTemperatureLUTColorBar.NanAnnotation = 'NaN'
finTemperatureLUTColorBar.TextPosition = 'Ticks right/top, annotations left/bottom'
finTemperatureLUTColorBar.ReverseLegend = 0
finTemperatureLUTColorBar.ScalarBarThickness = 16
finTemperatureLUTColorBar.ScalarBarLength = 0.33

# change scalar bar placement
finTemperatureLUTColorBar.Orientation = 'Horizontal'
finTemperatureLUTColorBar.WindowLocation = 'AnyLocation'
finTemperatureLUTColorBar.Position = [0.3247619047619048, 0.872504472271914]
finTemperatureLUTColorBar.ScalarBarLength = 0.3299999999999996

# Properties modified on finTemperatureLUTColorBar
finTemperatureLUTColorBar.AutoOrient = 0
finTemperatureLUTColorBar.TitleFontFamily = 'Times'

# Properties modified on finTemperatureLUTColorBar
finTemperatureLUTColorBar.ScalarBarThickness = 15

# Properties modified on finTemperatureLUTColorBar
finTemperatureLUTColorBar.ScalarBarLength = 0.6

# Rescale transfer function
finTemperatureLUT.RescaleTransferFunction(0.35, 0.55)

# Rescale transfer function
finTemperaturePWF.RescaleTransferFunction(0.35, 0.55)

# Properties modified on initTemperatureLUT
finTemperatureLUT.NumberOfTableValues = 20

# Properties modified on finTemperatureLUTColorBar
finTemperatureLUTColorBar.TitleFontSize = 8 
finTemperatureLUTColorBar.LabelFontSize = 8

# Properties modified on finTemperatureLUTColorBar
finTemperatureLUTColorBar.LabelFontFamily = 'Times'

# Properties modified on finTemperatureLUTColorBar
finTemperatureLUTColorBar.Position = [0.2, 0.9]


# update the view to ensure updated data information
renderView1.Update()

for i in range(int(animationScene1.EndTime+1)):
    # save screenshot
    print("Saving final figure %5.5i out of %5.5i" %(i, animationScene1.EndTime))
    SaveScreenshot(os.path.join(output_path, str('fig_fin_%2.2i.png' %i)), renderView1, ImageResolution=[1800, 2000],
        FontScaling='Scale fonts proportionally',
        OverrideColorPalette='',
        StereoMode='No change',
        TransparentBackground=0, 
        # PNG options
        CompressionLevel='5')
    animationScene1.GoToNext()
    
