# trace generated using paraview version 5.5.2

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

import glob

fi_names = glob.glob('/Users/sghelichkhani/Workplace/G-ADOPT_old/hessian_tests/hessian_17_MARCH/visual/rect_1e-5_1_to_0/*pvtu')
fi_names.sort()

# create a new 'XML Partitioned Unstructured Grid Reader'
opt_temperature_initial_ = XMLPartitionedUnstructuredGridReader(FileName=fi_names)
opt_temperature_initial_.CellArrayStatus = []
opt_temperature_initial_.PointArrayStatus = ['Temperature']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2150, 1166]

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [1.0, 0.5, 10000.0]
renderView1.CameraFocalPoint = [1.0, 0.5, 0.0]

# get display properties
opt_temperature_initial_Display = GetDisplayProperties(opt_temperature_initial_, view=renderView1)

# get color transfer function/color map for 'Temperature'
temperatureLUT = GetColorTransferFunction('Temperature')
temperatureLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
temperatureLUT.InterpretValuesAsCategories = 0
temperatureLUT.ShowCategoricalColorsinDataRangeOnly = 0
temperatureLUT.RescaleOnVisibilityChange = 0
temperatureLUT.EnableOpacityMapping = 0
temperatureLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.5, 0.865003, 0.865003, 0.865003, 1.0, 0.705882, 0.0156863, 0.14902]
temperatureLUT.UseLogScale = 0
temperatureLUT.ColorSpace = 'Diverging'
temperatureLUT.UseBelowRangeColor = 0
temperatureLUT.BelowRangeColor = [0.0, 0.0, 0.0]
temperatureLUT.UseAboveRangeColor = 0
temperatureLUT.AboveRangeColor = [0.5, 0.5, 0.5]
temperatureLUT.NanColor = [1.0, 1.0, 0.0]
temperatureLUT.Discretize = 1
temperatureLUT.NumberOfTableValues = 256
temperatureLUT.ScalarRangeInitialized = 1.0
temperatureLUT.HSVWrap = 0
temperatureLUT.VectorComponent = 0
temperatureLUT.VectorMode = 'Magnitude'
temperatureLUT.AllowDuplicateScalars = 1
temperatureLUT.Annotations = []
temperatureLUT.ActiveAnnotatedValues = []
temperatureLUT.IndexedColors = []
temperatureLUT.IndexedOpacities = []

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get color legend/bar for temperatureLUT in view renderView1
temperatureLUTColorBar = GetScalarBar(temperatureLUT, renderView1)
temperatureLUTColorBar.AutoOrient = 1
temperatureLUTColorBar.Orientation = 'Vertical'
temperatureLUTColorBar.WindowLocation = 'LowerRightCorner'
temperatureLUTColorBar.Position = [0.89, 0.02]
temperatureLUTColorBar.Title = 'Temperature'
temperatureLUTColorBar.ComponentTitle = ''
temperatureLUTColorBar.TitleJustification = 'Centered'
temperatureLUTColorBar.HorizontalTitle = 0
temperatureLUTColorBar.TitleColor = [1.0, 1.0, 1.0]
temperatureLUTColorBar.TitleOpacity = 1.0
temperatureLUTColorBar.TitleFontFamily = 'Arial'
temperatureLUTColorBar.TitleFontFile = ''
temperatureLUTColorBar.TitleBold = 0
temperatureLUTColorBar.TitleItalic = 0
temperatureLUTColorBar.TitleShadow = 0
temperatureLUTColorBar.TitleFontSize = 2 
temperatureLUTColorBar.LabelColor = [1.0, 1.0, 1.0]
temperatureLUTColorBar.LabelOpacity = 1.0
temperatureLUTColorBar.LabelFontFamily = 'Arial'
temperatureLUTColorBar.LabelFontFile = ''
temperatureLUTColorBar.LabelBold = 0
temperatureLUTColorBar.LabelItalic = 0
temperatureLUTColorBar.LabelShadow = 0
temperatureLUTColorBar.LabelFontSize = 2
temperatureLUTColorBar.AutomaticLabelFormat = 1
temperatureLUTColorBar.LabelFormat = '%-#6.3g'
temperatureLUTColorBar.DrawTickMarks = 1
temperatureLUTColorBar.DrawTickLabels = 1
temperatureLUTColorBar.UseCustomLabels = 0
temperatureLUTColorBar.CustomLabels = []
temperatureLUTColorBar.AddRangeLabels = 1
temperatureLUTColorBar.RangeLabelFormat = '%-#6.1e'
temperatureLUTColorBar.DrawAnnotations = 1
temperatureLUTColorBar.AddRangeAnnotations = 0
temperatureLUTColorBar.AutomaticAnnotations = 0
temperatureLUTColorBar.DrawNanAnnotation = 0
temperatureLUTColorBar.NanAnnotation = 'NaN'
temperatureLUTColorBar.TextPosition = 'Ticks right/top, annotations left/bottom'
temperatureLUTColorBar.ScalarBarThickness = 16
temperatureLUTColorBar.ScalarBarLength = 0.33

# change scalar bar placement
temperatureLUTColorBar.Orientation = 'Horizontal'
temperatureLUTColorBar.WindowLocation = 'AnyLocation'
temperatureLUTColorBar.Position = [0.4025475146198829, 0.09433962264150944]
temperatureLUTColorBar.ScalarBarLength = 0.33000000000000007

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
temperatureLUT.ApplyPreset('Cool to Warm (Extended)', True)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [1.0, 0.5, 10000.0]
renderView1.CameraFocalPoint = [1.0, 0.5, 0.0]
renderView1.CameraParallelScale = 1.118033988749895

for fi in fi_names:
    # save screenshot
    SaveScreenshot(fi.replace('.pvtu', '.png'), renderView1, ImageResolution=[1368, 1166],
        FontScaling='Scale fonts proportionally',
        OverrideColorPalette='',
        StereoMode='No change',
        TransparentBackground=0, 
        # PNG options
        CompressionLevel='5')
    animationScene1.GoToNext()


