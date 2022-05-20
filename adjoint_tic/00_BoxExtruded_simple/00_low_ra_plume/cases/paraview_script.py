"""
    Visualising the initial and final reconstruction states
    In addition to reference intial and final states

"""

# trace generated using paraview version 5.9.1
import os
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# 
init_fi_path = "/Users/sghelichkhani/Data/Remotes/kerguelen_data/Workplace/dev-geodynamic-firedrake/adjoint_tic/regularisation/SingleProc/cases/00_timesteps_20/visual_020/"
init_fi_name = "opt_temperature_int.pvd"

final_fi_path = init_fi_path 
final_fi_name = "opt_temperature_fin.pvd"

# create a new 'PVD Reader'
opt_temperature_finpvd = PVDReader(registrationName=final_fi_name, FileName=os.path.join(final_fi_path, final_fi_name))
opt_temperature_finpvd.CellArrays = []
opt_temperature_finpvd.PointArrays = ['Temperature', 'RefTemperature']
opt_temperature_finpvd.ColumnArrays = []

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
opt_temperature_finpvdDisplay = Show(opt_temperature_finpvd, renderView1, 'UnstructuredGridRepresentation')

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# get the material library
materialLibrary1 = GetMaterialLibrary()

# get color transfer function/color map for 'Temperature'
temperatureLUT = GetColorTransferFunction('Temperature')
temperatureLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
temperatureLUT.InterpretValuesAsCategories = 0
temperatureLUT.AnnotationsInitialized = 0
temperatureLUT.ShowCategoricalColorsinDataRangeOnly = 0
temperatureLUT.RescaleOnVisibilityChange = 0
temperatureLUT.EnableOpacityMapping = 0
temperatureLUT.RGBPoints = [0.43253340082384273, 0.231373, 0.298039, 0.752941, 0.46626670041192136, 0.865003, 0.865003, 0.865003, 0.5, 0.705882, 0.0156863, 0.14902]
temperatureLUT.UseLogScale = 0
temperatureLUT.UseOpacityControlPointsFreehandDrawing = 0
temperatureLUT.ShowDataHistogram = 0
temperatureLUT.AutomaticDataHistogramComputation = 0
temperatureLUT.DataHistogramNumberOfBins = 10
temperatureLUT.ColorSpace = 'Diverging'
temperatureLUT.UseBelowRangeColor = 0
temperatureLUT.BelowRangeColor = [0.0, 0.0, 0.0]
temperatureLUT.UseAboveRangeColor = 0
temperatureLUT.AboveRangeColor = [0.5, 0.5, 0.5]
temperatureLUT.NanColor = [1.0, 1.0, 0.0]
temperatureLUT.NanOpacity = 1.0
temperatureLUT.Discretize = 1
temperatureLUT.NumberOfTableValues = 10 
temperatureLUT.ScalarRangeInitialized = 1.0
temperatureLUT.HSVWrap = 0
temperatureLUT.VectorComponent = 0
temperatureLUT.VectorMode = 'Magnitude'
temperatureLUT.AllowDuplicateScalars = 1
temperatureLUT.Annotations = []
temperatureLUT.ActiveAnnotatedValues = []
temperatureLUT.IndexedColors = []
temperatureLUT.IndexedOpacities = []

# get opacity transfer function/opacity map for 'Temperature'
temperaturePWF = GetOpacityTransferFunction('Temperature')
temperaturePWF.Points = [0.43253340082384273, 0.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0]
temperaturePWF.AllowDuplicateScalars = 1
temperaturePWF.UseLogScale = 0
temperaturePWF.ScalarRangeInitialized = 1

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create a new 'PVD Reader'
opt_temperature_intpvd = PVDReader(registrationName=init_fi_name, FileName=os.path.join(init_fi_path, init_fi_name))
opt_temperature_intpvd.CellArrays = []
opt_temperature_intpvd.PointArrays = ['Temperature', 'InitTemperature_Ref']
opt_temperature_intpvd.ColumnArrays = []

# show data in view
opt_temperature_intpvdDisplay = Show(opt_temperature_intpvd, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
opt_temperature_intpvdDisplay.Selection = None
opt_temperature_intpvdDisplay.Representation = 'Surface'
opt_temperature_intpvdDisplay.ColorArrayName = ['POINTS', 'Temperature']
opt_temperature_intpvdDisplay.LookupTable = temperatureLUT
opt_temperature_intpvdDisplay.MapScalars = 1
opt_temperature_intpvdDisplay.MultiComponentsMapping = 0
opt_temperature_intpvdDisplay.InterpolateScalarsBeforeMapping = 1
opt_temperature_intpvdDisplay.Opacity = 1.0
opt_temperature_intpvdDisplay.PointSize = 2.0
opt_temperature_intpvdDisplay.LineWidth = 1.0
opt_temperature_intpvdDisplay.RenderLinesAsTubes = 0
opt_temperature_intpvdDisplay.RenderPointsAsSpheres = 0
opt_temperature_intpvdDisplay.Interpolation = 'Gouraud'
opt_temperature_intpvdDisplay.Specular = 0.0
opt_temperature_intpvdDisplay.SpecularColor = [1.0, 1.0, 1.0]
opt_temperature_intpvdDisplay.SpecularPower = 100.0
opt_temperature_intpvdDisplay.Luminosity = 0.0
opt_temperature_intpvdDisplay.Ambient = 0.0
opt_temperature_intpvdDisplay.Diffuse = 1.0
opt_temperature_intpvdDisplay.Roughness = 0.3
opt_temperature_intpvdDisplay.Metallic = 0.0
opt_temperature_intpvdDisplay.EdgeTint = [1.0, 1.0, 1.0]
opt_temperature_intpvdDisplay.SelectTCoordArray = 'None'
opt_temperature_intpvdDisplay.SelectNormalArray = 'None'
opt_temperature_intpvdDisplay.SelectTangentArray = 'None'
opt_temperature_intpvdDisplay.Texture = None
opt_temperature_intpvdDisplay.RepeatTextures = 1
opt_temperature_intpvdDisplay.InterpolateTextures = 0
opt_temperature_intpvdDisplay.SeamlessU = 0
opt_temperature_intpvdDisplay.SeamlessV = 0
opt_temperature_intpvdDisplay.UseMipmapTextures = 0
opt_temperature_intpvdDisplay.BaseColorTexture = None
opt_temperature_intpvdDisplay.NormalTexture = None
opt_temperature_intpvdDisplay.NormalScale = 1.0
opt_temperature_intpvdDisplay.MaterialTexture = None
opt_temperature_intpvdDisplay.OcclusionStrength = 1.0
opt_temperature_intpvdDisplay.EmissiveTexture = None
opt_temperature_intpvdDisplay.EmissiveFactor = [1.0, 1.0, 1.0]
opt_temperature_intpvdDisplay.FlipTextures = 0
opt_temperature_intpvdDisplay.BackfaceRepresentation = 'Follow Frontface'
opt_temperature_intpvdDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
opt_temperature_intpvdDisplay.BackfaceOpacity = 1.0
opt_temperature_intpvdDisplay.Position = [0.0, 0.0, 0.0]
opt_temperature_intpvdDisplay.Scale = [1.0, 1.0, 1.0]
opt_temperature_intpvdDisplay.Orientation = [0.0, 0.0, 0.0]
opt_temperature_intpvdDisplay.Origin = [0.0, 0.0, 0.0]
opt_temperature_intpvdDisplay.CoordinateShiftScaleMethod = 'Always Auto Shift Scale'
opt_temperature_intpvdDisplay.Pickable = 1
opt_temperature_intpvdDisplay.Triangulate = 0
opt_temperature_intpvdDisplay.UseShaderReplacements = 0
opt_temperature_intpvdDisplay.ShaderReplacements = ''
opt_temperature_intpvdDisplay.NonlinearSubdivisionLevel = 1
opt_temperature_intpvdDisplay.UseDataPartitions = 0
opt_temperature_intpvdDisplay.OSPRayUseScaleArray = 'All Approximate'
opt_temperature_intpvdDisplay.OSPRayScaleArray = 'Temperature'
opt_temperature_intpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
opt_temperature_intpvdDisplay.OSPRayMaterial = 'None'
opt_temperature_intpvdDisplay.Orient = 0
opt_temperature_intpvdDisplay.OrientationMode = 'Direction'
opt_temperature_intpvdDisplay.SelectOrientationVectors = 'None'
opt_temperature_intpvdDisplay.Scaling = 0
opt_temperature_intpvdDisplay.ScaleMode = 'No Data Scaling Off'
opt_temperature_intpvdDisplay.ScaleFactor = 0.1
opt_temperature_intpvdDisplay.SelectScaleArray = 'Temperature'
opt_temperature_intpvdDisplay.GlyphType = 'Arrow'
opt_temperature_intpvdDisplay.UseGlyphTable = 0
opt_temperature_intpvdDisplay.GlyphTableIndexArray = 'Temperature'
opt_temperature_intpvdDisplay.UseCompositeGlyphTable = 0
opt_temperature_intpvdDisplay.UseGlyphCullingAndLOD = 0
opt_temperature_intpvdDisplay.LODValues = []
opt_temperature_intpvdDisplay.ColorByLODIndex = 0
opt_temperature_intpvdDisplay.GaussianRadius = 0.005
opt_temperature_intpvdDisplay.ShaderPreset = 'Sphere'
opt_temperature_intpvdDisplay.CustomTriangleScale = 3
opt_temperature_intpvdDisplay.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
opt_temperature_intpvdDisplay.Emissive = 0
opt_temperature_intpvdDisplay.ScaleByArray = 0
opt_temperature_intpvdDisplay.SetScaleArray = ['POINTS', 'Temperature']
opt_temperature_intpvdDisplay.ScaleArrayComponent = ''
opt_temperature_intpvdDisplay.UseScaleFunction = 1
opt_temperature_intpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
opt_temperature_intpvdDisplay.OpacityByArray = 0
opt_temperature_intpvdDisplay.OpacityArray = ['POINTS', 'Temperature']
opt_temperature_intpvdDisplay.OpacityArrayComponent = ''
opt_temperature_intpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
opt_temperature_intpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
opt_temperature_intpvdDisplay.SelectionCellLabelBold = 0
opt_temperature_intpvdDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
opt_temperature_intpvdDisplay.SelectionCellLabelFontFamily = 'Arial'
opt_temperature_intpvdDisplay.SelectionCellLabelFontFile = ''
opt_temperature_intpvdDisplay.SelectionCellLabelFontSize = 18
opt_temperature_intpvdDisplay.SelectionCellLabelItalic = 0
opt_temperature_intpvdDisplay.SelectionCellLabelJustification = 'Left'
opt_temperature_intpvdDisplay.SelectionCellLabelOpacity = 1.0
opt_temperature_intpvdDisplay.SelectionCellLabelShadow = 0
opt_temperature_intpvdDisplay.SelectionPointLabelBold = 0
opt_temperature_intpvdDisplay.SelectionPointLabelColor = [1.0, 1.0, 0.0]
opt_temperature_intpvdDisplay.SelectionPointLabelFontFamily = 'Arial'
opt_temperature_intpvdDisplay.SelectionPointLabelFontFile = ''
opt_temperature_intpvdDisplay.SelectionPointLabelFontSize = 18
opt_temperature_intpvdDisplay.SelectionPointLabelItalic = 0
opt_temperature_intpvdDisplay.SelectionPointLabelJustification = 'Left'
opt_temperature_intpvdDisplay.SelectionPointLabelOpacity = 1.0
opt_temperature_intpvdDisplay.SelectionPointLabelShadow = 0
opt_temperature_intpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
opt_temperature_intpvdDisplay.ScalarOpacityFunction = temperaturePWF
opt_temperature_intpvdDisplay.ScalarOpacityUnitDistance = 0.06564197879454707
opt_temperature_intpvdDisplay.UseSeparateOpacityArray = 0
opt_temperature_intpvdDisplay.OpacityArrayName = ['POINTS', 'Temperature']
opt_temperature_intpvdDisplay.OpacityComponent = ''
opt_temperature_intpvdDisplay.ExtractedBlockIndex = 0
opt_temperature_intpvdDisplay.SelectMapper = 'Projected tetra'
opt_temperature_intpvdDisplay.SamplingDimensions = [128, 128, 128]
opt_temperature_intpvdDisplay.UseFloatingPointFrameBuffer = 1

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
opt_temperature_intpvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
opt_temperature_intpvdDisplay.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
opt_temperature_intpvdDisplay.GlyphType.TipResolution = 6
opt_temperature_intpvdDisplay.GlyphType.TipRadius = 0.1
opt_temperature_intpvdDisplay.GlyphType.TipLength = 0.35
opt_temperature_intpvdDisplay.GlyphType.ShaftResolution = 6
opt_temperature_intpvdDisplay.GlyphType.ShaftRadius = 0.03
opt_temperature_intpvdDisplay.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
opt_temperature_intpvdDisplay.ScaleTransferFunction.Points = [0.41918302117959044, 0.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0]
opt_temperature_intpvdDisplay.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
opt_temperature_intpvdDisplay.OpacityTransferFunction.Points = [0.41918302117959044, 0.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0]
opt_temperature_intpvdDisplay.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
opt_temperature_intpvdDisplay.DataAxesGrid.XTitle = 'X Axis'
opt_temperature_intpvdDisplay.DataAxesGrid.YTitle = 'Y Axis'
opt_temperature_intpvdDisplay.DataAxesGrid.ZTitle = 'Z Axis'
opt_temperature_intpvdDisplay.DataAxesGrid.XTitleFontFamily = 'Arial'
opt_temperature_intpvdDisplay.DataAxesGrid.XTitleFontFile = ''
opt_temperature_intpvdDisplay.DataAxesGrid.XTitleBold = 0
opt_temperature_intpvdDisplay.DataAxesGrid.XTitleItalic = 0
opt_temperature_intpvdDisplay.DataAxesGrid.XTitleFontSize = 12
opt_temperature_intpvdDisplay.DataAxesGrid.XTitleShadow = 0
opt_temperature_intpvdDisplay.DataAxesGrid.XTitleOpacity = 1.0
opt_temperature_intpvdDisplay.DataAxesGrid.YTitleFontFamily = 'Arial'
opt_temperature_intpvdDisplay.DataAxesGrid.YTitleFontFile = ''
opt_temperature_intpvdDisplay.DataAxesGrid.YTitleBold = 0
opt_temperature_intpvdDisplay.DataAxesGrid.YTitleItalic = 0
opt_temperature_intpvdDisplay.DataAxesGrid.YTitleFontSize = 12
opt_temperature_intpvdDisplay.DataAxesGrid.YTitleShadow = 0
opt_temperature_intpvdDisplay.DataAxesGrid.YTitleOpacity = 1.0
opt_temperature_intpvdDisplay.DataAxesGrid.ZTitleFontFamily = 'Arial'
opt_temperature_intpvdDisplay.DataAxesGrid.ZTitleFontFile = ''
opt_temperature_intpvdDisplay.DataAxesGrid.ZTitleBold = 0
opt_temperature_intpvdDisplay.DataAxesGrid.ZTitleItalic = 0
opt_temperature_intpvdDisplay.DataAxesGrid.ZTitleFontSize = 12
opt_temperature_intpvdDisplay.DataAxesGrid.ZTitleShadow = 0
opt_temperature_intpvdDisplay.DataAxesGrid.ZTitleOpacity = 1.0
opt_temperature_intpvdDisplay.DataAxesGrid.FacesToRender = 63
opt_temperature_intpvdDisplay.DataAxesGrid.CullBackface = 0
opt_temperature_intpvdDisplay.DataAxesGrid.CullFrontface = 1
opt_temperature_intpvdDisplay.DataAxesGrid.ShowGrid = 0
opt_temperature_intpvdDisplay.DataAxesGrid.ShowEdges = 1
opt_temperature_intpvdDisplay.DataAxesGrid.ShowTicks = 1
opt_temperature_intpvdDisplay.DataAxesGrid.LabelUniqueEdgesOnly = 1
opt_temperature_intpvdDisplay.DataAxesGrid.AxesToLabel = 63
opt_temperature_intpvdDisplay.DataAxesGrid.XLabelFontFamily = 'Arial'
opt_temperature_intpvdDisplay.DataAxesGrid.XLabelFontFile = ''
opt_temperature_intpvdDisplay.DataAxesGrid.XLabelBold = 0
opt_temperature_intpvdDisplay.DataAxesGrid.XLabelItalic = 0
opt_temperature_intpvdDisplay.DataAxesGrid.XLabelFontSize = 12
opt_temperature_intpvdDisplay.DataAxesGrid.XLabelShadow = 0
opt_temperature_intpvdDisplay.DataAxesGrid.XLabelOpacity = 1.0
opt_temperature_intpvdDisplay.DataAxesGrid.YLabelFontFamily = 'Arial'
opt_temperature_intpvdDisplay.DataAxesGrid.YLabelFontFile = ''
opt_temperature_intpvdDisplay.DataAxesGrid.YLabelBold = 0
opt_temperature_intpvdDisplay.DataAxesGrid.YLabelItalic = 0
opt_temperature_intpvdDisplay.DataAxesGrid.YLabelFontSize = 12
opt_temperature_intpvdDisplay.DataAxesGrid.YLabelShadow = 0
opt_temperature_intpvdDisplay.DataAxesGrid.YLabelOpacity = 1.0
opt_temperature_intpvdDisplay.DataAxesGrid.ZLabelFontFamily = 'Arial'
opt_temperature_intpvdDisplay.DataAxesGrid.ZLabelFontFile = ''
opt_temperature_intpvdDisplay.DataAxesGrid.ZLabelBold = 0
opt_temperature_intpvdDisplay.DataAxesGrid.ZLabelItalic = 0
opt_temperature_intpvdDisplay.DataAxesGrid.ZLabelFontSize = 12
opt_temperature_intpvdDisplay.DataAxesGrid.ZLabelShadow = 0
opt_temperature_intpvdDisplay.DataAxesGrid.ZLabelOpacity = 1.0
opt_temperature_intpvdDisplay.DataAxesGrid.XAxisNotation = 'Mixed'
opt_temperature_intpvdDisplay.DataAxesGrid.XAxisPrecision = 2
opt_temperature_intpvdDisplay.DataAxesGrid.XAxisUseCustomLabels = 0
opt_temperature_intpvdDisplay.DataAxesGrid.XAxisLabels = []
opt_temperature_intpvdDisplay.DataAxesGrid.YAxisNotation = 'Mixed'
opt_temperature_intpvdDisplay.DataAxesGrid.YAxisPrecision = 2
opt_temperature_intpvdDisplay.DataAxesGrid.YAxisUseCustomLabels = 0
opt_temperature_intpvdDisplay.DataAxesGrid.YAxisLabels = []
opt_temperature_intpvdDisplay.DataAxesGrid.ZAxisNotation = 'Mixed'
opt_temperature_intpvdDisplay.DataAxesGrid.ZAxisPrecision = 2
opt_temperature_intpvdDisplay.DataAxesGrid.ZAxisUseCustomLabels = 0
opt_temperature_intpvdDisplay.DataAxesGrid.ZAxisLabels = []
opt_temperature_intpvdDisplay.DataAxesGrid.UseCustomBounds = 0
opt_temperature_intpvdDisplay.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
opt_temperature_intpvdDisplay.PolarAxes.Visibility = 0
opt_temperature_intpvdDisplay.PolarAxes.Translation = [0.0, 0.0, 0.0]
opt_temperature_intpvdDisplay.PolarAxes.Scale = [1.0, 1.0, 1.0]
opt_temperature_intpvdDisplay.PolarAxes.Orientation = [0.0, 0.0, 0.0]
opt_temperature_intpvdDisplay.PolarAxes.EnableCustomBounds = [0, 0, 0]
opt_temperature_intpvdDisplay.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
opt_temperature_intpvdDisplay.PolarAxes.EnableCustomRange = 0
opt_temperature_intpvdDisplay.PolarAxes.CustomRange = [0.0, 1.0]
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisVisibility = 1
opt_temperature_intpvdDisplay.PolarAxes.RadialAxesVisibility = 1
opt_temperature_intpvdDisplay.PolarAxes.DrawRadialGridlines = 1
opt_temperature_intpvdDisplay.PolarAxes.PolarArcsVisibility = 1
opt_temperature_intpvdDisplay.PolarAxes.DrawPolarArcsGridlines = 1
opt_temperature_intpvdDisplay.PolarAxes.NumberOfRadialAxes = 0
opt_temperature_intpvdDisplay.PolarAxes.AutoSubdividePolarAxis = 1
opt_temperature_intpvdDisplay.PolarAxes.NumberOfPolarAxis = 0
opt_temperature_intpvdDisplay.PolarAxes.MinimumRadius = 0.0
opt_temperature_intpvdDisplay.PolarAxes.MinimumAngle = 0.0
opt_temperature_intpvdDisplay.PolarAxes.MaximumAngle = 90.0
opt_temperature_intpvdDisplay.PolarAxes.RadialAxesOriginToPolarAxis = 1
opt_temperature_intpvdDisplay.PolarAxes.Ratio = 1.0
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
opt_temperature_intpvdDisplay.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
opt_temperature_intpvdDisplay.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
opt_temperature_intpvdDisplay.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTitleVisibility = 1
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTitle = 'Radial Distance'
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTitleLocation = 'Bottom'
opt_temperature_intpvdDisplay.PolarAxes.PolarLabelVisibility = 1
opt_temperature_intpvdDisplay.PolarAxes.PolarLabelFormat = '%-#6.3g'
opt_temperature_intpvdDisplay.PolarAxes.PolarLabelExponentLocation = 'Labels'
opt_temperature_intpvdDisplay.PolarAxes.RadialLabelVisibility = 1
opt_temperature_intpvdDisplay.PolarAxes.RadialLabelFormat = '%-#3.1f'
opt_temperature_intpvdDisplay.PolarAxes.RadialLabelLocation = 'Bottom'
opt_temperature_intpvdDisplay.PolarAxes.RadialUnitsVisibility = 1
opt_temperature_intpvdDisplay.PolarAxes.ScreenSize = 10.0
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTitleOpacity = 1.0
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTitleBold = 0
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTitleItalic = 0
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTitleShadow = 0
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTitleFontSize = 12
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisLabelOpacity = 1.0
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisLabelBold = 0
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisLabelItalic = 0
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisLabelShadow = 0
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisLabelFontSize = 12
opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisTextOpacity = 1.0
opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisTextBold = 0
opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisTextItalic = 0
opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisTextShadow = 0
opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisTextFontSize = 12
opt_temperature_intpvdDisplay.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
opt_temperature_intpvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
opt_temperature_intpvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''
opt_temperature_intpvdDisplay.PolarAxes.SecondaryRadialAxesTextBold = 0
opt_temperature_intpvdDisplay.PolarAxes.SecondaryRadialAxesTextItalic = 0
opt_temperature_intpvdDisplay.PolarAxes.SecondaryRadialAxesTextShadow = 0
opt_temperature_intpvdDisplay.PolarAxes.SecondaryRadialAxesTextFontSize = 12
opt_temperature_intpvdDisplay.PolarAxes.EnableDistanceLOD = 1
opt_temperature_intpvdDisplay.PolarAxes.DistanceLODThreshold = 0.7
opt_temperature_intpvdDisplay.PolarAxes.EnableViewAngleLOD = 1
opt_temperature_intpvdDisplay.PolarAxes.ViewAngleLODThreshold = 0.7
opt_temperature_intpvdDisplay.PolarAxes.SmallestVisiblePolarAngle = 0.5
opt_temperature_intpvdDisplay.PolarAxes.PolarTicksVisibility = 1
opt_temperature_intpvdDisplay.PolarAxes.ArcTicksOriginToPolarAxis = 1
opt_temperature_intpvdDisplay.PolarAxes.TickLocation = 'Both'
opt_temperature_intpvdDisplay.PolarAxes.AxisTickVisibility = 1
opt_temperature_intpvdDisplay.PolarAxes.AxisMinorTickVisibility = 0
opt_temperature_intpvdDisplay.PolarAxes.ArcTickVisibility = 1
opt_temperature_intpvdDisplay.PolarAxes.ArcMinorTickVisibility = 0
opt_temperature_intpvdDisplay.PolarAxes.DeltaAngleMajor = 10.0
opt_temperature_intpvdDisplay.PolarAxes.DeltaAngleMinor = 5.0
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisMajorTickSize = 0.0
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTickRatioSize = 0.3
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisMajorTickThickness = 1.0
opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTickRatioThickness = 0.5
opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisMajorTickSize = 0.0
opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisTickRatioSize = 0.3
opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
opt_temperature_intpvdDisplay.PolarAxes.ArcMajorTickSize = 0.0
opt_temperature_intpvdDisplay.PolarAxes.ArcTickRatioSize = 0.3
opt_temperature_intpvdDisplay.PolarAxes.ArcMajorTickThickness = 1.0
opt_temperature_intpvdDisplay.PolarAxes.ArcTickRatioThickness = 0.5
opt_temperature_intpvdDisplay.PolarAxes.Use2DMode = 0
opt_temperature_intpvdDisplay.PolarAxes.UseLogAxis = 0

# show color bar/color legend
opt_temperature_intpvdDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
Hide(opt_temperature_intpvd, renderView1)

# Rescale transfer function
temperatureLUT.RescaleTransferFunction(0.4, 0.5)

# Rescale transfer function
temperaturePWF.RescaleTransferFunction(0.4, 0.5)

# get color legend/bar for temperatureLUT in view renderView1
temperatureLUTColorBar = GetScalarBar(temperatureLUT, renderView1)
temperatureLUTColorBar.AutoOrient = 0
temperatureLUTColorBar.Orientation = 'Horizontal'
temperatureLUTColorBar.WindowLocation = 'AnyLocation'
temperatureLUTColorBar.Position = [0.3, 0.9]
temperatureLUTColorBar.Title = 'Temperature'
temperatureLUTColorBar.ComponentTitle = ''
temperatureLUTColorBar.TitleJustification = 'Centered'
temperatureLUTColorBar.HorizontalTitle = 0
temperatureLUTColorBar.TitleOpacity = 1.0
temperatureLUTColorBar.TitleFontFamily = 'Times'
temperatureLUTColorBar.TitleFontFile = ''
temperatureLUTColorBar.TitleBold = 0
temperatureLUTColorBar.TitleItalic = 0
temperatureLUTColorBar.TitleShadow = 0
temperatureLUTColorBar.TitleFontSize = 16
temperatureLUTColorBar.LabelOpacity = 1.0
temperatureLUTColorBar.LabelFontFamily = 'Times'
temperatureLUTColorBar.LabelFontFile = ''
temperatureLUTColorBar.LabelBold = 0
temperatureLUTColorBar.LabelItalic = 0
temperatureLUTColorBar.LabelShadow = 0
temperatureLUTColorBar.LabelFontSize = 16
temperatureLUTColorBar.AutomaticLabelFormat = 0
temperatureLUTColorBar.LabelFormat = '%4.2f'
temperatureLUTColorBar.DrawTickMarks = 1
temperatureLUTColorBar.DrawTickLabels = 1
temperatureLUTColorBar.UseCustomLabels = 0
temperatureLUTColorBar.CustomLabels = []
temperatureLUTColorBar.AddRangeLabels = 1
temperatureLUTColorBar.RangeLabelFormat = '%4.1f'
temperatureLUTColorBar.DrawAnnotations = 1
temperatureLUTColorBar.AddRangeAnnotations = 0
temperatureLUTColorBar.AutomaticAnnotations = 0
temperatureLUTColorBar.DrawNanAnnotation = 0
temperatureLUTColorBar.NanAnnotation = 'NaN'
temperatureLUTColorBar.TextPosition = 'Ticks right/top, annotations left/bottom'
temperatureLUTColorBar.ReverseLegend = 0
temperatureLUTColorBar.ScalarBarThickness = 20
temperatureLUTColorBar.ScalarBarLength = 0.4


# set active source
SetActiveSource(opt_temperature_finpvd)

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(1673, 1164)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476

# updating the renderview
renderView1.Update()

# save screenshot
for i in range(int(animationScene1.EndTime)+1):
    SaveScreenshot(os.path.join(final_fi_path, str('./Rec_Fin_%3.3i.png' %(i))), renderView1, ImageResolution=[1673, 1164],
        FontScaling='Scale fonts proportionally',
        OverrideColorPalette='',
        StereoMode='No change',
        TransparentBackground=1, 
        # PNG options
        CompressionLevel='5')
    animationScene1.GoToNext()


# set active source
SetActiveSource(opt_temperature_finpvd)

# get display properties
opt_temperature_finpvdDisplay = GetDisplayProperties(opt_temperature_finpvd, view=renderView1)

# set scalar coloring
ColorBy(opt_temperature_finpvdDisplay, ('POINTS', 'RefTemperature'))

#Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(temperatureLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
opt_temperature_intpvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
opt_temperature_intpvdDisplay.SetScalarBarVisibility(renderView1, False)

# get color transfer function/color map for 'Temperature'
Temperature_RefLUT = GetColorTransferFunction('RefTemperature')
Temperature_RefLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
Temperature_RefLUT.InterpretValuesAsCategories = 0
Temperature_RefLUT.AnnotationsInitialized = 0
Temperature_RefLUT.ShowCategoricalColorsinDataRangeOnly = 0
Temperature_RefLUT.RescaleOnVisibilityChange = 0
Temperature_RefLUT.EnableOpacityMapping = 0
Temperature_RefLUT.RGBPoints = [0.43253340082384273, 0.231373, 0.298039, 0.752941, 0.46626670041192136, 0.865003, 0.865003, 0.865003, 0.5, 0.705882, 0.0156863, 0.14902]
Temperature_RefLUT.UseLogScale = 0
Temperature_RefLUT.UseOpacityControlPointsFreehandDrawing = 0
Temperature_RefLUT.ShowDataHistogram = 0
Temperature_RefLUT.AutomaticDataHistogramComputation = 0
Temperature_RefLUT.DataHistogramNumberOfBins = 10
Temperature_RefLUT.ColorSpace = 'Diverging'
Temperature_RefLUT.UseBelowRangeColor = 0
Temperature_RefLUT.BelowRangeColor = [0.0, 0.0, 0.0]
Temperature_RefLUT.UseAboveRangeColor = 0
Temperature_RefLUT.AboveRangeColor = [0.5, 0.5, 0.5]
Temperature_RefLUT.NanColor = [1.0, 1.0, 0.0]
Temperature_RefLUT.NanOpacity = 1.0
Temperature_RefLUT.Discretize = 1
Temperature_RefLUT.NumberOfTableValues = 10 
Temperature_RefLUT.ScalarRangeInitialized = 1.0
Temperature_RefLUT.HSVWrap = 0
Temperature_RefLUT.VectorComponent = 0
Temperature_RefLUT.VectorMode = 'Magnitude'
Temperature_RefLUT.AllowDuplicateScalars = 1
Temperature_RefLUT.Annotations = []
Temperature_RefLUT.ActiveAnnotatedValues = []
Temperature_RefLUT.IndexedColors = []
Temperature_RefLUT.IndexedOpacities = []

# get opacity transfer function/opacity map for 'Temperature'
Temperature_RefPWF = GetOpacityTransferFunction('RefTemperature')
Temperature_RefPWF.Points = [0.43253340082384273, 0.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0]
Temperature_RefPWF.AllowDuplicateScalars = 1
Temperature_RefPWF.UseLogScale = 0
Temperature_RefPWF.ScalarRangeInitialized = 1

# get animation scene
animationScene1 = GetAnimationScene()

ReftemperatureLUTColorBar = GetScalarBar(Temperature_RefLUT, renderView1)
ReftemperatureLUTColorBar.AutoOrient = 0
ReftemperatureLUTColorBar.Orientation = 'Horizontal'
ReftemperatureLUTColorBar.WindowLocation = 'AnyLocation'
ReftemperatureLUTColorBar.Position = [0.3, 0.9]
ReftemperatureLUTColorBar.Title = 'Temperature'
ReftemperatureLUTColorBar.ComponentTitle = ''
ReftemperatureLUTColorBar.TitleJustification = 'Centered'
ReftemperatureLUTColorBar.HorizontalTitle = 0
ReftemperatureLUTColorBar.TitleOpacity = 1.0
ReftemperatureLUTColorBar.TitleFontFamily = 'Times'
ReftemperatureLUTColorBar.TitleFontFile = ''
ReftemperatureLUTColorBar.TitleBold = 0
ReftemperatureLUTColorBar.TitleItalic = 0
ReftemperatureLUTColorBar.TitleShadow = 0
ReftemperatureLUTColorBar.TitleFontSize = 16
ReftemperatureLUTColorBar.LabelOpacity = 1.0
ReftemperatureLUTColorBar.LabelFontFamily = 'Times'
ReftemperatureLUTColorBar.LabelFontFile = ''
ReftemperatureLUTColorBar.LabelBold = 0
ReftemperatureLUTColorBar.LabelItalic = 0
ReftemperatureLUTColorBar.LabelShadow = 0
ReftemperatureLUTColorBar.LabelFontSize = 16
ReftemperatureLUTColorBar.AutomaticLabelFormat = 0
ReftemperatureLUTColorBar.LabelFormat = '%4.2f'
ReftemperatureLUTColorBar.DrawTickMarks = 1
ReftemperatureLUTColorBar.DrawTickLabels = 1
ReftemperatureLUTColorBar.UseCustomLabels = 0
ReftemperatureLUTColorBar.CustomLabels = []
ReftemperatureLUTColorBar.AddRangeLabels = 1
ReftemperatureLUTColorBar.RangeLabelFormat = '%4.1f'
ReftemperatureLUTColorBar.DrawAnnotations = 1
ReftemperatureLUTColorBar.AddRangeAnnotations = 0
ReftemperatureLUTColorBar.AutomaticAnnotations = 0
ReftemperatureLUTColorBar.DrawNanAnnotation = 0
ReftemperatureLUTColorBar.NanAnnotation = 'NaN'
ReftemperatureLUTColorBar.TextPosition = 'Ticks right/top, annotations left/bottom'
ReftemperatureLUTColorBar.ReverseLegend = 0
ReftemperatureLUTColorBar.ScalarBarThickness = 20
ReftemperatureLUTColorBar.ScalarBarLength = 0.4


# updating the renderview
renderView1.Update()


# Exporting Reference file
SaveScreenshot('./Ref_Fin.png', renderView1, ImageResolution=[1673, 1164],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=1, 
    # PNG options
    CompressionLevel='5')

# hide data in view
Hide(opt_temperature_finpvd, renderView1)

# set active source
SetActiveSource(opt_temperature_intpvd)


# show data in view
opt_temperature_intpvdDisplay = Show(opt_temperature_intpvd, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
opt_temperature_intpvdDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

animationScene1.GoToFirst()

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(1673, 1164)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476

# Updating the renderview
renderView1.Update()

# save screenshot
for i in range(int(animationScene1.EndTime)+1):
    SaveScreenshot(os.path.join(init_fi_path ,str('./Rec_Int_%2.2i.png' %(i))), renderView1, ImageResolution=[1673, 1164],
        FontScaling='Scale fonts proportionally',
        OverrideColorPalette='',
        StereoMode='No change',
        TransparentBackground=1, 
        # PNG options
        CompressionLevel='5')
    animationScene1.GoToNext()



# set active source
SetActiveSource(opt_temperature_intpvd)

# get display properties
opt_temperature_intpvdDisplay = GetDisplayProperties(opt_temperature_intpvd, view=renderView1)

# set scalar coloring
ColorBy(opt_temperature_intpvdDisplay, ('POINTS', 'InitTemperature_Ref'))

#Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(temperatureLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
opt_temperature_intpvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
opt_temperature_intpvdDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Temperature'
initTemperature_RefLUT = GetColorTransferFunction('InitTemperature_Ref')
initTemperature_RefLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
initTemperature_RefLUT.InterpretValuesAsCategories = 0
initTemperature_RefLUT.AnnotationsInitialized = 0
initTemperature_RefLUT.ShowCategoricalColorsinDataRangeOnly = 0
initTemperature_RefLUT.RescaleOnVisibilityChange = 0
initTemperature_RefLUT.EnableOpacityMapping = 0
initTemperature_RefLUT.RGBPoints = [0.43253340082384273, 0.231373, 0.298039, 0.752941, 0.46626670041192136, 0.865003, 0.865003, 0.865003, 0.5, 0.705882, 0.0156863, 0.14902]
initTemperature_RefLUT.UseLogScale = 0
initTemperature_RefLUT.UseOpacityControlPointsFreehandDrawing = 0
initTemperature_RefLUT.ShowDataHistogram = 0
initTemperature_RefLUT.AutomaticDataHistogramComputation = 0
initTemperature_RefLUT.DataHistogramNumberOfBins = 10
initTemperature_RefLUT.ColorSpace = 'Diverging'
initTemperature_RefLUT.UseBelowRangeColor = 0
initTemperature_RefLUT.BelowRangeColor = [0.0, 0.0, 0.0]
initTemperature_RefLUT.UseAboveRangeColor = 0
initTemperature_RefLUT.AboveRangeColor = [0.5, 0.5, 0.5]
initTemperature_RefLUT.NanColor = [1.0, 1.0, 0.0]
initTemperature_RefLUT.NanOpacity = 1.0
initTemperature_RefLUT.Discretize = 1
initTemperature_RefLUT.NumberOfTableValues = 10 
initTemperature_RefLUT.ScalarRangeInitialized = 1.0
initTemperature_RefLUT.HSVWrap = 0
initTemperature_RefLUT.VectorComponent = 0
initTemperature_RefLUT.VectorMode = 'Magnitude'
initTemperature_RefLUT.AllowDuplicateScalars = 1
initTemperature_RefLUT.Annotations = []
initTemperature_RefLUT.ActiveAnnotatedValues = []
initTemperature_RefLUT.IndexedColors = []
initTemperature_RefLUT.IndexedOpacities = []

# get opacity transfer function/opacity map for 'Temperature'
initTemperature_RefPWF = GetOpacityTransferFunction('Temperature')
initTemperature_RefPWF.Points = [0.43253340082384273, 0.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0]
initTemperature_RefPWF.AllowDuplicateScalars = 1
initTemperature_RefPWF.UseLogScale = 0
initTemperature_RefPWF.ScalarRangeInitialized = 1

# get animation scene
animationScene1 = GetAnimationScene()

# Updating the renderview
renderView1.Update()

# Exporting Reference file
SaveScreenshot('./Ref_Int.png', renderView1, ImageResolution=[1673, 1164],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=1, 
    # PNG options
    CompressionLevel='5')


