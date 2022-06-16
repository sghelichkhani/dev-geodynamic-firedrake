# trace generated using paraview version 5.9.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
import os

def __main__():
    final_fi_name = '/Users/sghelichkhani/Data/Remotes/kerguelen_data/Workplace/dev-geodynamic-firedrake/adjoint_tic/00_BoxExtruded_simple/01_low_ra_sink/inverse/visual_00_fslip/opt_temperature_fin.pvd' 
    initial_fi_name = '/Users/sghelichkhani/Data/Remotes/kerguelen_data/Workplace/dev-geodynamic-firedrake/adjoint_tic/00_BoxExtruded_simple/01_low_ra_sink/inverse/visual_00_fslip/opt_temperature_int.pvd'
    output_dir = './test/'

    visualise_recontruction(final_fi_name, initial_fi_name, output_dir)


def visualise_recontruction(final_fi_name, initial_fi_name, output_dir):

    
    
    my_image_resolution = [1514, 1228]
    my_parallel_scale = 0.7071067811865476
    my_z_camcoord =  2.7320508075688776
    number_of_table_values = 20
    
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    
    # create a new 'PVD Reader'
    opt_temperature_finpvd = PVDReader(registrationName='opt_temperature_fin.pvd', FileName=final_fi_name)
    opt_temperature_finpvd.CellArrays = []
    opt_temperature_finpvd.PointArrays = ['Temperature', 'final_state']
    opt_temperature_finpvd.ColumnArrays = []
    
    # get animation scene
    animationScene1 = GetAnimationScene()
    
    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()
    
    # create a new 'PVD Reader'
    opt_temperature_intpvd = PVDReader(registrationName='opt_temperature_int.pvd', FileName=initial_fi_name)
    opt_temperature_intpvd.CellArrays = []
    opt_temperature_intpvd.PointArrays = ['Temperature', 'ref_initial_state']
    opt_temperature_intpvd.ColumnArrays = []
    
    # set active source
    SetActiveSource(opt_temperature_finpvd)
    
    # set active source
    SetActiveSource(opt_temperature_finpvd)
    
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    
    # show data in view
    opt_temperature_finpvdDisplay = Show(opt_temperature_finpvd, renderView1, 'UnstructuredGridRepresentation')
    
    # get color transfer function/color map for 'Temperature'
    temperatureLUT = GetColorTransferFunction('Temperature')
    temperatureLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
    temperatureLUT.InterpretValuesAsCategories = 0
    temperatureLUT.AnnotationsInitialized = 0
    temperatureLUT.ShowCategoricalColorsinDataRangeOnly = 0
    temperatureLUT.RescaleOnVisibilityChange = 0
    temperatureLUT.EnableOpacityMapping = 0
    temperatureLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.5, 0.865003, 0.865003, 0.865003, 1.0, 0.705882, 0.0156863, 0.14902]
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
    temperatureLUT.NumberOfTableValues = number_of_table_values 
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
    temperaturePWF.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
    temperaturePWF.AllowDuplicateScalars = 1
    temperaturePWF.UseLogScale = 0
    temperaturePWF.ScalarRangeInitialized = 1
    
    # trace defaults for the display properties.
    opt_temperature_finpvdDisplay.Selection = None
    opt_temperature_finpvdDisplay.Representation = 'Surface'
    opt_temperature_finpvdDisplay.ColorArrayName = ['POINTS', 'Temperature']
    opt_temperature_finpvdDisplay.LookupTable = temperatureLUT
    opt_temperature_finpvdDisplay.MapScalars = 1
    opt_temperature_finpvdDisplay.MultiComponentsMapping = 0
    opt_temperature_finpvdDisplay.InterpolateScalarsBeforeMapping = 1
    opt_temperature_finpvdDisplay.Opacity = 1.0
    opt_temperature_finpvdDisplay.PointSize = 2.0
    opt_temperature_finpvdDisplay.LineWidth = 1.0
    opt_temperature_finpvdDisplay.RenderLinesAsTubes = 0
    opt_temperature_finpvdDisplay.RenderPointsAsSpheres = 0
    opt_temperature_finpvdDisplay.Interpolation = 'Gouraud'
    opt_temperature_finpvdDisplay.Specular = 0.0
    opt_temperature_finpvdDisplay.SpecularColor = [1.0, 1.0, 1.0]
    opt_temperature_finpvdDisplay.SpecularPower = 100.0
    opt_temperature_finpvdDisplay.Luminosity = 0.0
    opt_temperature_finpvdDisplay.Ambient = 0.0
    opt_temperature_finpvdDisplay.Diffuse = 1.0
    opt_temperature_finpvdDisplay.Roughness = 0.3
    opt_temperature_finpvdDisplay.Metallic = 0.0
    opt_temperature_finpvdDisplay.EdgeTint = [1.0, 1.0, 1.0]
    opt_temperature_finpvdDisplay.SelectTCoordArray = 'None'
    opt_temperature_finpvdDisplay.SelectNormalArray = 'None'
    opt_temperature_finpvdDisplay.SelectTangentArray = 'None'
    opt_temperature_finpvdDisplay.Texture = None
    opt_temperature_finpvdDisplay.RepeatTextures = 1
    opt_temperature_finpvdDisplay.InterpolateTextures = 0
    opt_temperature_finpvdDisplay.SeamlessU = 0
    opt_temperature_finpvdDisplay.SeamlessV = 0
    opt_temperature_finpvdDisplay.UseMipmapTextures = 0
    opt_temperature_finpvdDisplay.BaseColorTexture = None
    opt_temperature_finpvdDisplay.NormalTexture = None
    opt_temperature_finpvdDisplay.NormalScale = 1.0
    opt_temperature_finpvdDisplay.MaterialTexture = None
    opt_temperature_finpvdDisplay.OcclusionStrength = 1.0
    opt_temperature_finpvdDisplay.EmissiveTexture = None
    opt_temperature_finpvdDisplay.EmissiveFactor = [1.0, 1.0, 1.0]
    opt_temperature_finpvdDisplay.FlipTextures = 0
    opt_temperature_finpvdDisplay.BackfaceRepresentation = 'Follow Frontface'
    opt_temperature_finpvdDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
    opt_temperature_finpvdDisplay.BackfaceOpacity = 1.0
    opt_temperature_finpvdDisplay.Position = [0.0, 0.0, 0.0]
    opt_temperature_finpvdDisplay.Scale = [1.0, 1.0, 1.0]
    opt_temperature_finpvdDisplay.Orientation = [0.0, 0.0, 0.0]
    opt_temperature_finpvdDisplay.Origin = [0.0, 0.0, 0.0]
    opt_temperature_finpvdDisplay.CoordinateShiftScaleMethod = 'Always Auto Shift Scale'
    opt_temperature_finpvdDisplay.Pickable = 1
    opt_temperature_finpvdDisplay.Triangulate = 0
    opt_temperature_finpvdDisplay.UseShaderReplacements = 0
    opt_temperature_finpvdDisplay.ShaderReplacements = ''
    opt_temperature_finpvdDisplay.NonlinearSubdivisionLevel = 1
    opt_temperature_finpvdDisplay.UseDataPartitions = 0
    opt_temperature_finpvdDisplay.OSPRayUseScaleArray = 'All Approximate'
    opt_temperature_finpvdDisplay.OSPRayScaleArray = 'Temperature'
    opt_temperature_finpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    opt_temperature_finpvdDisplay.OSPRayMaterial = 'None'
    opt_temperature_finpvdDisplay.Orient = 0
    opt_temperature_finpvdDisplay.OrientationMode = 'Direction'
    opt_temperature_finpvdDisplay.SelectOrientationVectors = 'None'
    opt_temperature_finpvdDisplay.Scaling = 0
    opt_temperature_finpvdDisplay.ScaleMode = 'No Data Scaling Off'
    opt_temperature_finpvdDisplay.ScaleFactor = 0.1
    opt_temperature_finpvdDisplay.SelectScaleArray = 'Temperature'
    opt_temperature_finpvdDisplay.GlyphType = 'Arrow'
    opt_temperature_finpvdDisplay.UseGlyphTable = 0
    opt_temperature_finpvdDisplay.GlyphTableIndexArray = 'Temperature'
    opt_temperature_finpvdDisplay.UseCompositeGlyphTable = 0
    opt_temperature_finpvdDisplay.UseGlyphCullingAndLOD = 0
    opt_temperature_finpvdDisplay.LODValues = []
    opt_temperature_finpvdDisplay.ColorByLODIndex = 0
    opt_temperature_finpvdDisplay.GaussianRadius = 0.005
    opt_temperature_finpvdDisplay.ShaderPreset = 'Sphere'
    opt_temperature_finpvdDisplay.CustomTriangleScale = 3
    opt_temperature_finpvdDisplay.CustomShader = """ // This custom shader code define a gaussian blur
     // Please take a look into vtkSMPointGaussianRepresentation.cxx
     // for other custom shader examples
     //VTK::Color::Impl
       float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
       float gaussian = exp(-0.5*dist2);
       opacity = opacity*gaussian;
    """
    opt_temperature_finpvdDisplay.Emissive = 0
    opt_temperature_finpvdDisplay.ScaleByArray = 0
    opt_temperature_finpvdDisplay.SetScaleArray = ['POINTS', 'Temperature']
    opt_temperature_finpvdDisplay.ScaleArrayComponent = ''
    opt_temperature_finpvdDisplay.UseScaleFunction = 1
    opt_temperature_finpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    opt_temperature_finpvdDisplay.OpacityByArray = 0
    opt_temperature_finpvdDisplay.OpacityArray = ['POINTS', 'Temperature']
    opt_temperature_finpvdDisplay.OpacityArrayComponent = ''
    opt_temperature_finpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    opt_temperature_finpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
    opt_temperature_finpvdDisplay.SelectionCellLabelBold = 0
    opt_temperature_finpvdDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
    opt_temperature_finpvdDisplay.SelectionCellLabelFontFamily = 'Times'
    opt_temperature_finpvdDisplay.SelectionCellLabelFontFile = ''
    opt_temperature_finpvdDisplay.SelectionCellLabelFontSize = 18
    opt_temperature_finpvdDisplay.SelectionCellLabelItalic = 0
    opt_temperature_finpvdDisplay.SelectionCellLabelJustification = 'Left'
    opt_temperature_finpvdDisplay.SelectionCellLabelOpacity = 1.0
    opt_temperature_finpvdDisplay.SelectionCellLabelShadow = 0
    opt_temperature_finpvdDisplay.SelectionPointLabelBold = 0
    opt_temperature_finpvdDisplay.SelectionPointLabelColor = [1.0, 1.0, 0.0]
    opt_temperature_finpvdDisplay.SelectionPointLabelFontFamily = 'Times'
    opt_temperature_finpvdDisplay.SelectionPointLabelFontFile = ''
    opt_temperature_finpvdDisplay.SelectionPointLabelFontSize = 18
    opt_temperature_finpvdDisplay.SelectionPointLabelItalic = 0
    opt_temperature_finpvdDisplay.SelectionPointLabelJustification = 'Left'
    opt_temperature_finpvdDisplay.SelectionPointLabelOpacity = 1.0
    opt_temperature_finpvdDisplay.SelectionPointLabelShadow = 0
    opt_temperature_finpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
    opt_temperature_finpvdDisplay.ScalarOpacityFunction = temperaturePWF
    opt_temperature_finpvdDisplay.ScalarOpacityUnitDistance = 0.05009420535700122
    opt_temperature_finpvdDisplay.UseSeparateOpacityArray = 0
    opt_temperature_finpvdDisplay.OpacityArrayName = ['POINTS', 'Temperature']
    opt_temperature_finpvdDisplay.OpacityComponent = ''
    opt_temperature_finpvdDisplay.ExtractedBlockIndex = 0
    opt_temperature_finpvdDisplay.SelectMapper = 'Projected tetra'
    opt_temperature_finpvdDisplay.SamplingDimensions = [128, 128, 128]
    opt_temperature_finpvdDisplay.UseFloatingPointFrameBuffer = 1
    
    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    opt_temperature_finpvdDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
    opt_temperature_finpvdDisplay.OSPRayScaleFunction.UseLogScale = 0
    
    # init the 'Arrow' selected for 'GlyphType'
    opt_temperature_finpvdDisplay.GlyphType.TipResolution = 6
    opt_temperature_finpvdDisplay.GlyphType.TipRadius = 0.1
    opt_temperature_finpvdDisplay.GlyphType.TipLength = 0.35
    opt_temperature_finpvdDisplay.GlyphType.ShaftResolution = 6
    opt_temperature_finpvdDisplay.GlyphType.ShaftRadius = 0.03
    opt_temperature_finpvdDisplay.GlyphType.Invert = 0
    
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    opt_temperature_finpvdDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
    opt_temperature_finpvdDisplay.ScaleTransferFunction.UseLogScale = 0
    
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    opt_temperature_finpvdDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
    opt_temperature_finpvdDisplay.OpacityTransferFunction.UseLogScale = 0
    
    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    opt_temperature_finpvdDisplay.DataAxesGrid.XTitle = 'X Axis'
    opt_temperature_finpvdDisplay.DataAxesGrid.YTitle = 'Y Axis'
    opt_temperature_finpvdDisplay.DataAxesGrid.ZTitle = 'Z Axis'
    opt_temperature_finpvdDisplay.DataAxesGrid.XTitleFontFamily = 'Times'
    opt_temperature_finpvdDisplay.DataAxesGrid.XTitleFontFile = ''
    opt_temperature_finpvdDisplay.DataAxesGrid.XTitleBold = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.XTitleItalic = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.XTitleFontSize = 12
    opt_temperature_finpvdDisplay.DataAxesGrid.XTitleShadow = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.XTitleOpacity = 1.0
    opt_temperature_finpvdDisplay.DataAxesGrid.YTitleFontFamily = 'Times'
    opt_temperature_finpvdDisplay.DataAxesGrid.YTitleFontFile = ''
    opt_temperature_finpvdDisplay.DataAxesGrid.YTitleBold = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.YTitleItalic = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.YTitleFontSize = 12
    opt_temperature_finpvdDisplay.DataAxesGrid.YTitleShadow = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.YTitleOpacity = 1.0
    opt_temperature_finpvdDisplay.DataAxesGrid.ZTitleFontFamily = 'Times'
    opt_temperature_finpvdDisplay.DataAxesGrid.ZTitleFontFile = ''
    opt_temperature_finpvdDisplay.DataAxesGrid.ZTitleBold = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.ZTitleItalic = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.ZTitleFontSize = 12
    opt_temperature_finpvdDisplay.DataAxesGrid.ZTitleShadow = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.ZTitleOpacity = 1.0
    opt_temperature_finpvdDisplay.DataAxesGrid.FacesToRender = 63
    opt_temperature_finpvdDisplay.DataAxesGrid.CullBackface = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.CullFrontface = 1
    opt_temperature_finpvdDisplay.DataAxesGrid.ShowGrid = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.ShowEdges = 1
    opt_temperature_finpvdDisplay.DataAxesGrid.ShowTicks = 1
    opt_temperature_finpvdDisplay.DataAxesGrid.LabelUniqueEdgesOnly = 1
    opt_temperature_finpvdDisplay.DataAxesGrid.AxesToLabel = 63
    opt_temperature_finpvdDisplay.DataAxesGrid.XLabelFontFamily = 'Times'
    opt_temperature_finpvdDisplay.DataAxesGrid.XLabelFontFile = ''
    opt_temperature_finpvdDisplay.DataAxesGrid.XLabelBold = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.XLabelItalic = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.XLabelFontSize = 12
    opt_temperature_finpvdDisplay.DataAxesGrid.XLabelShadow = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.XLabelOpacity = 1.0
    opt_temperature_finpvdDisplay.DataAxesGrid.YLabelFontFamily = 'Times'
    opt_temperature_finpvdDisplay.DataAxesGrid.YLabelFontFile = ''
    opt_temperature_finpvdDisplay.DataAxesGrid.YLabelBold = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.YLabelItalic = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.YLabelFontSize = 12
    opt_temperature_finpvdDisplay.DataAxesGrid.YLabelShadow = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.YLabelOpacity = 1.0
    opt_temperature_finpvdDisplay.DataAxesGrid.ZLabelFontFamily = 'Times'
    opt_temperature_finpvdDisplay.DataAxesGrid.ZLabelFontFile = ''
    opt_temperature_finpvdDisplay.DataAxesGrid.ZLabelBold = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.ZLabelItalic = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.ZLabelFontSize = 12
    opt_temperature_finpvdDisplay.DataAxesGrid.ZLabelShadow = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.ZLabelOpacity = 1.0
    opt_temperature_finpvdDisplay.DataAxesGrid.XAxisNotation = 'Mixed'
    opt_temperature_finpvdDisplay.DataAxesGrid.XAxisPrecision = 2
    opt_temperature_finpvdDisplay.DataAxesGrid.XAxisUseCustomLabels = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.XAxisLabels = []
    opt_temperature_finpvdDisplay.DataAxesGrid.YAxisNotation = 'Mixed'
    opt_temperature_finpvdDisplay.DataAxesGrid.YAxisPrecision = 2
    opt_temperature_finpvdDisplay.DataAxesGrid.YAxisUseCustomLabels = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.YAxisLabels = []
    opt_temperature_finpvdDisplay.DataAxesGrid.ZAxisNotation = 'Mixed'
    opt_temperature_finpvdDisplay.DataAxesGrid.ZAxisPrecision = 2
    opt_temperature_finpvdDisplay.DataAxesGrid.ZAxisUseCustomLabels = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.ZAxisLabels = []
    opt_temperature_finpvdDisplay.DataAxesGrid.UseCustomBounds = 0
    opt_temperature_finpvdDisplay.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    
    # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
    opt_temperature_finpvdDisplay.PolarAxes.Visibility = 0
    opt_temperature_finpvdDisplay.PolarAxes.Translation = [0.0, 0.0, 0.0]
    opt_temperature_finpvdDisplay.PolarAxes.Scale = [1.0, 1.0, 1.0]
    opt_temperature_finpvdDisplay.PolarAxes.Orientation = [0.0, 0.0, 0.0]
    opt_temperature_finpvdDisplay.PolarAxes.EnableCustomBounds = [0, 0, 0]
    opt_temperature_finpvdDisplay.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    opt_temperature_finpvdDisplay.PolarAxes.EnableCustomRange = 0
    opt_temperature_finpvdDisplay.PolarAxes.CustomRange = [0.0, 1.0]
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisVisibility = 1
    opt_temperature_finpvdDisplay.PolarAxes.RadialAxesVisibility = 1
    opt_temperature_finpvdDisplay.PolarAxes.DrawRadialGridlines = 1
    opt_temperature_finpvdDisplay.PolarAxes.PolarArcsVisibility = 1
    opt_temperature_finpvdDisplay.PolarAxes.DrawPolarArcsGridlines = 1
    opt_temperature_finpvdDisplay.PolarAxes.NumberOfRadialAxes = 0
    opt_temperature_finpvdDisplay.PolarAxes.AutoSubdividePolarAxis = 1
    opt_temperature_finpvdDisplay.PolarAxes.NumberOfPolarAxis = 0
    opt_temperature_finpvdDisplay.PolarAxes.MinimumRadius = 0.0
    opt_temperature_finpvdDisplay.PolarAxes.MinimumAngle = 0.0
    opt_temperature_finpvdDisplay.PolarAxes.MaximumAngle = 90.0
    opt_temperature_finpvdDisplay.PolarAxes.RadialAxesOriginToPolarAxis = 1
    opt_temperature_finpvdDisplay.PolarAxes.Ratio = 1.0
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
    opt_temperature_finpvdDisplay.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
    opt_temperature_finpvdDisplay.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
    opt_temperature_finpvdDisplay.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
    opt_temperature_finpvdDisplay.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisTitleVisibility = 1
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisTitle = 'Radial Distance'
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisTitleLocation = 'Bottom'
    opt_temperature_finpvdDisplay.PolarAxes.PolarLabelVisibility = 1
    opt_temperature_finpvdDisplay.PolarAxes.PolarLabelFormat = '%-#6.3g'
    opt_temperature_finpvdDisplay.PolarAxes.PolarLabelExponentLocation = 'Labels'
    opt_temperature_finpvdDisplay.PolarAxes.RadialLabelVisibility = 1
    opt_temperature_finpvdDisplay.PolarAxes.RadialLabelFormat = '%-#3.1f'
    opt_temperature_finpvdDisplay.PolarAxes.RadialLabelLocation = 'Bottom'
    opt_temperature_finpvdDisplay.PolarAxes.RadialUnitsVisibility = 1
    opt_temperature_finpvdDisplay.PolarAxes.ScreenSize = 10.0
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisTitleOpacity = 1.0
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisTitleFontFamily = 'Times'
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisTitleBold = 0
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisTitleItalic = 0
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisTitleShadow = 0
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisTitleFontSize = 12
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisLabelOpacity = 1.0
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisLabelFontFamily = 'Times'
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisLabelBold = 0
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisLabelItalic = 0
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisLabelShadow = 0
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisLabelFontSize = 12
    opt_temperature_finpvdDisplay.PolarAxes.LastRadialAxisTextOpacity = 1.0
    opt_temperature_finpvdDisplay.PolarAxes.LastRadialAxisTextFontFamily = 'Times'
    opt_temperature_finpvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
    opt_temperature_finpvdDisplay.PolarAxes.LastRadialAxisTextBold = 0
    opt_temperature_finpvdDisplay.PolarAxes.LastRadialAxisTextItalic = 0
    opt_temperature_finpvdDisplay.PolarAxes.LastRadialAxisTextShadow = 0
    opt_temperature_finpvdDisplay.PolarAxes.LastRadialAxisTextFontSize = 12
    opt_temperature_finpvdDisplay.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
    opt_temperature_finpvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Times'
    opt_temperature_finpvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''
    opt_temperature_finpvdDisplay.PolarAxes.SecondaryRadialAxesTextBold = 0
    opt_temperature_finpvdDisplay.PolarAxes.SecondaryRadialAxesTextItalic = 0
    opt_temperature_finpvdDisplay.PolarAxes.SecondaryRadialAxesTextShadow = 0
    opt_temperature_finpvdDisplay.PolarAxes.SecondaryRadialAxesTextFontSize = 12
    opt_temperature_finpvdDisplay.PolarAxes.EnableDistanceLOD = 1
    opt_temperature_finpvdDisplay.PolarAxes.DistanceLODThreshold = 0.7
    opt_temperature_finpvdDisplay.PolarAxes.EnableViewAngleLOD = 1
    opt_temperature_finpvdDisplay.PolarAxes.ViewAngleLODThreshold = 0.7
    opt_temperature_finpvdDisplay.PolarAxes.SmallestVisiblePolarAngle = 0.5
    opt_temperature_finpvdDisplay.PolarAxes.PolarTicksVisibility = 1
    opt_temperature_finpvdDisplay.PolarAxes.ArcTicksOriginToPolarAxis = 1
    opt_temperature_finpvdDisplay.PolarAxes.TickLocation = 'Both'
    opt_temperature_finpvdDisplay.PolarAxes.AxisTickVisibility = 1
    opt_temperature_finpvdDisplay.PolarAxes.AxisMinorTickVisibility = 0
    opt_temperature_finpvdDisplay.PolarAxes.ArcTickVisibility = 1
    opt_temperature_finpvdDisplay.PolarAxes.ArcMinorTickVisibility = 0
    opt_temperature_finpvdDisplay.PolarAxes.DeltaAngleMajor = 10.0
    opt_temperature_finpvdDisplay.PolarAxes.DeltaAngleMinor = 5.0
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisMajorTickSize = 0.0
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisTickRatioSize = 0.3
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisMajorTickThickness = 1.0
    opt_temperature_finpvdDisplay.PolarAxes.PolarAxisTickRatioThickness = 0.5
    opt_temperature_finpvdDisplay.PolarAxes.LastRadialAxisMajorTickSize = 0.0
    opt_temperature_finpvdDisplay.PolarAxes.LastRadialAxisTickRatioSize = 0.3
    opt_temperature_finpvdDisplay.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
    opt_temperature_finpvdDisplay.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
    opt_temperature_finpvdDisplay.PolarAxes.ArcMajorTickSize = 0.0
    opt_temperature_finpvdDisplay.PolarAxes.ArcTickRatioSize = 0.3
    opt_temperature_finpvdDisplay.PolarAxes.ArcMajorTickThickness = 1.0
    opt_temperature_finpvdDisplay.PolarAxes.ArcTickRatioThickness = 0.5
    opt_temperature_finpvdDisplay.PolarAxes.Use2DMode = 0
    opt_temperature_finpvdDisplay.PolarAxes.UseLogAxis = 0
    
    # show color bar/color legend
    opt_temperature_finpvdDisplay.SetScalarBarVisibility(renderView1, True)
    
    # reset view to fit data
    renderView1.ResetCamera()
    
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
    temperatureLUTColorBar.ReverseLegend = 0
    temperatureLUTColorBar.ScalarBarThickness = 16
    temperatureLUTColorBar.ScalarBarLength = 0.33
    
    
    # Properties modified on temperatureLUTColorBar
    temperatureLUTColorBar.AutoOrient = 0
    temperatureLUTColorBar.Orientation = 'Horizontal'
    temperatureLUTColorBar.WindowLocation = 'AnyLocation'
    temperatureLUTColorBar.Position = [0.3, 0.86]
    temperatureLUTColorBar.Title = 'Temperature [K]'
    temperatureLUTColorBar.TitleFontFamily = 'Times'
    temperatureLUTColorBar.LabelFontFamily = 'Times'
    temperatureLUTColorBar.ReverseLegend = 1
    temperatureLUTColorBar.ScalarBarLength = 0.4
    
    # Properties modified on temperatureLUTColorBar
    temperatureLUTColorBar.Position = [0.3, 0.9]
    
    # Properties modified on temperatureLUTColorBar
    temperatureLUTColorBar.Position = [0.3, 0.86]
    
    # Hide orientation axes
    renderView1.OrientationAxesVisibility = 0
    
    # get the material library
    materialLibrary1 = GetMaterialLibrary()
    
    # Properties modified on temperatureLUT
    temperatureLUT.NumberOfTableValues = number_of_table_values
    
    # Properties modified on temperatureLUT
    temperatureLUT.NumberOfTableValues = number_of_table_values
    
    # get layout
    layout1 = GetLayout()
    
    # layout/tab size in pixels
    layout1.SetSize(1514, 1228)
    
    # current camera placement for renderView1
    renderView1.CameraPosition = [0.5, 0.5, my_z_camcoord]
    renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
    renderView1.CameraParallelScale = my_parallel_scale 
    
    # save animation
    SaveAnimation(os.path.join(output_dir,'reconstruction_final.png'), renderView1, ImageResolution=my_image_resolution,
        FontScaling='Scale fonts proportionally',
        OverrideColorPalette='',
        StereoMode='No change',
        TransparentBackground=1,
        FrameRate=1,
        FrameWindow=[0, 50], 
        # PNG options
        CompressionLevel='5',
        SuffixFormat='_%4.4i')
    
    # set scalar coloring
    ColorBy(opt_temperature_finpvdDisplay, ('POINTS', 'final_state'))
    
    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(temperatureLUT, renderView1)
    
    # rescale color and/or opacity maps used to include current data range
    opt_temperature_finpvdDisplay.RescaleTransferFunctionToDataRange(True, False)
    
    # show color bar/color legend
    opt_temperature_finpvdDisplay.SetScalarBarVisibility(renderView1, True)
    
    # get color transfer function/color map for 'final_state'
    final_stateLUT = GetColorTransferFunction('final_state')
    final_stateLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
    final_stateLUT.InterpretValuesAsCategories = 0
    final_stateLUT.AnnotationsInitialized = 0
    final_stateLUT.ShowCategoricalColorsinDataRangeOnly = 0
    final_stateLUT.RescaleOnVisibilityChange = 0
    final_stateLUT.EnableOpacityMapping = 0
    final_stateLUT.RGBPoints = [1.9319120686716416e-19, 0.231373, 0.298039, 0.752941, 0.5000000000000006, 0.865003, 0.865003, 0.865003, 1.000000000000001, 0.705882, 0.0156863, 0.14902]
    final_stateLUT.UseLogScale = 0
    final_stateLUT.UseOpacityControlPointsFreehandDrawing = 0
    final_stateLUT.ShowDataHistogram = 0
    final_stateLUT.AutomaticDataHistogramComputation = 0
    final_stateLUT.DataHistogramNumberOfBins = 10
    final_stateLUT.ColorSpace = 'Diverging'
    final_stateLUT.UseBelowRangeColor = 0
    final_stateLUT.BelowRangeColor = [0.0, 0.0, 0.0]
    final_stateLUT.UseAboveRangeColor = 0
    final_stateLUT.AboveRangeColor = [0.5, 0.5, 0.5]
    final_stateLUT.NanColor = [1.0, 1.0, 0.0]
    final_stateLUT.NanOpacity = 1.0
    final_stateLUT.Discretize = 1
    final_stateLUT.NumberOfTableValues = number_of_table_values
    final_stateLUT.ScalarRangeInitialized = 1.0
    final_stateLUT.HSVWrap = 0
    final_stateLUT.VectorComponent = 0
    final_stateLUT.VectorMode = 'Magnitude'
    final_stateLUT.AllowDuplicateScalars = 1
    final_stateLUT.Annotations = []
    final_stateLUT.ActiveAnnotatedValues = []
    final_stateLUT.IndexedColors = []
    final_stateLUT.IndexedOpacities = []
    
    # get opacity transfer function/opacity map for 'final_state'
    final_statePWF = GetOpacityTransferFunction('final_state')
    final_statePWF.Points = [1.9319120686716416e-19, 0.0, 0.5, 0.0, 1.000000000000001, 1.0, 0.5, 0.0]
    final_statePWF.AllowDuplicateScalars = 1
    final_statePWF.UseLogScale = 0
    final_statePWF.ScalarRangeInitialized = 1
    
    # get color legend/bar for final_stateLUT in view renderView1
    final_stateLUTColorBar = GetScalarBar(final_stateLUT, renderView1)
    final_stateLUTColorBar.AutoOrient = 1
    final_stateLUTColorBar.Orientation = 'Vertical'
    final_stateLUTColorBar.WindowLocation = 'LowerRightCorner'
    final_stateLUTColorBar.Position = [0.89, 0.02]
    final_stateLUTColorBar.Title = 'final_state'
    final_stateLUTColorBar.ComponentTitle = ''
    final_stateLUTColorBar.TitleJustification = 'Centered'
    final_stateLUTColorBar.HorizontalTitle = 0
    final_stateLUTColorBar.TitleOpacity = 1.0
    final_stateLUTColorBar.TitleFontFamily = 'Times'
    final_stateLUTColorBar.TitleFontFile = ''
    final_stateLUTColorBar.TitleBold = 0
    final_stateLUTColorBar.TitleItalic = 0
    final_stateLUTColorBar.TitleShadow = 0
    final_stateLUTColorBar.TitleFontSize = 16
    final_stateLUTColorBar.LabelOpacity = 1.0
    final_stateLUTColorBar.LabelFontFamily = 'Times'
    final_stateLUTColorBar.LabelFontFile = ''
    final_stateLUTColorBar.LabelBold = 0
    final_stateLUTColorBar.LabelItalic = 0
    final_stateLUTColorBar.LabelShadow = 0
    final_stateLUTColorBar.LabelFontSize = 16
    final_stateLUTColorBar.AutomaticLabelFormat = 1
    final_stateLUTColorBar.LabelFormat = '%-#6.3g'
    final_stateLUTColorBar.DrawTickMarks = 1
    final_stateLUTColorBar.DrawTickLabels = 1
    final_stateLUTColorBar.UseCustomLabels = 0
    final_stateLUTColorBar.CustomLabels = []
    final_stateLUTColorBar.AddRangeLabels = 1
    final_stateLUTColorBar.RangeLabelFormat = '%-#6.1e'
    final_stateLUTColorBar.DrawAnnotations = 1
    final_stateLUTColorBar.AddRangeAnnotations = 0
    final_stateLUTColorBar.AutomaticAnnotations = 0
    final_stateLUTColorBar.DrawNanAnnotation = 0
    final_stateLUTColorBar.NanAnnotation = 'NaN'
    final_stateLUTColorBar.TextPosition = 'Ticks right/top, annotations left/bottom'
    final_stateLUTColorBar.ReverseLegend = 0
    final_stateLUTColorBar.ScalarBarThickness = 16
    final_stateLUTColorBar.ScalarBarLength = 0.33
    
    # Properties modified on final_stateLUTColorBar
    final_stateLUTColorBar.AutoOrient = 0
    final_stateLUTColorBar.Orientation = 'Horizontal'
    final_stateLUTColorBar.WindowLocation = 'AnyLocation'
    final_stateLUTColorBar.Position = [0.3, 0.86]
    final_stateLUTColorBar.TitleFontFamily = 'Times'
    final_stateLUTColorBar.LabelFontFamily = 'Times'
    final_stateLUTColorBar.ReverseLegend = 1
    final_stateLUTColorBar.ScalarBarLength = 0.4
    
    # Properties modified on final_stateLUTColorBar
    final_stateLUTColorBar.Title = 'Temperature [K]'
    
    # Properties modified on final_stateLUT
    final_stateLUT.NumberOfTableValues = number_of_table_values
    
    # layout/tab size in pixels
    layout1.SetSize(1514, 1228)
    
    # current camera placement for renderView1
    renderView1.CameraPosition = [0.5, 0.5, my_z_camcoord]
    renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
    renderView1.CameraParallelScale = my_parallel_scale 
    
    # save screenshot
    SaveScreenshot(os.path.join(output_dir, 'reference_final.png'), renderView1, ImageResolution=my_image_resolution,
        FontScaling='Scale fonts proportionally',
        OverrideColorPalette='',
        StereoMode='No change',
        TransparentBackground=1, 
        # PNG options
        CompressionLevel='5')
    
    # set active source
    SetActiveSource(opt_temperature_intpvd)
    
    # set active source
    SetActiveSource(opt_temperature_intpvd)
    
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
    opt_temperature_intpvdDisplay.SelectionCellLabelFontFamily = 'Times'
    opt_temperature_intpvdDisplay.SelectionCellLabelFontFile = ''
    opt_temperature_intpvdDisplay.SelectionCellLabelFontSize = 18
    opt_temperature_intpvdDisplay.SelectionCellLabelItalic = 0
    opt_temperature_intpvdDisplay.SelectionCellLabelJustification = 'Left'
    opt_temperature_intpvdDisplay.SelectionCellLabelOpacity = 1.0
    opt_temperature_intpvdDisplay.SelectionCellLabelShadow = 0
    opt_temperature_intpvdDisplay.SelectionPointLabelBold = 0
    opt_temperature_intpvdDisplay.SelectionPointLabelColor = [1.0, 1.0, 0.0]
    opt_temperature_intpvdDisplay.SelectionPointLabelFontFamily = 'Times'
    opt_temperature_intpvdDisplay.SelectionPointLabelFontFile = ''
    opt_temperature_intpvdDisplay.SelectionPointLabelFontSize = 18
    opt_temperature_intpvdDisplay.SelectionPointLabelItalic = 0
    opt_temperature_intpvdDisplay.SelectionPointLabelJustification = 'Left'
    opt_temperature_intpvdDisplay.SelectionPointLabelOpacity = 1.0
    opt_temperature_intpvdDisplay.SelectionPointLabelShadow = 0
    opt_temperature_intpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
    opt_temperature_intpvdDisplay.ScalarOpacityFunction = temperaturePWF
    opt_temperature_intpvdDisplay.ScalarOpacityUnitDistance = 0.05009420535700122
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
    opt_temperature_intpvdDisplay.ScaleTransferFunction.Points = [1.9319120686716416e-19, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
    opt_temperature_intpvdDisplay.ScaleTransferFunction.UseLogScale = 0
    
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    opt_temperature_intpvdDisplay.OpacityTransferFunction.Points = [1.9319120686716416e-19, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
    opt_temperature_intpvdDisplay.OpacityTransferFunction.UseLogScale = 0
    
    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    opt_temperature_intpvdDisplay.DataAxesGrid.XTitle = 'X Axis'
    opt_temperature_intpvdDisplay.DataAxesGrid.YTitle = 'Y Axis'
    opt_temperature_intpvdDisplay.DataAxesGrid.ZTitle = 'Z Axis'
    opt_temperature_intpvdDisplay.DataAxesGrid.XTitleFontFamily = 'Times'
    opt_temperature_intpvdDisplay.DataAxesGrid.XTitleFontFile = ''
    opt_temperature_intpvdDisplay.DataAxesGrid.XTitleBold = 0
    opt_temperature_intpvdDisplay.DataAxesGrid.XTitleItalic = 0
    opt_temperature_intpvdDisplay.DataAxesGrid.XTitleFontSize = 12
    opt_temperature_intpvdDisplay.DataAxesGrid.XTitleShadow = 0
    opt_temperature_intpvdDisplay.DataAxesGrid.XTitleOpacity = 1.0
    opt_temperature_intpvdDisplay.DataAxesGrid.YTitleFontFamily = 'Times'
    opt_temperature_intpvdDisplay.DataAxesGrid.YTitleFontFile = ''
    opt_temperature_intpvdDisplay.DataAxesGrid.YTitleBold = 0
    opt_temperature_intpvdDisplay.DataAxesGrid.YTitleItalic = 0
    opt_temperature_intpvdDisplay.DataAxesGrid.YTitleFontSize = 12
    opt_temperature_intpvdDisplay.DataAxesGrid.YTitleShadow = 0
    opt_temperature_intpvdDisplay.DataAxesGrid.YTitleOpacity = 1.0
    opt_temperature_intpvdDisplay.DataAxesGrid.ZTitleFontFamily = 'Times'
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
    opt_temperature_intpvdDisplay.DataAxesGrid.XLabelFontFamily = 'Times'
    opt_temperature_intpvdDisplay.DataAxesGrid.XLabelFontFile = ''
    opt_temperature_intpvdDisplay.DataAxesGrid.XLabelBold = 0
    opt_temperature_intpvdDisplay.DataAxesGrid.XLabelItalic = 0
    opt_temperature_intpvdDisplay.DataAxesGrid.XLabelFontSize = 12
    opt_temperature_intpvdDisplay.DataAxesGrid.XLabelShadow = 0
    opt_temperature_intpvdDisplay.DataAxesGrid.XLabelOpacity = 1.0
    opt_temperature_intpvdDisplay.DataAxesGrid.YLabelFontFamily = 'Times'
    opt_temperature_intpvdDisplay.DataAxesGrid.YLabelFontFile = ''
    opt_temperature_intpvdDisplay.DataAxesGrid.YLabelBold = 0
    opt_temperature_intpvdDisplay.DataAxesGrid.YLabelItalic = 0
    opt_temperature_intpvdDisplay.DataAxesGrid.YLabelFontSize = 12
    opt_temperature_intpvdDisplay.DataAxesGrid.YLabelShadow = 0
    opt_temperature_intpvdDisplay.DataAxesGrid.YLabelOpacity = 1.0
    opt_temperature_intpvdDisplay.DataAxesGrid.ZLabelFontFamily = 'Times'
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
    opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTitleFontFamily = 'Times'
    opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
    opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTitleBold = 0
    opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTitleItalic = 0
    opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTitleShadow = 0
    opt_temperature_intpvdDisplay.PolarAxes.PolarAxisTitleFontSize = 12
    opt_temperature_intpvdDisplay.PolarAxes.PolarAxisLabelOpacity = 1.0
    opt_temperature_intpvdDisplay.PolarAxes.PolarAxisLabelFontFamily = 'Times'
    opt_temperature_intpvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
    opt_temperature_intpvdDisplay.PolarAxes.PolarAxisLabelBold = 0
    opt_temperature_intpvdDisplay.PolarAxes.PolarAxisLabelItalic = 0
    opt_temperature_intpvdDisplay.PolarAxes.PolarAxisLabelShadow = 0
    opt_temperature_intpvdDisplay.PolarAxes.PolarAxisLabelFontSize = 12
    opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisTextOpacity = 1.0
    opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisTextFontFamily = 'Times'
    opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
    opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisTextBold = 0
    opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisTextItalic = 0
    opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisTextShadow = 0
    opt_temperature_intpvdDisplay.PolarAxes.LastRadialAxisTextFontSize = 12
    opt_temperature_intpvdDisplay.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
    opt_temperature_intpvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Times'
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
    
    # hide data in view
    Hide(opt_temperature_finpvd, renderView1)
    
    # Properties modified on temperatureLUTColorBar
    temperatureLUTColorBar.WindowLocation = 'AnyLocation'
    
    # change scalar bar placement
    temperatureLUTColorBar.Position = [0.5561426684280053, 0.8925692157558975]
    
    # change scalar bar placement
    temperatureLUTColorBar.Position = [0.3091149273447821, 0.8974552092412396]
    
    # Properties modified on temperatureLUTColorBar
    temperatureLUTColorBar.Position = [0.3, 0.86]
    
    # layout/tab size in pixels
    layout1.SetSize(1514, 1228)
    
    # current camera placement for renderView1
    renderView1.CameraPosition = [0.5, 0.5, my_z_camcoord]
    renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
    renderView1.CameraParallelScale = my_parallel_scale 
    
    # save animation
    SaveAnimation(os.path.join(output_dir, 'reconstruction_initial.png'), renderView1, ImageResolution=my_image_resolution,
        FontScaling='Scale fonts proportionally',
        OverrideColorPalette='',
        StereoMode='No change',
        TransparentBackground=1,
        FrameRate=1,
        FrameWindow=[0, 50], 
        # PNG options
        CompressionLevel='5',
        SuffixFormat='_%4.4i')
    
    # set active source
    SetActiveSource(opt_temperature_intpvd)
    
    # set scalar coloring
    ColorBy(opt_temperature_intpvdDisplay, ('POINTS', 'ref_initial_state'))
    
    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(temperatureLUT, renderView1)
    
    # rescale color and/or opacity maps used to include current data range
    opt_temperature_intpvdDisplay.RescaleTransferFunctionToDataRange(True, False)
    
    # show color bar/color legend
    opt_temperature_intpvdDisplay.SetScalarBarVisibility(renderView1, True)
    
    # get color transfer function/color map for 'ref_initial_state'
    ref_initial_stateLUT = GetColorTransferFunction('ref_initial_state')
    ref_initial_stateLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
    ref_initial_stateLUT.InterpretValuesAsCategories = 0
    ref_initial_stateLUT.AnnotationsInitialized = 0
    ref_initial_stateLUT.ShowCategoricalColorsinDataRangeOnly = 0
    ref_initial_stateLUT.RescaleOnVisibilityChange = 0
    ref_initial_stateLUT.EnableOpacityMapping = 0
    ref_initial_stateLUT.RGBPoints = [1.4914345609058255e-18, 0.231373, 0.298039, 0.752941, 0.5000000000000006, 0.865003, 0.865003, 0.865003, 1.000000000000001, 0.705882, 0.0156863, 0.14902]
    ref_initial_stateLUT.UseLogScale = 0
    ref_initial_stateLUT.UseOpacityControlPointsFreehandDrawing = 0
    ref_initial_stateLUT.ShowDataHistogram = 0
    ref_initial_stateLUT.AutomaticDataHistogramComputation = 0
    ref_initial_stateLUT.DataHistogramNumberOfBins = 10
    ref_initial_stateLUT.ColorSpace = 'Diverging'
    ref_initial_stateLUT.UseBelowRangeColor = 0
    ref_initial_stateLUT.BelowRangeColor = [0.0, 0.0, 0.0]
    ref_initial_stateLUT.UseAboveRangeColor = 0
    ref_initial_stateLUT.AboveRangeColor = [0.5, 0.5, 0.5]
    ref_initial_stateLUT.NanColor = [1.0, 1.0, 0.0]
    ref_initial_stateLUT.NanOpacity = 1.0
    ref_initial_stateLUT.Discretize = 1
    ref_initial_stateLUT.NumberOfTableValues = number_of_table_values
    ref_initial_stateLUT.ScalarRangeInitialized = 1.0
    ref_initial_stateLUT.HSVWrap = 0
    ref_initial_stateLUT.VectorComponent = 0
    ref_initial_stateLUT.VectorMode = 'Magnitude'
    ref_initial_stateLUT.AllowDuplicateScalars = 1
    ref_initial_stateLUT.Annotations = []
    ref_initial_stateLUT.ActiveAnnotatedValues = []
    ref_initial_stateLUT.IndexedColors = []
    ref_initial_stateLUT.IndexedOpacities = []
    
    # get opacity transfer function/opacity map for 'ref_initial_state'
    ref_initial_statePWF = GetOpacityTransferFunction('ref_initial_state')
    ref_initial_statePWF.Points = [1.4914345609058255e-18, 0.0, 0.5, 0.0, 1.000000000000001, 1.0, 0.5, 0.0]
    ref_initial_statePWF.AllowDuplicateScalars = 1
    ref_initial_statePWF.UseLogScale = 0
    ref_initial_statePWF.ScalarRangeInitialized = 1
    
    # get color legend/bar for ref_initial_stateLUT in view renderView1
    ref_initial_stateLUTColorBar = GetScalarBar(ref_initial_stateLUT, renderView1)
    ref_initial_stateLUTColorBar.AutoOrient = 1
    ref_initial_stateLUTColorBar.Orientation = 'Vertical'
    ref_initial_stateLUTColorBar.WindowLocation = 'LowerRightCorner'
    ref_initial_stateLUTColorBar.Position = [0.89, 0.02]
    ref_initial_stateLUTColorBar.Title = 'ref_initial_state'
    ref_initial_stateLUTColorBar.ComponentTitle = ''
    ref_initial_stateLUTColorBar.TitleJustification = 'Centered'
    ref_initial_stateLUTColorBar.HorizontalTitle = 0
    ref_initial_stateLUTColorBar.TitleOpacity = 1.0
    ref_initial_stateLUTColorBar.TitleFontFamily = 'Times'
    ref_initial_stateLUTColorBar.TitleFontFile = ''
    ref_initial_stateLUTColorBar.TitleBold = 0
    ref_initial_stateLUTColorBar.TitleItalic = 0
    ref_initial_stateLUTColorBar.TitleShadow = 0
    ref_initial_stateLUTColorBar.TitleFontSize = 16
    ref_initial_stateLUTColorBar.LabelOpacity = 1.0
    ref_initial_stateLUTColorBar.LabelFontFamily = 'Times'
    ref_initial_stateLUTColorBar.LabelFontFile = ''
    ref_initial_stateLUTColorBar.LabelBold = 0
    ref_initial_stateLUTColorBar.LabelItalic = 0
    ref_initial_stateLUTColorBar.LabelShadow = 0
    ref_initial_stateLUTColorBar.LabelFontSize = 16
    ref_initial_stateLUTColorBar.AutomaticLabelFormat = 1
    ref_initial_stateLUTColorBar.LabelFormat = '%-#6.3g'
    ref_initial_stateLUTColorBar.DrawTickMarks = 1
    ref_initial_stateLUTColorBar.DrawTickLabels = 1
    ref_initial_stateLUTColorBar.UseCustomLabels = 0
    ref_initial_stateLUTColorBar.CustomLabels = []
    ref_initial_stateLUTColorBar.AddRangeLabels = 1
    ref_initial_stateLUTColorBar.RangeLabelFormat = '%-#6.1e'
    ref_initial_stateLUTColorBar.DrawAnnotations = 1
    ref_initial_stateLUTColorBar.AddRangeAnnotations = 0
    ref_initial_stateLUTColorBar.AutomaticAnnotations = 0
    ref_initial_stateLUTColorBar.DrawNanAnnotation = 0
    ref_initial_stateLUTColorBar.NanAnnotation = 'NaN'
    ref_initial_stateLUTColorBar.TextPosition = 'Ticks right/top, annotations left/bottom'
    ref_initial_stateLUTColorBar.ReverseLegend = 0
    ref_initial_stateLUTColorBar.ScalarBarThickness = 16
    ref_initial_stateLUTColorBar.ScalarBarLength = 0.33
    
    # Properties modified on ref_initial_stateLUTColorBar
    ref_initial_stateLUTColorBar.AutoOrient = 0
    ref_initial_stateLUTColorBar.Orientation = 'Horizontal'
    ref_initial_stateLUTColorBar.WindowLocation = 'AnyLocation'
    ref_initial_stateLUTColorBar.Position = [0.3, 0.86]
    ref_initial_stateLUTColorBar.Title = 'Temperature [K]'
    ref_initial_stateLUTColorBar.TitleFontFamily = 'Times'
    ref_initial_stateLUTColorBar.LabelFontFamily = 'Times'
    ref_initial_stateLUTColorBar.ReverseLegend = 1
    ref_initial_stateLUTColorBar.ScalarBarLength = 0.4
    
    
    # Properties modified on ref_initial_stateLUT
    ref_initial_stateLUT.NumberOfTableValues = number_of_table_values 
    
    # layout/tab size in pixels
    layout1.SetSize(1514, 1228)
    
    # current camera placement for renderView1
    renderView1.CameraPosition = [0.5, 0.5, my_z_camcoord]
    renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
    renderView1.CameraParallelScale = my_parallel_scale
    
    # save screenshot
    SaveScreenshot(os.path.join(output_dir, 'reference_initial.png'), renderView1, ImageResolution=my_image_resolution,
        FontScaling='Scale fonts proportionally',
        OverrideColorPalette='',
        StereoMode='No change',
        TransparentBackground=1, 
        # PNG options
        CompressionLevel='5')


if __name__=="__main__":
    __main__()
