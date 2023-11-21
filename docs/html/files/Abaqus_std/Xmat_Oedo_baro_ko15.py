# -*- coding: mbcs -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *


## Configurable Parameters
#
sample_height = 0.02  # [m]
sample_radius = 0.035 # [m]

oedo_pressure = [-10.0, -100.0, -50.0, -200.0] # [kPa]

steptime = 0.1
numoutputs = 2000

modelname = 'Xmat_Oedo_Abaqus'
userroutine = "xmat.f"

materialname = "Baro-Ko15_Hostun Sand 01"
voidratio = 0.7
density_material = 2.64/(voidratio+1) # [t/m^3]
materialparameters = (0.5899, 1.0, -2.5076, 1.0, 40.0, 0.91, 0.5)
statevariables = [voidratio, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


## Automated model creation
#
mdb.Model(name=modelname, modelType=STANDARD_EXPLICIT)
mymodel = mdb.models[modelname]

# Part: Element
#
sketch = mymodel.ConstrainedSketch(name='__profile__', sheetSize=1.0)
sketch.sketchOptions.setValues(viewStyle=AXISYM)
sketch.ConstructionLine(point1=(0.0, -0.5), point2=(0.0, 0.5))
sketch.FixedConstraint(entity=sketch.geometry.findAt((0.0, 0.0)))
sketch.rectangle(point1=(0.0, 0.0), point2=(sample_radius, sample_height))
sketch.ObliqueDimension(textPoint=(sample_radius/2, -sample_height/2), value=sample_radius,
   vertex1=sketch.vertices.findAt((0.0, 0.0)),
   vertex2=sketch.vertices.findAt((sample_radius, 0.0)))
sketch.ObliqueDimension(textPoint=(-sample_radius/2, sample_height/2), value=sample_height,
   vertex1=sketch.vertices.findAt((0.0, 0.0)),
   vertex2=sketch.vertices.findAt((0.0, sample_height)))
partElement = mymodel.Part(dimensionality=AXISYMMETRIC, name='Element', type=DEFORMABLE_BODY)
partElement.BaseShell(sketch=sketch)
del sketch

point_mu = partElement.vertices.findAt((0.0, sample_height, 0.0))
point_ou = partElement.vertices.findAt((sample_radius, sample_height, 0.0))
point_ml = partElement.vertices.findAt((0.0, 0.0, 0.0))
point_ol = partElement.vertices.findAt((sample_radius, 0.0, 0.0))

partElement.Set(name='setLowerEdge', vertices=[
   partElement.vertices[point_ml.index:point_ml.index+1] \
   + partElement.vertices[point_ol.index:point_ol.index+1]])
partElement.Set(name='setUpperEdge', vertices=[
   partElement.vertices[point_mu.index:point_mu.index+1] \
   + partElement.vertices[point_ou.index:point_ou.index+1]])
partElement.Set(faces=partElement.faces, name='setElement')

upperEdge = partElement.edges.findAt((sample_radius/2, sample_height, 0.0))
partElement.Surface(name='surfUpperEdge',
   side1Edges=partElement.edges[upperEdge.index:upperEdge.index+1])

# Material
#
mymodel.Material(name=materialname)
mymodel.materials[materialname].Density(table=((density_material, ), ))
mymodel.materials[materialname].Depvar(n=len(statevariables))
mymodel.materials[materialname].UserMaterial(mechanicalConstants=materialparameters)
mymodel.HomogeneousSolidSection(material=materialname, name='secElement', thickness=None)
partElement.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE,
   region=partElement.sets['setElement'], sectionName='secElement', thicknessAssignment=FROM_SECTION)

# Mesh
#
partElement.setElementType(elemTypes=(ElemType(elemCode=CAX4, elemLibrary=STANDARD),
   ElemType(elemCode=CAX3, elemLibrary=STANDARD)), regions=(partElement.sets['setElement']))
partElement.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=max(sample_height, sample_radius))
partElement.generateMesh()

# Assembly
#
mymodel.rootAssembly.DatumCsysByThreePoints(coordSysType=CYLINDRICAL,
   origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), point2=(0.0, 0.0, -1.0))
instElement = mymodel.rootAssembly.Instance(dependent=ON, name='instElement', part=partElement)

# Initial Conditions
#
mymodel.DisplacementBC(amplitude=UNSET, createStepName='Initial', distributionType=UNIFORM,
   fieldName='', localCsys=None, name='bcLowerEdge', region=instElement.sets['setLowerEdge'],
   u1=UNSET, u2=SET, ur3=UNSET)
mymodel.DisplacementBC(amplitude=UNSET, createStepName='Initial', distributionType=UNIFORM,
   fieldName='', localCsys=None, name='bcHorizontal', region=instElement.sets['setElement'],
   u1=SET, u2=UNSET, ur3=UNSET)

mymodel.GeostaticStress(lateralCoeff1=0.5, lateralCoeff2=None, name='StressState',
   region=instElement.sets['setElement'], stressMag1=oedo_pressure[0],
   stressMag2=oedo_pressure[0], vCoord1=0.0, vCoord2=sample_height)

# Step: InitialStressState
#
mymodel.GeostaticStep(name='InitialStressState', previous='Initial')
mymodel.steps['InitialStressState'].Restart(overlay=ON, timeMarks=OFF)
mymodel.fieldOutputRequests.changeKey(fromName='F-Output-1', toName='FieldOutput')
mymodel.fieldOutputRequests['FieldOutput'].setValues(variables=('U', 'E', 'S', 'SDV'))
mymodel.historyOutputRequests.changeKey(fromName='H-Output-1', toName='ElementHistoryOutput')
mymodel.historyOutputRequests['ElementHistoryOutput'].setValues(rebar=EXCLUDE,
   region=instElement.sets['setElement'], sectionPoints=DEFAULT,
   variables=('E11', 'E22', 'S11', 'S22', 'SDV'))
mymodel.HistoryOutputRequest(createStepName='InitialStressState', name='OutputUpperEdge',
   rebar=EXCLUDE, region=instElement.sets['setUpperEdge'], sectionPoints=DEFAULT,
   variables=('U1', 'U2'))
mymodel.Pressure(amplitude=UNSET, createStepName='InitialStressState', distributionType=UNIFORM,
   field='', name='VerticalLoad', magnitude=-oedo_pressure[0],
   region=instElement.surfaces['surfUpperEdge'])

# Step: LoadingStep1
#
thisstep = 'LoadingStep1'
mymodel.StaticStep(timePeriod=steptime, initialInc=0.0001, matrixSolver=DIRECT,
   matrixStorage=UNSYMMETRIC, maxNumInc=10000, name=thisstep, noStop=OFF, previous='InitialStressState',
   timeIncrementationMethod=AUTOMATIC)
mymodel.steps[thisstep].Restart(overlay=ON, timeMarks=OFF)
mymodel.fieldOutputRequests['FieldOutput'].setValuesInStep(stepName=thisstep,
   numIntervals=int(numoutputs/len(oedo_pressure)))
mymodel.historyOutputRequests['ElementHistoryOutput'].setValuesInStep(stepName=thisstep,
   numIntervals=int(numoutputs/len(oedo_pressure)))
mymodel.historyOutputRequests['OutputUpperEdge'].setValuesInStep(stepName=thisstep,
   numIntervals=int(numoutputs/len(oedo_pressure)))
mymodel.loads['VerticalLoad'].setValuesInStep(magnitude=-oedo_pressure[1], stepName=thisstep)

# Step: LoadingStep 2-n
#
for ilast in range(2, len(oedo_pressure)):
   laststep = thisstep
   thisstep = 'LoadingStep' + str(ilast)
   mymodel.StaticStep(timePeriod=steptime, initialInc=0.0001, matrixSolver=DIRECT,
      matrixStorage=UNSYMMETRIC, maxNumInc=10000, name=thisstep, noStop=OFF, previous=laststep,
      timeIncrementationMethod=AUTOMATIC)
   mymodel.steps[thisstep].Restart(overlay=ON, timeMarks=OFF)
   mymodel.loads['VerticalLoad'].setValuesInStep(magnitude=-oedo_pressure[ilast], stepName=thisstep)

# Job
#
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, explicitPrecision=SINGLE,
   getMemoryFromAnalysis=True, historyPrint=OFF, memory=90, memoryUnits=PERCENTAGE, model=modelname,
   modelPrint=OFF, multiprocessingMode=DEFAULT, name=modelname, nodalOutputPrecision=FULL,
   numCpus=2, numDomains=2, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=ANALYSIS,
   userSubroutine=userroutine, waitHours=0, waitMinutes=0)
myjob = mdb.jobs[modelname]

# Postprocessing the key word editor (last operation before writing the input file)
from math import ceil

output_block_size = 8
num_blocks = int(ceil(len(statevariables[(output_block_size-1):])/float(output_block_size)))
blocked_statevariables = [statevariables[(output_block_size*(n+1)-1):(output_block_size*(n+2)-1)] for n in range(num_blocks)]
output_statevar = ', ' + ', '.join([str(x) for x in statevariables[:(output_block_size-1)]])
for statevar_block in blocked_statevariables:
   output_statevar += ',\n' + ', '.join([str(x) for x in statevar_block])

mymodel.keywordBlock.synchVersions(storeNodesAndElements=False)
for idx, text in enumerate(mymodel.keywordBlock.sieBlocks):
   if (text[0:19] == '*Initial Conditions'):
      # The given statevariables are set as initial conditions
      # (Abaqus requires input files to have not more than eight values per line)
      newBlock = '**\n*Initial Conditions, type=SOLUTION\ninstElement.setElement' + output_statevar
      mymodel.keywordBlock.insert(idx, newBlock)

myjob.writeInput(consistencyChecking=OFF)
