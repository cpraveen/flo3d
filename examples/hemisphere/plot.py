#!/usr/bin/env python
"""
Generates cut.vtp and line.dat files
cut.vtp is a VTK XMLPolyData file -> Use paraview to visualize this
line.dat is a plain ascii file -> use gnuplot
"""

import sys

#Check if a vtk file argument is given
if len(sys.argv) == 1:
   print "Please specify a vtk file name."
   print "Example:"
   print "  ", sys.argv[0], " vol_200.vtk"
   sys.exit(1)

#Set the vtk input file name
vtkfile = sys.argv[1]

import vtk
import math

# Flow properties
gamma = 1.4
mach  = 3.0

#Normal to the cut plane
Nx      = 1
Ny      = 0
Nz      = 0

#Origin, some point on the cut plane
Ox      = 0
Oy      = 0
Oz      = 0

#Read the unstructured grid data
reader = vtk.vtkUnstructuredGridReader()
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.SetFileName(vtkfile)
reader.Update()
data = reader.GetOutput()

#Define cut plane normal to z-axis
p = vtk.vtkPlane()
p.SetNormal(Nx, Ny, Nz)
p.SetOrigin(Ox, Oy, Oz)

#Make cutter and apply it to data
cutter=vtk.vtkCutter()
cutter.SetCutFunction(p)
cutter.SetInput(data)
cutter.Update()
out = cutter.GetOutput()

#Save cut data to file
#Comment this if you dont want the cut-plane data
cfile = 'cut.vtp'
w = vtk.vtkXMLPolyDataWriter()
w.SetFileName(cfile)
w.SetInput(out)
w.Write()

# Extract solution along stagnation streamline

# Create three the line source to use for the probe lines.
line = vtk.vtkLineSource()
line.SetPoint1(0.0,0.0,2.0)
line.SetPoint2(0.0,0.0,1.0)
line.SetResolution(100)

# Move the line into place and create the probe filter.  For
# vtkProbeFilter, the probe line is the input, and the underlying data
# set is the source.
probe = vtk.vtkProbeFilter()
probe.SetInputConnection(line.GetOutputPort())
probe.SetSource(reader.GetOutput())
probe.Update()
data=probe.GetOutput()

#Extract velocity from point data
ptdata   = data.GetPointData()

arrayid  = ptdata.SetActiveVectors("velocity")
velocity = ptdata.GetArray(arrayid)

arrayid  = ptdata.SetActiveScalars("density")
density  = ptdata.GetArray(arrayid)

arrayid  = ptdata.SetActiveScalars("pressure")
pressure = ptdata.GetArray(arrayid)

# Writing following variables
# z, density, velocity, pressure, entropy
f = open('line.dat','w')
for i in range(velocity.GetNumberOfTuples()):
   pt = data.GetPoint(i)
   z  = pt[2]
   u  = velocity.GetTuple3(i)
   r  = density.GetTuple1(i)
   p  = pressure.GetTuple1(i)
   ent= math.log(p/r**gamma)
   s  = str(z) + " " + str(r) + " " + str(u[2]) + " " + str(p)
   s += " " + str(ent) + "\n"
   f.write(s)
