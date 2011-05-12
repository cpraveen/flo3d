#!/usr/bin/env python
"""
Generates cut?.vtp and section?.dat files
cut?.vtp is a VTK XMLPolyData file -> Use paraview to visualize this
section?.dat is a plain ascii file which contains x, y, z, -Cp
"""

import sys
import math

#Check if a vtk file argument is given
if len(sys.argv) == 1:
   print "Please specify a vtk file name."
   print "Example:"
   print "  ", sys.argv[0], " vol_200.vtk"
   sys.exit(1)

#Set the vtk input file name
vtkfile = sys.argv[1]

import vtk

#Free stream values, q_inf = 0.5 * rho * (velocity)**2
mu      = 8.0e-4
rho_inf = 0.1
q_inf   = 86.797

nu      = mu/rho_inf

#Read the unstructured grid data
reader = vtk.vtkUnstructuredGridReader()
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.SetFileName(vtkfile)
reader.Update()

x = 0.8
Rex = q_inf * x / nu

 # Create three the line source to use for the probe lines.
line = vtk.vtkLineSource()
line.SetPoint1(x,0.0,0.0)
line.SetPoint2(x,0.5,0.0)
line.SetResolution(500)

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

f = open('x.dat','w')
for i in range(velocity.GetNumberOfTuples()):
   p  = data.GetPoint(i)
   a  = velocity.GetTuple3(i)
   y  = p[1]
   ux = a[0] / q_inf
   uy = a[1] * math.sqrt(2*Rex) / q_inf
   eta = y*math.sqrt(0.5*Rex)/x
   s  = str(eta) + " " + str(ux) + " " + str(uy) + "\n"
   f.write(s)
