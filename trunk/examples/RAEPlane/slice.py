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
GAMMA   = 1.4
Mach    = 0.9
pre_inf = 29765
q_inf   = 16650

#Bounding box containing the geometry
#should not contain outer boundary
xmin, xmax = -0.1, 2.00
ymin, ymax = -0.08, 0.08
zmin, zmax = -0.5, 0.5

deg = 15.0 * math.pi / 180.0

#Normal to the cut plane
Nx = [0,  0,             0, 0, 0, 0, 0, 0]
Ny = [0, -math.cos(deg), 0, 0, 0, 0, 0, 0]
Nz = [1,  math.sin(deg), 1, 1, 1, 1, 1, 1]

#Origin, some point on the cut plane
Ox      = 0
Oy      = 0

#List of cut planes to extract
#these are section for which experimental data is available
span  = 0.457
zcuts = [0, 0, 0.25*span, 0.4*span, 0.6*span, 0.75*span, 0.85*span, 0.925*span]

#Fine cell/line which has p as first point
def NextCell(p,d):
   ncells = d.GetNumberOfCells()
   for i in range(0,ncells):
      cell = d.GetCell(i)
      p1   = cell.GetPointId(0)
      p2   = cell.GetPointId(1)
      if p == p1:
         return i

#Read the unstructured grid data
reader = vtk.vtkUnstructuredGridReader()
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.SetFileName(vtkfile)
reader.Update()
data = reader.GetOutput()

#Loop over cut-planes
ncut = 0
for Oz in zcuts:
   print "Extracting data for z =", Oz

   #Define cut plane normal to z-axis
   p = vtk.vtkPlane()
   p.SetNormal(Nx[ncut], Ny[ncut], Nz[ncut])
   p.SetOrigin(Ox, Oy, Oz)

   #Make cutter and apply it to data
   cutter=vtk.vtkCutter()
   cutter.SetCutFunction(p)
   cutter.SetInput(data)
   cutter.Update()
   out = cutter.GetOutput()

   #Save cut data to file
   #Comment this if you dont want the cut-plane data
   cfile = 'cut' + str(int(ncut)) + '.vtp'
   w = vtk.vtkXMLPolyDataWriter()
   w.SetFileName(cfile)
   w.SetInput(out)
   w.Write()

   #Extract boundary, includes outer boundary
   bd = vtk.vtkFeatureEdges()
   bd.BoundaryEdgesOn()
   bd.ColoringOff()
   bd.FeatureEdgesOff()
   bd.ManifoldEdgesOff()
   bd.AddInput(out)
   bd.Update()
   shape = bd.GetOutput()
   if shape.GetNumberOfPoints() == 0:
      print "No points found, something wrong !!!"
   
   # Now extract only the shape not the outer boundary
   g = vtk.vtkGeometryFilter()
   g.ExtentClippingOn()
   g.SetExtent(xmin,xmax,ymin,ymax,zmin,zmax)
   g.SetInput(shape)
   g.Update()
   geom = g.GetOutput()

   #Clean geom to remove unused points
   #Not necessary, can be removed I think
   cleaner = vtk.vtkCleanPolyData()
   cleaner.SetInput(geom)
   cleaner.Update()
   geom = cleaner.GetOutput()
   if geom.GetNumberOfPoints() == 0:
      print "No points found, something wrong !!!"

   #Extract pressure from point data
   ptdata   = geom.GetPointData()
   arrayid  = ptdata.SetActiveScalars("pressure")
   pressure = ptdata.GetArray(arrayid)

   #Write shape points into file
   sfile = 'section' + str(int(ncut)) + '.dat'
   numCells = geom.GetNumberOfCells()
   f = open(sfile, 'w')
   cellid = 0
   for i in range(0,numCells+1):
      cell = geom.GetCell(cellid)
      p1   = cell.GetPointId(0)
      p2   = cell.GetPointId(1)
      pt   = geom.GetPoint(p1)
      pre  = pressure.GetValue(p1)
      Cp   = -(pre - pre_inf)/q_inf
      spt  = str(pt[0]) + " " + str(pt[1]) + " " + str(pt[2])
      spt  = spt + " " + str(Cp) + "\n"
      f.write(spt)
      cellid = NextCell(p2,geom)

   f.close()
   ncut = ncut + 1
