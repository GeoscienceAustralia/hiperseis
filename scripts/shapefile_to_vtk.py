#!/bin/env python

"""
Description:
    Script for converting a shapefile to vtk format
References:

CreationDate:   22/06/18
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     22/06/18   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import shapefile
import vtk
import os
import sys
import math
import getopt, os, string, sys, math, time, commands
import datetime
import pprint
import numpy as np
from scipy import spatial
import numpy

msg = '''usage: shapefileToVtk.py <shapefile.shp>'''

def rtp2xyz(r, theta, phi):
    rst = r * math.sin(theta)
    xout = [0]*3
    xout[0] = rst * math.cos(phi)	    # x 
    xout[1] = rst * math.sin(phi) 	    # y 
    xout[2] = r * math.cos(theta)       # z

    return xout
#end

def processFile(shpFileName):
    # Open the extracted islands
    r = shapefile.Reader(shpFileName)
    # Setup the world to pixels conversion
    xdist = r.bbox[2] - r.bbox[0]
    ydist = r.bbox[3] - r.bbox[1]
    iwidth = 800
    iheight = 600
    xratio = iwidth/xdist
    yratio = iheight/ydist
    polygons = []
    
    # Find field-idx of 'PLATE_ID'
    fields = r.fields

    # Loop through all shapes
    for shape, shaperec in zip(r.shapes(), r.shapeRecords()):
        # Loop through all parts to catch
        # polygon holes!
        for i in range(len(shape.parts)):
            pixels=[]
            pt = None

            if i < len(shape.parts)-1:
                pt = shape.points[shape.parts[i]:shape.parts[i+1]]
            else:
                pt = shape.points[shape.parts[i]:]
            #end if
            for x,y in pt:
                pixels.append([x,y])
            #end for
            polygons.append(pixels)
        #end for
    #end for
    
    #process points
    xyzPolygons=[]

    for p in polygons:
        xyzPolygon=[]
        for pt in p:
            r = 6371;
            t = (90-pt[1])/180*math.pi
            p = (360+pt[0])%360/180*math.pi

            xyz=rtp2xyz(r,t,p)
            xyzPolygon.append(xyz)
        #end for
        xyzPolygons.append(xyzPolygon)
    #end for

    #write VTK
    count = 0
    masterVTKPolygons = vtk.vtkAppendPolyData()
    for xyzp in xyzPolygons:
        dArray = vtk.vtkFloatArray()
        dArray.SetName('d')
        dArray.SetNumberOfValues(len(xyzp))

        points = vtk.vtkPoints()
        xyzptPrev = None
        for i, xyzpt in enumerate(xyzp):
            points.InsertNextPoint(xyzpt)

            if(xyzptPrev is None):
                dArray.SetValue(i, 0)
            else:
                dArray.SetValue(i, np.linalg.norm(np.array(xyzptPrev) - np.array(xyzpt)))
            #end if

            xyzptPrev = xyzpt
        #end for
        npoints = points.GetNumberOfPoints()
        lines = vtk.vtkCellArray()
        lines.InsertNextCell(npoints)
        
        points = vtk.vtkPoints()
        xyzptPrev = None
        dlist = []
        for i, xyzpt in enumerate(xyzp):
            points.InsertNextPoint(xyzpt)

            if(xyzptPrev is None):
                dlist.append(0)
            else:
                dlist.append(np.linalg.norm(np.array(xyzptPrev) - np.array(xyzpt)))
            #end if

            xyzptPrev = xyzpt
        #end for
        
        dlist = np.array(dlist)
        for i, xyzpt in enumerate(xyzp):
            dArray.SetValue(i, dlist[i])
        #end for
        
        npoints = points.GetNumberOfPoints()
        for i in range(npoints):
            lines.InsertCellPoint(i)
        #end for
        
        polyline = vtk.vtkPolyData()
        polyline.SetPoints(points)
        polyline.SetLines(lines)
        polyline.GetPointData().SetScalars(dArray)
        
        masterVTKPolygons.AddInputData(polyline)
        count += 1
    #end for

    masterVTKPolygons.Update()

    writer = vtk.vtkXMLPolyDataWriter()
    outFileName = os.path.splitext(shpFileName)[0]+'.vtp';
    writer.SetFileName(outFileName)
    writer.SetDataModeToAscii()
    o = masterVTKPolygons.GetOutput()
    #print o
    writer.SetInputData(o)
    print ('Writing %s..') % (outFileName)
    writer.Write()
#end function

if __name__=="__main__":
    if (len(sys.argv) != 2):
        print msg
        exit(0);
    #end if

    shpFileName = sys.argv[1];

    processFile(shpFileName)
#end if
