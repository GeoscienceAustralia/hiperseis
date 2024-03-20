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
import click
import os
import sys
import math
import getopt, os, string, sys, math, time
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

    return np.array(xout)
#end

def processFile(shpFileName, outputFileStem=None,
                min_lon=None, min_lat=None,
                max_lon=None, max_lat=None,
                salvus_domain=None):
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

    dom = None
    if(salvus_domain is not None):
        from seismic.ASDFdatabase.domain import Domain
        dom = Domain(salvus_domain)
    # end if

    #process points
    xyzPolygons=[]

    for p in polygons:
        xyzPolygon=[]
        for pt in p:
            lon = pt[0]
            lat = pt[1]

            if(min_lon!=None and min_lat!=None and max_lon!=None and max_lat!=None):
                if (lon > max_lon or lon < min_lon): continue
                if (lat > max_lat or lat < min_lat): continue
            elif(dom is not None):
                if(not dom.contains(lon, lat)): continue
            # end if
            r = 6371e3 # Earth radius in km
            t = np.radians(90-lat)
            p = np.radians(lon)

            xyz = rtp2xyz(r, t, p)

            xyzPolygon.append(xyz)
        #end for
        if(len(xyzPolygon)):
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

    if(outputFileStem): outFileName = outputFileStem+'.vtp'
    else: outFileName = os.path.splitext(shpFileName)[0]+'.vtp'

    writer.SetFileName(outFileName)
    writer.SetDataModeToAscii()
    o = masterVTKPolygons.GetOutput()

    writer.SetInputData(o)
    print ('Writing %s..' % (outFileName))
    writer.Write()
#end function

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('input', required=True,
                type=click.Path(exists=True))
@click.option('--output-file-stem', default=None, type=str,
              help="Stem of output file; the default, 'none', uses the input file name without the extension")
@click.option('--min-lat', default=None, type=float,
              help="Minimum latitude for clipping shape file")
@click.option('--max-lat', default=None, type=float,
              help="Minimum latitude for clipping shape file")
@click.option('--min-lon', default=None, type=float,
              help="Minimum longitude for clipping shape file")
@click.option('--max-lon', default=None, type=float,
              help="Maximum longitude for clipping shape file")
@click.option('--salvus-domain', default=None, type=click.Path(exists=True),
              help="Clip shapefile based on Salvus chunk domain in json format")
def process(input, output_file_stem, min_lat, max_lat, min_lon, max_lon, salvus_domain):
    """
    Script for converting a shapefile to VTK format

    INPUT: Path to shapefile

    Example Usage:
        python shapefile_to_vtk.py sample.shp
    Example Usage (clip geometry):
        python shapefile_to_vtk.py sample.shp --min-lon 100.5 --min-lat -53.5 --max-lon 189.5 --max-lat -0.5
    """
    clip = np.array([min_lon, min_lat, max_lon, max_lat])

    if(np.any(clip) and salvus_domain is not None):
        raise RuntimeError('Clipping is performed using either --min/max[lon/lat] or --salvus-domain, '
                           'but not both simultaneously. Aborting..')
    # end if

    if(np.sum(clip==None)>0 and np.sum(clip==None)!=4):
        raise RuntimeError('To clip geometry, all four parameters '
                           '(minlon, minlat, maxlon, maxlat) must be passed in')
    # end if

    processFile(input, output_file_stem, min_lon, min_lat, max_lon, max_lat, salvus_domain)
# end func

if __name__=="__main__":
    process()
#end if
