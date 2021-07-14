# -*- coding: utf-8 -*-
from xml.etree import ElementTree
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
from xml.dom import minidom

from urllib import request

from datetime import datetime, timedelta

from PyQt5.QtCore import QDate, QTime, QDateTime, Qt, QVariant
from PyQt5.QtGui import QColor
from qgis.core import QgsProject, QgsVectorLayer, QgsFeature, QgsField, QgsGeometry, QgsMessageLog, Qgis

import os

import re

import ogr

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from netCDF4 import Dataset, num2date

from collections import OrderedDict

import argparse

PLUGIN_NAME = 'FloodTool2'

class THREDDS_parser:

    URL_XML = {
    'artabro (MOHID)'     : 'http://193.144.35.143/thredds/catalog/MyCoast/MOHID/artabro/catalog.xml',
    'arousa (MOHID)'      : 'http://193.144.35.143/thredds/catalog/MyCoast/MOHID/arousa/catalog.xml',
    'vigo (MOHID)'        : 'http://193.144.35.143/thredds/catalog/MyCoast/MOHID/vigo/catalog.xml',
    'noia (MOHID)'        : 'http://193.144.35.143/thredds/catalog/MyCoast/MOHID/noia/catalog.xml',
    'iberia (ROMS)'       : 'http://193.144.35.143/thredds/catalog/MyCoast/ROMS/iberia/catalog.xml',
    'cantabrico (ROMS)'   : 'http://193.144.35.143/thredds/catalog/MyCoast/ROMS/cantabrico/catalog.xml',
    'wrf12km (WRF)'       : 'http://193.144.35.143/thredds/catalog/MyCoast/WRF/iberia/catalog.xml',
    'wrf04km (WRF)'       : 'http://193.144.35.143/thredds/catalog/MyCoast/WRF/galicia/catalog.xml',
    'galicia (SWAN)'      : 'http://193.144.35.143/thredds/catalog/MyCoast/SWAN/galicia/catalog.xml',
    'galicia (WW3)'       : 'http://193.144.35.143/thredds/catalog/MyCoast/WW3/galicia/catalog.xml',
    'iberia (WW3)'        : 'http://193.144.35.143/thredds/catalog/MyCoast/WW3/iberia/catalog.xml',
    }


    def __init__(self, model):

        self.URL = self.URL_XML[model]

    def parse_dates(self):

        request.urlopen(self.URL)

        contenido = ''.join([linea.decode("utf-8") for linea in request.urlopen(self.URL).readlines()])

        XML  = ElementTree.fromstring(contenido)

        filtered = XML.findall('{http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0}dataset/{http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0}dataset')

        dates = [datetime.strptime(re.findall(r'\d{10}',element.attrib['name'])[0],'%Y%m%d%H') for element in filtered if 'MyCOAST' in element.attrib['name']]

        return dates

class modelGrid:

    '''An abstraction of different model grids that will be used'''

    url_templates = {
    'artabro (MOHID)'     : 'http://193.144.35.143/thredds/dodsC/MyCoast/MOHID/artabro/MyCOAST_V1_MeteoGalicia_MOHID_artabro_01hr_%Y%m%d00_PR.ncml',
    'arousa (MOHID)'      : 'http://193.144.35.143/thredds/dodsC/MyCoast/MOHID/arousa/MyCOAST_V1_MeteoGalicia_MOHID_arousa_01hr_%Y%m%d00_PR.ncml',
    'vigo (MOHID)'        : 'http://193.144.35.143/thredds/dodsC/MyCoast/MOHID/vigo/MyCOAST_V1_MeteoGalicia_MOHID_vigo_01hr_%Y%m%d00_PR.ncml',
    'noia (MOHID)'        : 'http://193.144.35.143/thredds/dodsC/MyCoast/MOHID/noia/MyCOAST_V1_MeteoGalicia_MOHID_noia_01hr_%Y%m%d00_PR.ncml',
    'iberia (ROMS)'       : 'http://193.144.35.143/thredds/dodsC/MyCoast/ROMS/iberia/MyCOAST_V1_MeteoGalicia_ROMS_iberia_01hr_%Y%m%d00_PR.ncml',
    'cantabrico (ROMS)'   : 'http://193.144.35.143/thredds/dodsC/MyCoast/ROMS/cantabrico/MyCOAST_V1_AZTI_ROMS_cantabrico_01hr_%Y%m%d00_PR.ncml',
    'wrf12km (WRF)'       : 'http://193.144.35.143/thredds/dodsC/MyCoast/WRF/iberia/MyCOAST_V1_MeteoGalicia_WRF_iberia_01hr_%Y%m%d00_PR.ncml',
    'wrf04km (WRF)'       : 'http://193.144.35.143/thredds/dodsC/MyCoast/WRF/galicia/MyCOAST_V1_MeteoGalicia_WRF_galicia_01hr_%Y%m%d00_PR.ncml',
    'galicia (SWAN)'      : 'http://193.144.35.143/thredds/dodsC/MyCoast/SWAN/galicia/MyCOAST_V1_MeteoGalicia_SWAN_galicia_01hr_%Y%m%d00_PR.ncml',
    'galicia (WW3)'       : 'http://193.144.35.143/thredds/dodsC/MyCoast/WW3/galicia/MyCOAST_V1_MeteoGalicia_WW3_galicia_01hr_%Y%m%d00_ANPR.ncml',
    'iberia (WW3)'        : 'http://193.144.35.143/thredds/dodsC/MyCoast/WW3/iberia/MyCOAST_V1_MeteoGalicia_WW3_iberia_01hr_%Y%m%d00_ANPR.ncml',
    }

    unstructured_flags = {
    'artabro (MOHID)'     : False,
    'arousa (MOHID)'      : False,
    'vigo (MOHID)'        : False,
    'noia (MOHID)'        : False,
    'iberia (ROMS)'       : False,
    'cantabrico (ROMS)'   : False,
    'wrf12km (WRF)'       : False,
    'wrf04km (WRF)'       : False,
    'galicia (SWAN)'      : True,
    'galicia (WW3)'       : False,
    'iberia (WW3)'        : False,
    }

    tide_solution_flags = {
    'artabro (MOHID)'     : False,
    'arousa (MOHID)'      : False,
    'vigo (MOHID)'        : False,
    'noia (MOHID)'        : False,
    'iberia (ROMS)'       : True,
    'cantabrico (ROMS)'   : False,
    'wrf12km (WRF)'       : False,
    'wrf04km (WRF)'       : False,
    'galicia (SWAN)'      : False,
    'galicia (WW3)'       : False,
    'iberia (WW3)'        : False,
    }

    def __init__(self, model):

        self.gridName = model

        self.THREDDS_parser = THREDDS_parser(model)

        if model not in self.url_templates.keys():

            QgsMessageLog.logMessage('No template for %s' % model, PLUGIN_NAME, level=Qgis.Critical)
            exit()

        else:

            self.template      = self.url_templates[model]
            self.unstructured  = self.unstructured_flags[model]
            self.tide_solution = self.tide_solution_flags[model]

        # Last element of available dates in THREDDS server for grid:
        origen = Dataset(self.THREDDS_parser.parse_dates()[-1].strftime(self.template) )

        self.variables = origen.variables.keys()

        # This translation is based on standard_name attribute (based on THREDDS data standardisation):
        standard_names_to_var = {}

        for key in origen.variables.keys():

            try:
                standard_names_to_var[origen.variables[key].standard_name] = key

            except:
                pass

        self.standard_names_to_var = standard_names_to_var

        # Search by standard_name attribute:
        self.lon   = origen.variables[standard_names_to_var['longitude']][:].astype('double')
        self.lat   = origen.variables[standard_names_to_var['latitude' ]][:].astype('double')

        # Time span of grid file predictions:
        time          = origen.variables[standard_names_to_var['time']]
        time          = num2date(time[:],time.units)
        self.timespan = time[-1]-time[0]

        if not self.unstructured:
            if len(self.lon.shape) == 1:
                self.lon, self.lat = np.meshgrid(self.lon, self.lat)

        self.Xmin,self.Ymin = self.lon.min(), self.lat.min()
        self.Xmax,self.Ymax = self.lon.max(), self.lat.max()

    def get_vectorLayer(self):

        '''Get model grid contour as a layer for QGIS'''

        vectorlayer = QgsVectorLayer("Linestring?crs=EPSG:4326", "Bounding box", "memory")

        segmento = ogr.Geometry(ogr.wkbLineString)

        for X,Y in zip(self.lon[0,:], self.lat[0,:]):
            segmento.AddPoint(X,Y)

        for X,Y in zip(self.lon[:,-1], self.lat[:,-1]):
            segmento.AddPoint(X,Y)

        for X,Y in zip(self.lon[-1,::-1], self.lat[-1,::-1]):
            segmento.AddPoint(X,Y)

        for X,Y in zip(self.lon[::-1,0], self.lat[::-1,0]):
            segmento.AddPoint(X,Y)

        geom = QgsGeometry.fromWkt(segmento.ExportToWkt())

        feature = QgsFeature()
        feature.setGeometry(geom)

        pr = vectorlayer.dataProvider()
        pr.addAttributes([QgsField("id"    , QVariant.Int)])

        vectorlayer.updateFields()

        feature.setAttributes([int(0)])
        pr.addFeature(feature)

        vectorlayer.renderer().symbol().setWidth(0.7)
        vectorlayer.renderer().symbol().setColor(QColor.fromRgb(0,137,0))

        proyecto = QgsProject.instance()
        proyecto.addMapLayer(vectorlayer)

    def get_boundingBox(self):

        ''' Get aproximate bounding box'''

        return self.Xmin, self.Ymin, self.Xmax, self.Ymax

    def get_dates(self):

        return list(self.THREDDS_parser.parse_dates())


class tidalSolution:

    def __init__(self):

        self.DATA_PATH = os.path.join(os.path.dirname(__file__), 'data')

        QgsMessageLog.logMessage('DATA_PATH: %s' % self.DATA_PATH, PLUGIN_NAME, level=Qgis.Info)

    def read_grid(self, nodos='z'):

        grid_file = os.path.join(self.DATA_PATH, 'RAIA_Iberia_002Grid_17Set2013.nc')

        datos     = Dataset(grid_file)

        if nodos=='z':

            # Dimensiones de las matrices:
            self.ni = len( datos.dimensions['xi_rho'] )
            self.nj = len( datos.dimensions['eta_rho'] )

            self.lon = datos.variables['lon_rho'][:]
            self.lat = datos.variables['lat_rho'][:]

        elif nodos=='u':

            # Dimensiones de las matrices:
            self.ni = len( datos.dimensions['xi_u'] )
            self.nj = len( datos.dimensions['eta_u'] )

            self.lon = datos.variables['lon_u'][:]
            self.lat = datos.variables['lat_u'][:]

        else:

            # Dimensiones de las matrices:
            self.ni = len( datos.dimensions['xi_v'] )
            self.nj = len( datos.dimensions['eta_v'] )

            self.lon = datos.variables['lon_v'][:]
            self.lat = datos.variables['lat_v'][:]

    def least_squares(self, nodos='z'):

        data_file = os.path.join(self.DATA_PATH, 'marea_operativo_365dias.nc')

        # Least square from a detiding model run:
        datos         = Dataset(data_file)

       # Número de constituyentes:
        NTC       = len( datos.dimensions['tide_period'] )

        # Número de componentes armónicas:
        harmonics = len( datos.dimensions['harmonics'] )

        # Dimensionamos:
        C  = np.empty((harmonics,harmonics))*np.nan
        Ak = np.empty((harmonics,self.nj,self.ni))

        # Cargamos los elementos de la matríz C:
        CosW     = datos.variables['CosW'][:]
        SinW     = datos.variables['SinW'][:]
        CosWCosW = datos.variables['CosWCosW'][:]
        SinWSinW = datos.variables['SinWSinW'][:]
        SinWCosW = datos.variables['SinWCosW'][:]

        # Número de datos acumulados:
        Hcount = datos.variables['Hcount'][0]

        # Este factor no vale para nada, solo para escalar los cálculos, o eso creo:
        fac1   = 1.0/Hcount
        QgsMessageLog.logMessage('Total de datos acumulados: %i' % Hcount, PLUGIN_NAME, level=Qgis.Info)

        # Datos acumulados:

        if nodos=='z':
            variable = datos.variables['zeta_tide'][:]
        elif nodos=='u':
            variable = datos.variables['ubar_tide'][:]
        else:
            variable = datos.variables['vbar_tide'][:]

        # Matriz de coeficientes. Depende del tiempo pero es única para cada I,J:      
        ## Coeficiente para Z0:
        C[0,0]=1.0  

        for nk in range(1,NTC+1):

            C[0,nk    ] = fac1*SinW[nk-1]
            C[0,nk+NTC] = fac1*CosW[nk-1]
            C[nk,0    ] = C[0,nk]             # symmetric
            C[nk+NTC,0] = C[0,nk+NTC]         # symmetric

            for mk in range(1,NTC+1):

              C[mk,nk]         = fac1*SinWSinW[mk-1,nk-1]
              C[mk,nk+NTC]     = fac1*SinWCosW[mk-1,nk-1]
              C[mk+NTC,nk]     = fac1*SinWCosW[nk-1,mk-1]
              C[mk+NTC,nk+NTC] = fac1*CosWCosW[mk-1,nk-1]

        # La matriz ya es simétrica. Para que este código?:
        for nk in range(NTC):
            for mk in range(NTC):
                C[nk,mk] = C[mk,nk]


        # Invertimos la matriz para obtener la solución del sistema en las variables Z0, Ak y Bk:
        Y  = np.linalg.inv(C)
        Y *= fac1

        for i in range(self.ni):
            for j in range(self.nj):
                Ak[:,j,i] = np.dot(Y,variable[:,j,i])

        self.Z0         = Ak[0,:]
        self.Amplitudes = np.sqrt(Ak[1:9,:]*Ak[1:9,:] + Ak[9:18,:]*Ak[9:18,:])
        self.Fases      = np.arctan2(Ak[1:9,:],Ak[9:18,:])


    def tideSynthesis(self, time_index, nodes):

        # Frecuencias fundamentales:
        f = np.array( [14.49205211, 0.54901653, 0.04106864, 0.00464184, -0.00220641] )/360

        ## #################################################################################
        # Armonicos del ROMS:
        constituyentes  = ['k2','s2','m2','n2','k1','p1','o1','q1']

        # Números de Doodson:
        doodson = np.array( [ 
                              [2,  2,  0,  0,  0], # K2
                              [2,  2, -2,  0,  0], # S2
                              [2,  0,  0,  0,  0], # M2
                              [2, -1,  0,  1,  0], # N2
                              [1,  1,  0,  0,  0], # K1
                              [1,  1, -2,  0,  0], # P1
                              [1, -1,  0,  0,  0], # O1
                              [1, -2,  0,  1,  0]  # Q1
                            ] )

        frecuencia = np.dot(doodson,f)

        ## #################################################################################
        # Prueba de sintesis:

        # Calculamos las fases de cada armónico:
        t     = (time_index - datetime(2005,1,1,0,0,0)).astype('timedelta64[h]')
        x     = pd.DataFrame(data=np.outer(2*np.pi*t,frecuencia), columns=constituyentes)
        table = pd.DataFrame(index=time_index)

        for node in np.unique(nodes):

            # Translate from node index to coordinates:
            j, i = node//self.ni, node%self.ni

            # Sumamos todas las series:
            offset = 2.08
            tmp = self.Amplitudes[:,j,i]*np.cos(x - self.Fases[:,j,i])
            tmp = tmp.sum(axis=1) + offset
            
            table[node] = tmp.values

        return table
