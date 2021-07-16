# -*- coding: utf-8 -*-
"""
/***************************************************************************
 FloodTool2DockWidget
                                 A QGIS plugin
 MyCoast FloodTool
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                             -------------------
        begin                : 2019-02-27
        git sha              : $Format:%H$
        copyright            : (C) 2019 by Meteogalicia
        email                : pablo.enrique.carracedo.garcia@xunta.es
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

import os

from PyQt5 import QtGui, QtWidgets, uic
from PyQt5.QtCore import pyqtSignal

from PyQt5.QtCore import QDate, QTime, QDateTime, Qt, QVariant

from qgis.core import QgsProject, QgsVectorLayer, QgsFeature, QgsField, QgsMessageLog, Qgis

from datetime import datetime
import re
from netCDF4 import Dataset, num2date
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt

from .utils import THREDDS_parser, modelGrid, tidalSolution

'''
FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'flood_tool2_dockwidget_base.ui'))
'''

from .flood_tool2_dockwidget_base import Ui_FloodTool2DockWidgetBase

class FloodTool2DockWidget(QtWidgets.QDockWidget, Ui_FloodTool2DockWidgetBase):

    closingPlugin = pyqtSignal()

    def __init__(self, iface, parent=None):

        """Constructor."""

        super(FloodTool2DockWidget, self).__init__(parent)

        # Set up the user interface from Designer.
        # After setupUI you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://doc.qt.io/qt-5/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect

        self.setupUi(self)

        self.progressBar.setValue(0)

        # Setting up combo boxes and connecting signals:
        self.hydro_grid_comboBox.addItems(['artabro (MOHID)', 'arousa (MOHID)', 'vigo (MOHID)', 'noia (MOHID)', 'iberia (ROMS)', 'cantabrico (ROMS)'])
        self.hydro_grid_comboBox.currentIndexChanged.connect(self.enable_calendar_hydro)

        self.wave_grid_comboBox.addItems(['galicia (SWAN)', 'galicia (WW3)', 'iberia (WW3)'])
        self.wave_grid_comboBox.currentIndexChanged.connect(self.enable_calendar_wave)

        # Triggers initial evaluation of dates:
        self.enable_calendar_hydro()
        self.enable_calendar_wave()

        self.tideSolution_checkBox.stateChanged.connect(self.checkBoxState)

        self.run_button.clicked.connect(self.run)

        # The plugin acts over active layer:
        campos = [field.name() for field in iface.activeLayer().fields()]

        self.wave_field_comboBox.addItems(campos)
        self.hydro_field_comboBox.addItems(campos)

        # Access to iface from plugin:
        self.iface = iface

    def checkBoxState(self):

        self.hydro_variable_comboBox.setEnabled( not self.tideSolution_checkBox.isChecked() )

    def closeEvent(self, event):

        self.closingPlugin.emit()
        event.accept()

    def run(self):

        # QRadioButton state checker:
        _, self.time_control = [button.isChecked() for button in self.groupBox_time.findChildren(QtWidgets.QRadioButton)]
        _, self.color_codes  = [button.isChecked() for button in self.groupBox_codes.findChildren(QtWidgets.QRadioButton)]

        self.progressBar.setValue(0)

        text = ''

        text += 'Wave variable: %s\n' % self.wave_variable_comboBox.currentText()
        text += 'Wave layer field: %s\n' % self.wave_field_comboBox.currentText()
        text += 'Selected model date: %s\n' %  self.wave_calendarWidget.selectedDate().toString("yyyy-MM-dd")

        text += '----------------------------------------\n'
        
        text += 'Hydro variable: %s\n' % self.hydro_variable_comboBox.currentText()
        text += 'Hydro layer field: %s\n' % self.hydro_field_comboBox.currentText()
        text += 'Selected model date: %s\n' %  self.hydro_calendarWidget.selectedDate().toString("yyyy-MM-dd")

        text += '\nRetrieving active layer fields...\n'

        self.results_textBrowser.setText(text)

        features = self.iface.activeLayer().getFeatures()
        
        wave_id  = []
        hydro_id = []

        # Accedemos a la linea para obtener los vertices:
        for current, feature in enumerate(features):
          
            wave_id.append( feature[self.wave_field_comboBox.currentText()] )
            hydro_id.append( feature[self.hydro_field_comboBox.currentText()] )

        wave_id  = np.array(wave_id)
        hydro_id = np.array(hydro_id)

        self.progressBar.setValue(20)

        # Processing of wave model data:
        self.results_textBrowser.append('Retrieving THREDDS wave data...\n')

        wave_date  = self.wave_calendarWidget.selectedDate().toPyDate()
        wave_url   = wave_date.strftime(self.waveGrid.template)

        QgsMessageLog.logMessage('Wave URL: %s' % wave_url, 'DEBUG', level=Qgis.Info)

        datos = Dataset(wave_url)

        var_name   = self.waveGrid.standard_names_to_var[self.wave_variable_comboBox.currentText()]
        wave_data  = datos.variables[var_name][:]
        wave_dates = num2date(datos.variables['time'][:], datos.variables['time'].units)

        datos.close()

        # Reshaping of matrices in order to translate i,j to nodes ids in active layer:
        if len(wave_data.shape)==3:
            (nt, j, i) = wave_data.shape
            wave_data = wave_data.reshape(nt,j*i)

        self.progressBar.setValue(40)

        # Processing of hydrodynamic model data:
        self.results_textBrowser.append('Retrieving THREDDS hydrodynamic data...\n')

        hydro_date = self.hydro_calendarWidget.selectedDate().toPyDate()
        hydro_url  = hydro_date.strftime(self.hydroGrid.template)

        if self.tideSolution_checkBox.isChecked():

            romsTide = tidalSolution()
            romsTide.read_grid()
            romsTide.least_squares()

            hydro_dates = wave_dates

        else:

            datos = Dataset(hydro_url)
            var_name    = self.hydroGrid.standard_names_to_var[self.hydro_variable_comboBox.currentText()]
            hydro_data  = datos.variables[var_name][:]
            hydro_dates = num2date(datos.variables['time'][:], datos.variables['time'].units)
            datos.close()

            # Reshaping of matrices in order to translate i,j to nodes ids in active layer:
            if len(hydro_data.shape)==3:
                (nt, j, i) = hydro_data.shape
                hydro_data = hydro_data.reshape(nt,j*i)

        self.progressBar.setValue(60)

        self.results_textBrowser.append('Performing date selection over data...\n')

        # Date selection: Only process data where dates overlap:
        start = np.max((hydro_dates.min(), wave_dates.min()))
        end   = np.min((hydro_dates.max(), wave_dates.max()))

        self.results_textBrowser.append('Start: %s, End: %s\n' % (start.strftime('%Y/%m/%d %H:%M'), end.strftime('%Y/%m/%d %H:%M')))

        condition = (wave_dates>=start) & (wave_dates<=end)
        wave_data = wave_data[condition][:, wave_id]

        condition = (hydro_dates>=start) & (hydro_dates<=end)

        '''
        plt.plot(hydro_dates[condition], hydro_data[condition][:, 4239]+2.08,'r--')
        tideSerie = romsTide.tideSynthesis(pd.date_range(start,end,freq='1H'), [128325])
        tideSerie[128325].plot()
        plt.show()
        '''

        if self.tideSolution_checkBox.isChecked():
            hydro_data = romsTide.tideSynthesis(pd.date_range(start,end,freq='1H'), hydro_id)
            hydro_data = hydro_data[hydro_id].values # Needs a reshape cause tideSynthesis only calculate in unique id
        else:
            hydro_data = hydro_data[condition][:, hydro_id]

        QgsMessageLog.logMessage('hydro_data shape: (%i,%i)' % hydro_data.shape, 'DEBUG', level=Qgis.Info)

        self.progressBar.setValue(80)

        self.results_textBrowser.append('Running FL calculations...\n')

        # Very simple parametrization of flood level based on Vousdoukas et al. 2017:
        if self.time_control:
            FL          = hydro_data + 0.2*wave_data
            fechas      = pd.date_range(start,end,freq='1H')
            self.fechas = [fecha.strftime('%Y/%m/%d %H:%M') for fecha in fechas]

        else:
            FL = np.max(hydro_data + 0.2*wave_data, axis=0)

        # Processing output layer:
        features = self.iface.activeLayer().getFeatures()
        
        new_features = []

        # Accedemos a la linea para obtener los vertices:
        for current, feature in enumerate(features):

            field_names = [field.name() for field in feature.fields()]
          

            if self.time_control:
                for i,fecha in enumerate(self.fechas):
                    new_feature = QgsFeature(feature)
                    new_feature.setAttributes([current, float(FL[i,current]), fecha])
                    new_features.append(new_feature)
            else:
                new_feature = QgsFeature(feature)
                new_feature.setAttributes([current, float(FL[current])])
                new_features.append(new_feature)

        # Needs code for CRS:
        Crc_source_id = int(self.iface.activeLayer().sourceCrs().authid().split(':')[-1])
        vectorlayer = QgsVectorLayer("MultiPolygon?crs=epsg:%i" % Crc_source_id, "Output", "memory")

        pr = vectorlayer.dataProvider()

        if self.time_control:
            atributos = [QgsField("Id", QVariant.Int), QgsField("FL", QVariant.Double), QgsField("time"  , QVariant.String)]
        else:
            atributos = [QgsField("Id", QVariant.Int), QgsField("FL", QVariant.Double)]

        pr.addAttributes(atributos)

        vectorlayer.updateFields()
        pr.addFeatures(new_features)

        # Add layer to project:
        proyecto = QgsProject.instance()
        proyecto.addMapLayer(vectorlayer)

        self.progressBar.setValue(100)

    def enable_calendar_hydro(self, **kwargs):

        calendarWidget    = self.hydro_calendarWidget
        variable_comboBox = self.hydro_variable_comboBox

        self.hydroGrid = modelGrid(self.hydro_grid_comboBox.currentText())

        self.tideSolution_checkBox.setEnabled(self.hydroGrid.tide_solution)
    
        fechas = [fecha.strftime('%Y%m%d') for fecha in self.hydroGrid.THREDDS_parser.parse_dates()]

        inicio = datetime.strptime(fechas[ 0],'%Y%m%d')
        fin    = datetime.strptime(fechas[-1],'%Y%m%d')

        calendarWidget.setDateRange(QDate(inicio), QDate(fin))

        variable_comboBox.clear()
        #variable_comboBox.addItems(self.hydroGrid.variables)
        variable_comboBox.addItems(self.hydroGrid.standard_names_to_var.keys())

    def enable_calendar_wave(self, **kwargs):

        calendarWidget    = self.wave_calendarWidget
        variable_comboBox = self.wave_variable_comboBox

        self.waveGrid = modelGrid(self.wave_grid_comboBox.currentText())

        fechas = [fecha.strftime('%Y%m%d') for fecha in self.waveGrid.THREDDS_parser.parse_dates()]

        inicio = datetime.strptime(fechas[ 0],'%Y%m%d')
        fin    = datetime.strptime(fechas[-1],'%Y%m%d')

        calendarWidget.setDateRange(QDate(inicio), QDate(fin))

        variable_comboBox.clear()
        #variable_comboBox.addItems(self.waveGrid.variables)
        variable_comboBox.addItems(self.waveGrid.standard_names_to_var.keys())
