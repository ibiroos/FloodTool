# -*- coding: utf-8 -*-

"""
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

from PyQt5.QtCore import QCoreApplication, QVariant

from qgis.core import (QgsProcessing,
                       QgsFeatureSink,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsCoordinateReferenceSystem,
                       QgsCoordinateTransform, 
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterFeatureSink,
                       QgsGeometry,
                       QgsPointXY,
                       QgsFields,
                       QgsField,
                       QgsFeature)
import processing

import numpy as np
from scipy.ndimage import gaussian_filter1d
from osgeo import ogr


class ExampleProcessingAlgorithm(QgsProcessingAlgorithm):
    """
    This is an example algorithm that takes a vector layer and
    creates a new identical one.

    It is meant to be used as an example of how to create your own
    algorithms and explain methods and variables used to do it. An
    algorithm like this will be available in all elements, and there
    is not need for additional work.

    All Processing algorithms should extend the QgsProcessingAlgorithm
    class.
    """

    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    INPUT  = 'INPUT'
    OUTPUT = 'OUTPUT'

    def tr(self, string):
        """
        Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return ExampleProcessingAlgorithm()

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'bufferbylines'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr('Buffer by perpendicular lines')

    def group(self):
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return self.tr('MyCoast')

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'mycoast'

    def shortHelpString(self):
        """
        Returns a localised short helper string for the algorithm. This string
        should provide a basic description about what the algorithm does and the
        parameters and outputs associated with it..
        """
        return self.tr("This algorithm generates polygon buffer from simplified line")

    def initAlgorithm(self, config=None):
        """
        Here we define the inputs and output of the algorithm, along
        with some other properties.
        """

        # We add the input vector features source. It can have any kind of
        # geometry.
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT,
                self.tr('Input layer'),
                [QgsProcessing.TypeVectorAnyGeometry]
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                name='Radius',
                description=self.tr('Buffer radius (m)'),
                type=QgsProcessingParameterNumber.Integer,
                defaultValue=20,
                optional=False
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                name='Length',
                description=self.tr('Length of each polygon of the buffer'),
                type=QgsProcessingParameterNumber.Integer,
                defaultValue=100,
                optional=False
            )
        )

        # We add a feature sink in which to store our processed features (this
        # usually takes the form of a newly created vector layer when the
        # algorithm is run in QGIS).
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT,
                self.tr('Output layer')
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        """
        Here is where the processing itself takes place.
        """

        # Retrieve the feature source and sink. The 'dest_id' variable is used
        # to uniquely identify the feature sink, and must be included in the
        # dictionary returned by the processAlgorithm function.
        source = self.parameterAsSource(
            parameters,
            self.INPUT,
            context
        )


        radio = self.parameterAsDouble(
            parameters,
            'Radius',
            context
        )

        feedback.pushInfo('Radius: %f' % radio)
        
        longitud = self.parameterAsInt(
            parameters,
            'Length',
            context
        )        
        
        feedback.pushInfo('Length: %i' % longitud)

        # If source was not found, throw an exception to indicate that the algorithm
        # encountered a fatal error. The exception text can be any string, but in this
        # case we use the pre-built invalidSourceError method to return a standard
        # helper text for when a source cannot be evaluated
        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT))

        # Fields to add to the resulting layer:
        campos = QgsFields()
        campos.append( QgsField('id',QVariant.Int) )
        campos.append( QgsField('X_centroid', QVariant.Double) )
        campos.append( QgsField('Y_centroid', QVariant.Double) )
        campos.append( QgsField('Lon_centroid', QVariant.Double) )
        campos.append( QgsField('Lat_centroid', QVariant.Double) )

        (sink, dest_id) = self.parameterAsSink(
            parameters,
            self.OUTPUT,
            context,
            campos, # source.fields(),
            3, # = QGis.WKBPolygon (estaba a source.wkbType())
            source.sourceCrs()
        )

        # Send some information to the user
        crs_id = int(source.sourceCrs().authid().split(':')[1])
        feedback.pushInfo('CRS is {}'.format(crs_id))

        #proyector = QgsCoordinateTransform(QgsCoordinateReferenceSystem(23029), QgsCoordinateReferenceSystem(4326), 23029, 4326)
        proyector = QgsCoordinateTransform(QgsCoordinateReferenceSystem(crs_id), QgsCoordinateReferenceSystem(4326), crs_id, 4326)

        # If sink was not created, throw an exception to indicate that the algorithm
        # encountered a fatal error. The exception text can be any string, but in this
        # case we use the pre-built invalidSinkError method to return a standard
        # helper text for when a sink cannot be evaluated
        if sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, self.OUTPUT))

        features = source.getFeatures()
        
        x = []
        y = []

        # Accedemos a la linea para obtener los vertices:
        for current, feature in enumerate(features):
            for punto in feature.geometry().vertices():
                x.append( punto.x() )
                y.append( punto.y() )

        x = np.array(x)
        y = np.array(y)
        
        feedback.pushInfo('Got coordinates')
        
        # Número de vertices que contiene la linea:
        n = len(x)
        
        feedback.pushInfo('Number of line vertices: %i' % n)

        lineas = []

        R = radio

        for i in range(0,n-1):

            # Test perpendicular:
            # Dos puntos
            x0, y0 = x[i]  ,y[i]
            x1, y1 = x[i+1],y[i+1]

            # El punto medio del segmento:
            x2, y2 = (x0+x1)/2, (y0+y1)/2

            # feedback.pushInfo('Punto medio del segmento: (%f, %f)' % (x2,y2))

            # Pendiente de la recta perpendicular (-1/m de la original m):
            d   = np.sqrt((y1-y0)**2 + (x1-x0)**2)
            sin = (y1-y0)/d
            cos = (x1-x0)/d
            # m   = -(x1-x0)/(y1-y0) # tan

            # Intercept para que pase por el punto medio del segemento:
            # b = y2 - m*x2

            # X = np.linspace(-10,10,2000) + x2
            # Y = m*X + b

            # Coordenadas de los puntos extremos:
            # lineas.append( [(sin*R + x2, -cos*R + y2), (-sin*R + x2, cos*R + y2)] )
            # lineas.append( [(sin*R + x2, -cos*R + y2), (x2, y2)] )  # Interior
            lineas.append( [(x2, y2), (-sin*R + x2, cos*R + y2)] )  # Exterior

        feedback.pushInfo('Number of perpendicular lines: %i' % len(lineas))

        # Construimos los poligonos:
        nl        = longitud
        poligonos = [lineas[i*nl:(i+1)*nl+1] for i in range(len(lineas)//nl)]

        # Compute the number of steps to display within the progress bar and
        # get features from source
        total = 100.0 / len(poligonos) if len(poligonos) else 0
        
        for i, poligono in enumerate(poligonos):
            
            # Stop the algorithm if cancel button has been clicked
            if feedback.isCanceled():
                break

            puntos = []

            for P1, P2 in poligono:
                puntos.append(P1)

            for P1, P2 in poligono[::-1]:
                puntos.append(P2)

            ring = ogr.Geometry(ogr.wkbLinearRing)

            for Xp,Yp in puntos:
            
                ring.AddPoint(Xp,Yp)

            poligono = ogr.Geometry(ogr.wkbPolygon)
            poligono.AddGeometry(ring)
                
            geom = QgsGeometry.fromWkt(poligono.ExportToWkt())
            
            feature = QgsFeature()
            feature.setGeometry(geom)

            centroide_x = feature.geometry().centroid().asPoint().x()
            centroide_y = feature.geometry().centroid().asPoint().y()

            proyectado = proyector.transform(centroide_x, centroide_y)

            feature.setAttributes([int(i), centroide_x, centroide_y, proyectado.x(), proyectado.y() ])
            
            sink.addFeature(feature)
            
            feedback.setProgress(int(i * total))
            
        if False:
            for current, feature in enumerate(features):
                # Stop the algorithm if cancel button has been clicked
                if feedback.isCanceled():
                    break

                # Add a feature in the sink
                sink.addFeature(feature, QgsFeatureSink.FastInsert)

                # Update the progress bar
                feedback.setProgress(int(current * total))

        # To run another Processing algorithm as part of this algorithm, you can use
        # processing.run(...). Make sure you pass the current context and feedback
        # to processing.run to ensure that all temporary layer outputs are available
        # to the executed algorithm, and that the executed algorithm can send feedback
        # reports to the user (and correctly handle cancelation and progress reports!)
        if False:
            buffered_layer = processing.run("native:buffer", {
                'INPUT': dest_id,
                'DISTANCE': 1.5,
                'SEGMENTS': 5,
                'END_CAP_STYLE': 0,
                'JOIN_STYLE': 0,
                'MITER_LIMIT': 2,
                'DISSOLVE': False,
                'OUTPUT': 'memory:'
            }, context=context, feedback=feedback)['OUTPUT']

        # Return the results of the algorithm. In this case our only result is
        # the feature sink which contains the processed features, but some
        # algorithms may return multiple feature sinks, calculated numeric
        # statistics, etc. These should all be included in the returned
        # dictionary, with keys matching the feature corresponding parameter
        # or output names.
        return {self.OUTPUT: dest_id}
