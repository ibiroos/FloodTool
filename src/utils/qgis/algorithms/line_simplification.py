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
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterFeatureSink,
                       QgsGeometry,
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

    INPUT = 'INPUT'
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
        return 'simplelinefilter'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr('Simple line filter')

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
        return self.tr("This algorithm uses a gaussian filter to simplify coast lines and contours")

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
                name='Submuestreo',
                description=self.tr('Subsampling resolution for the output line'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.5,
                optional=False
            )
        )        
        
        self.addParameter(
            QgsProcessingParameterNumber(
                name='Gaussiano',
                description=self.tr('Number of points for the gaussian filter (sigma=points*resolution)'),
                type=QgsProcessingParameterNumber.Integer,
                defaultValue=20,
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
        

        submuestreo = self.parameterAsDouble(
            parameters,
            'Submuestreo',
            context
        )

        feedback.pushInfo('Submuestreo: %f' % submuestreo)
        
        gaussiano = self.parameterAsInt(
            parameters,
            'Gaussiano',
            context
        )        
        
        feedback.pushInfo('Gaussiano: %i' % gaussiano)


        # If source was not found, throw an exception to indicate that the algorithm
        # encountered a fatal error. The exception text can be any string, but in this
        # case we use the pre-built invalidSourceError method to return a standard
        # helper text for when a source cannot be evaluated
        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT))

        campos = QgsFields()
        campos.append( QgsField('id',QVariant.Int) )

        (sink, dest_id) = self.parameterAsSink(
            parameters,
            self.OUTPUT,
            context,
            campos, #QgsFields(), # source.fields(),
            source.wkbType(),
            source.sourceCrs()
        )

        # Send some information to the user
        feedback.pushInfo('CRS is {}'.format(source.sourceCrs().authid()))

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

        # Distancias:
        d = np.zeros(n)

        # Calculamos la variable distancia que nos permitirá hacer el filtrado opcional de las lineas:
        dx = x[1:n] - x[0:n-1]
        dy = y[1:n] - y[0:n-1]

        d[1:] = (np.sqrt(dx**2 + dy**2)).cumsum()

        # El parámetro lo hacemos mas denso para representar correctamente la linea:        
        longitud = d[::-1][0]

        l     = 0.5
        sigma = 20 # Filtro a l*sigma metros
        
        l     = submuestreo
        sigma = gaussiano   # Filtro a l*sigma metros        

        divisiones = round(longitud/l)

        feedback.pushInfo( 'l old: %f , l new: %f' % (l, longitud/divisiones) )

        d1 = np.linspace(0,longitud,divisiones+1) # Numero de puntos y no de intervalos

        # El nuevo tamaño es:
        n1 = len(d1)

        # Interpolamos:
        x1 = np.interp(d1,d,x)
        y1 = np.interp(d1,d,y)

        # Filtrado gaussiano (Conservamos el primer/ultimo punto):
        N = len(x1) + 2

        # Distancias:
        X = np.zeros(N)
        Y = np.zeros(N)

        D = np.zeros(N)

        X[0] = x1[0]
        Y[0] = y1[0]

        X[-1] = x1[-1]
        Y[-1] = y1[-1]

        # Filtramos la linea para que las transiciones sean mas suaves:
        X[1:-1] = gaussian_filter1d(x1, sigma)
        Y[1:-1] = gaussian_filter1d(y1, sigma)

        ## Necesitamos volver a interpolar:
        # Recalculamos la distancia
        dX = X[1:N] - X[0:N-1]
        dY = Y[1:N] - Y[0:N-1]

        D[1:] = (np.sqrt(dX**2 + dY**2)).cumsum()

        # Ajustamos de nuevo d1 para reflejar los cambios:
        longitud = D[::-1][0]

        divisiones = round(longitud/l)

        feedback.pushInfo('l old: %f , l new: %f' % (l, longitud/divisiones) )
                
        d1 = np.linspace(0,longitud,divisiones+1)

        x = np.interp(d1,D,X)[:]
        y = np.interp(d1,D,Y)[:]
        d = d1[:]
        
        if False:
            x = np.interp(d1,D,X)[::-1]
            y = np.interp(d1,D,Y)[::-1]
            d = d1[::-1]

        indices = range(1) # Una única linea
        
        # Compute the number of steps to display within the progress bar and
        # get features from source
        total = 100.0 / len(indices) if len(indices) else 0        
        
        #for i, (inicio, fin) in enumerate(indices):
        for i in indices:
            
            # Stop the algorithm if cancel button has been clicked
            if feedback.isCanceled():
                break            
                
            segmento = ogr.Geometry(ogr.wkbLineString)
            
            for X,Y in zip(x,y):
                
                segmento.AddPoint(X,Y)
                
            geom = QgsGeometry.fromWkt(segmento.ExportToWkt())
            
            feature = QgsFeature()
            feature.setGeometry(geom)

            feature.setAttributes([int(i)])
            
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
