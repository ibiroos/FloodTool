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
                       QgsProcessingParameterFeatureSink,
                       QgsGeometry,
                       QgsPointXY,
                       QgsFields,
                       QgsField,
                       QgsFeature)
import processing

import numpy as np
from scipy.spatial import cKDTree, KDTree
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
        return 'attachlayer'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr('Attach layer fields by nearest neighbors')

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
        return self.tr("This algorithm searchs for nearest neighbors in other layer")

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
                self.tr('Initial layer'),
                [QgsProcessing.TypeVectorAnyGeometry]
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSource(
                'Attach',
                self.tr('Layer to link'),
                [QgsProcessing.TypeVectorAnyGeometry]
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

        attach = self.parameterAsSource(
            parameters,
            'Attach',
            context
        )

        # If source was not found, throw an exception to indicate that the algorithm
        # encountered a fatal error. The exception text can be any string, but in this
        # case we use the pre-built invalidSourceError method to return a standard
        # helper text for when a source cannot be evaluated
        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT))

        # Campos a incluir en el resultado:  
        campos = QgsFields()
        
        for field in source.fields():
            campos.append( field )
            
        campos.append( QgsField('%s_distance' % attach.sourceName(), QVariant.Double) )
        campos.append( QgsField('%s_neighbor' % attach.sourceName(), QVariant.Int) )

        (sink, dest_id) = self.parameterAsSink(
            parameters,
            self.OUTPUT,
            context,
            campos, # source.fields(),
            source.wkbType(),
            source.sourceCrs()
        )

        # Send some information to the user
        feedback.pushInfo('CRS is {}'.format(source.sourceCrs().authid()))

        Crc_source_id = int(source.sourceCrs().authid().split(':')[-1])
        Crc_attach_id = int(attach.sourceCrs().authid().split(':')[-1])

        if Crc_source_id!=Crc_attach_id:
            
            feedback.pushInfo('Ojo! Las proyecciones no son iguales. Source: %i, Attach_ %i' % (Crc_source_id, Crc_attach_id))
            feedback.pushInfo('Voy a tratar de hacer la proyección de Attach %i a Source: %i' % (Crc_attach_id, Crc_source_id))
            
            proyector = QgsCoordinateTransform(attach.sourceCrs(), source.sourceCrs(), Crc_attach_id, Crc_source_id)

        # If sink was not created, throw an exception to indicate that the algorithm
        # encountered a fatal error. The exception text can be any string, but in this
        # case we use the pre-built invalidSinkError method to return a standard
        # helper text for when a sink cannot be evaluated
        if sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, self.OUTPUT))
            
        features_source = source.getFeatures()
        
        source_x  = []
        source_y  = []
        source_id = []

        # Accedemos a la linea para obtener los vertices:
        for current, feature in enumerate(features_source):

            centroide_x = feature.geometry().centroid().asPoint().x()
            centroide_y = feature.geometry().centroid().asPoint().y()
            source_x.append( centroide_x )
            source_y.append( centroide_y )
                    
            source_id.append( feature['id'] )

        source_x = np.array(source_x)
        source_y = np.array(source_y)
        source_id = np.array(source_id)

        feedback.pushInfo('Got coordinates of initial layer')

        features_attach = attach.getFeatures()
        
        attach_x  = []
        attach_y  = []
        attach_id = []

        # Accedemos a la linea para obtener los vertices:
        for current, feature in enumerate(features_attach):

            point_x, point_y = feature.geometry().asPoint().x(), feature.geometry().asPoint().y()
            
            # Chech: Only reproject if needed:            
            if Crc_source_id!=Crc_attach_id:
                
                proyectado = proyector.transform(point_x, point_y)
                point_x, point_y = proyectado.x(), proyectado.y()

            attach_x.append( point_x )
            attach_y.append( point_y )
            
            attach_id.append( feature['id'] ) # Id from table field not from .id() method

        attach_x  = np.array(attach_x)
        attach_y  = np.array(attach_y)
        attach_id = np.array(attach_id)

        feedback.pushInfo('Got coordinates to append')

        # Busqueda cKDTree:
        nearest = cKDTree(np.column_stack([attach_x, attach_y]))
        distancias, vecinos = nearest.query(np.column_stack([source_x, source_y]))
        
        # Compute the number of steps to display within the progress bar and
        # get features from source
        total = 100.0 / len(source_x) if len(source_x) else 0
        
        # Si se usa el iterador, este queda vacío :-(
        features_source = source.getFeatures()        
        
        for i, feature in enumerate(features_source):
            
            # Stop the algorithm if cancel button has been clicked
            if feedback.isCanceled():
                break
            
            atributos = feature.attributes()
            
            # feedback.pushInfo("i=%i --> id=%i" % (i, feature['id']))

            # Añadimos los nuevos campos:
            atributos.append(float(distancias[i]))
            atributos.append(int(attach_id[vecinos[i]]))

            # Actualizamos la lista de campos de la feature:
            feature.setFields(campos)

            feature.setAttributes(atributos)

            sink.addFeature(feature, QgsFeatureSink.FastInsert)
            
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
                'INPUT': dest2_id,
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

