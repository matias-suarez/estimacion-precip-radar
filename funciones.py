#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 1, 2021
@author: msuarez
"""

#-----------------------------------------------------------------
# Este modulo contiene las siguientes funciones:
#             1) get_nearest_gate_azimuth
#             2) radar_variable_lat_lon
#             3) radar_variable_window_lat_lon
#             4) radar_variable_window_lat_lon_list
#             5) create_df_from_radar_windows_vars
#             6) db_to_linear
#             7) linear_to_db
#             8) lat_lon_to_range_azimuth
#             9) calculate_R_from_Z_R
#             9) calculate_R_from_Z_R
#            10) calculate_R_from_Z_ZDR_R
#-----------------------------------------------------------------

from dateutil.relativedelta import relativedelta
import pandas as pd
import math
import pyart
from datetime import datetime
import numpy as np



def get_nearest_gate_azimuth(radar, longitud, latitud, verbose=False):
    """
    Funcion para calcular el angulo donde se ubica una estación contando desde el norte
    en sentido horario, la distancia de la estación al radar (en km), el nro de gate y el azimuth del
    del gate mas cercano. Y finalmente las coordenadas geograficas (lat, lon y altitud) del
    gate mas cercano.
    Para la ejecucion de esta funcion son necesarios los paquetes PyART y Math
    
    Parameters:
              radar (radar obj): objeto radar
              longitud (float): longitud del punto de interes (en grados decimales)
              latitud (float): latitud del punto de interes (en grados decimales)
              verbose (Bool): True para obtener los print de pantalla
    
    Returns:
              gate (int): es el numero del gate mas cercano al punto
              alfa (int): es el numero de rayo correspondiente al gate mas cercano
    """
    # Calculo coordenadas x,y de la estacion
    x,y = pyart.core.geographic_to_cartesian_aeqd(longitud,
                                                  latitud,
                                                  radar.longitude['data'][0],
                                                  radar.latitude['data'][0]) 

    

    # CALCULO ANGULO Y DISTANCIA AL PLUVIOMETRO
    r = math.sqrt(x*x+y*y)    
    if x > 0 and y > 0: # Si el punto x,y esta en el 1er cuadrante
        alfa = round(abs(math.degrees(math.atan(x/y)))) # Angulo redondeado al entero mas cercano
        theta = abs(math.degrees(math.atan(x/y)))       # Angulo exacto
    elif x > 0 and y < 0: # Si el punto x,y esta en el 2do cuadrante
        alfa = round(abs(math.degrees(math.atan(abs(y)/x)) + 90)) # Angulo redondeado al entero mas cercano
        theta = abs(math.degrees(math.atan(abs(y)/x)) + 90)       # Angulo exacto
    elif x < 0 and y > 0: # Si el punto x,y esta en el 4to cuadrante
        alfa = round(abs(math.degrees(math.atan(y/abs(x))) + 270)) # Angulo redondeado al entero mas cercano
        theta = abs(math.degrees(math.atan(y/abs(x))) + 270)       # Angulo exacto
    elif x < 0 and y < 0: # Si el punto x,y esta en el 3er cuadrante
        alfa = round(abs(math.degrees(math.atan(x/y)) + 180)) # Angulo redondeado al entero mas cercano
        theta = abs(math.degrees(math.atan(x/y)) + 180)       # Angulo exacto


    # DISTANCIA AL PLUVIOMETRO EN KM. Alfa contando desde el Norte geografico en sentido horario
    r_km = r/1000 # Paso r en metros a r_km en kilometros

    # CALCULO DE LA DISTANCIA DEL 1ER GATE AL RADAR
    # para el calculo de la distancia al radar
    first_gate_latitude = radar.gate_latitude['data'][0, 0] # Latitud del primer gate con azimuth 0°
    first_gate_longitude = radar.gate_longitude['data'][0, 0] # Longitud del primer gate con azimuth 0°
    first_x,first_y = pyart.core.geographic_to_cartesian_aeqd(first_gate_longitude,
                                                              first_gate_latitude,
                                                              radar.longitude['data'][0],
                                                              radar.latitude['data'][0]) # Coord. x y del 
                                                                            # primer gate con
                                                                            # azimuth 0
    r_first = math.sqrt(first_x*first_x+first_y*first_y) # Distancia al primer gate

    # CALCULO DE LA DISTANCIA DEL 2DO GATE AL RADAR
    # para el calculo de la distancia al radar
    second_gate_latitude = radar.gate_latitude['data'][0, 1] # Latitud del segundo gate con azimuth 0°
    second_gate_longitude = radar.gate_longitude['data'][0, 1] # Latitud del segundo gate con azimuth 0°
    second_x,second_y = pyart.core.geographic_to_cartesian_aeqd(second_gate_longitude,
                                                                second_gate_latitude,
                                                                radar.longitude['data'][0],
                                                                radar.latitude['data'][0]) # Coord. x y del
                                                                              # segundo gate
                                                                              # con azimuth 0
    r_second = math.sqrt(second_x*second_x+second_y*second_y) # Distancia al segundo gate

    # Distancia entre gates (dist 2do gate<->radar - dist 1er gate<->radar )
    dist_bet_gates = r_second - r_first 

    # (dist. pluv. - dist. 1st gate) / dist. entre gates = numero de gate
    gate = round((r-r_first)/dist_bet_gates) 

    # LATITUD, LONGITUD Y ALTITUD DEL GATE MAS CERCANO AL PLUVIOMETRO
    gate_latitude = radar.gate_latitude['data'][alfa, gate]
    gate_longitude = radar.gate_longitude['data'][alfa, gate]
    gate_altitude = radar.gate_altitude['data'][alfa, gate]

    if verbose:
        print('-----------------------------------------------')
        print('Angulo donde se ubica la estacion:',theta,'°','contando desde el norte geografico')
        print('El pluviometro esta a',r_km,'km del radar')
        print(alfa,'° es el angulo redondeado')
        print ('gate más cercano:',gate)
        print('Coordenadas del centro del gate', gate,'con azimuth',alfa,':')
        print('Latidud  ','  Longitud','Altitud')
        print(gate_latitude, gate_longitude, gate_altitude)
        print('-----------------------------------------------')
    
    return gate,alfa



def radar_variable_lat_lon(radarfilepath,
                           fields,
                           lat,
                           lon,
                           rhohv_field='RHOHV',
                           rhohv_threshold=0.8,
                           mask=False):
    """
    Función para obtener los valores de variables de radar sobre un punto con coordenadas
    lat y lon. La salida es una tupla: Fecha y lista con los valores extraidos.
    Parameters:
            radarfilepath (str): Path al archivo volumen de radar.
            fields (str o list): Nombres de los campos de radar a extraer el valor.
            lat (float): Latitud en grados decimales.
            lon (float): Longitud en grados decimales.
            rhohv_field (str): Nombre del campo RHOHV. Default 'RHOHV'.
            rhohv_threshold (float): Valor de RHOHV (0 a 1) para aplicar mascara. Default 0.8.
            mask (bool): True para aplicar mascara. Default False.

    Returns:
            fecha (DateTime object): Fecha y hora del volumen de radar
            result (list): Lista con los valores de interes. Estan ordenados de la misma
                           manera que se ingresan en la variable fields.
    """
    
    # Creamos el objeto "radar" con los campos definidos en include_fields
    radar = pyart.io.read(radarfilepath)

    # Extraemos la fecha
    fecha = radar.time['units'][14:]

    # La convertimos en un objeto datetime con el módulo datetime 
    fecha = datetime.strptime(fecha, '%Y-%m-%dT%H:%M:%SZ')
    
    field_dict = {}
    
    if mask:
        for field in fields:
            try:
                # Enmascaro los valores de los campos indicados
                field_dict[field] = np.ma.masked_where(radar.fields[rhohv_field]['data'] < rhohv_threshold,
                                                       radar.fields[field]['data'])
            except KeyError:
                print('Error. Campo no encontrado. Los campos disponibles en el volumen son',radar.fields.keys())

    else:
        for field in fields:
            try:
                field_dict[field] = radar.fields[field]['data']
            except KeyError:
                print('Error. Campo no encontrado. Los campos disponibles son:',radar.fields.keys())
            
            
    # Funcion para calcular el gate y azimuth sobre las coordenadas
    gate,alfa = get_nearest_gate_azimuth(radar, lon, lat)
    
    result = []
    
    for field in fields:
        # Extraigo los valores de los campos enmascarados
        value = field_dict[field][alfa,gate]
        result.append(value)
        
    del radar,field_dict
    return fecha,result



def radar_variable_window_lat_lon(radarfilepath,
                                  fields,
                                  lat,
                                  lon,
                                  rhohv_field='RHOHV',
                                  rhohv_threshold=0.8,
                                  mask=False):
    
    """
    Función para obtener los valores de variables de radar en una ventana de 3x3 celdas
    centradas sobre un punto con coordenadas lat y lon.
    La salida es una tupla: Fecha y lista de arrays 3x3 con los valores extraidos.
    La disposición de la malla de celdas es la siguiente:
    
            | alfa-1 | alfa  | alfa+1 |
     ----------------------------------
     gate+1 | [0,0]  | [0,1] |  [0,2] |
     gate   | [1,0]  | [1,1] |  [1,2] |
     gate-1 | [2,0]  | [2,1] |  [2,2] |
     ----------------------------------
    Parameters:
            radarfilepath (str): Path al archivo volumen de radar.
            fields (str o list): Nombres de los campos de radar a extraer el valor.
            lat (float): Latitud en grados decimales.
            lon (float): Longitud en grados decimales.
            rhohv_field (str): Nombre del campo RHOHV. Default 'RHOHV'.
            rhohv_threshold (float): Valor de RHOHV (0 a 1) para aplicar mascara. Default 0.8.
            mask (bool): True para aplicar mascara. Default False.

    Returns:
            fecha (DateTime object): Fecha y hora del volumen de radar
            result (list): Lista de los arrays con los valores de interes.
                           Estan ordenados de la misma manera que se ingresan en la variable fields.
    """
    
    # Creamos el objeto "radar" con los campos definidos en include_fields
    radar = pyart.io.read(radarfilepath)

    # Extraemos la fecha
    fecha = radar.time['units'][14:]

    # La convertimos en un objeto datetime con el módulo datetime 
    fecha = datetime.strptime(fecha, '%Y-%m-%dT%H:%M:%SZ')
    
    field_dict = {}
    
    if mask:
        for field in fields:
            try:
                # Enmascaro los valores de los campos indicados
                field_dict[field] = np.ma.masked_where(radar.fields[rhohv_field]['data'] < rhohv_threshold,
                                                       radar.fields[field]['data'])
            except KeyError:
                print('Error. Campo no encontrado. Los campos disponibles en el volumen son',radar.fields.keys())

    else:
        for field in fields:
            try:
                field_dict[field] = radar.fields[field]['data']
            except KeyError:
                print('Error. Campo no encontrado. Los campos disponibles son:',radar.fields.keys())
            
            
    # Funcion para calcular el gate y azimuth sobre las coordenadas
    gate,alfa = get_nearest_gate_azimuth(radar, lon, lat)
    
    result = []
    
    for field in fields:
        # np.full([shape], fillvalue)
        # shape = [number rows, number cols]
        result_array = np.full([3, 3], None)
        # Extraigo los valores de los campos enmascarados
        result_array[0,0] = field_dict[field][alfa-1,gate+1]
        result_array[0,1] = field_dict[field][alfa  ,gate+1]
        result_array[0,2] = field_dict[field][alfa+1,gate+1]
        result_array[1,0] = field_dict[field][alfa-1,gate]
        result_array[1,1] = field_dict[field][alfa  ,gate]
        result_array[1,2] = field_dict[field][alfa+1,gate]
        result_array[2,0] = field_dict[field][alfa-1,gate-1]
        result_array[2,1] = field_dict[field][alfa  ,gate-1]
        result_array[2,2] = field_dict[field][alfa+1,gate-1]

        result.append(result_array)
        
    del radar,field_dict
    return fecha,result



def radar_variable_window_lat_lon_list(radarfilepath,
                                       fields,
                                       coords_lst,
                                       lat_lst,
                                       lon_lst,
                                       rhohv_field='RHOHV',
                                       rhohv_threshold=0.8,
                                       mask=False):
    
    """
    Función para obtener los valores de variables de radar en una ventana de 3x3 celdas
    centradas sobre una serie de puntos con coordenadas lat y lon.
    La salida es un diccionario: Fecha y lista de arrays 3x3 con los valores extraidos.
    La disposición de la malla de celdas es la siguiente:
    
            | alfa-1 | alfa  | alfa+1 |
     ----------------------------------
     gate+1 | [0,0]  | [0,1] |  [0,2] |
     gate   | [1,0]  | [1,1] |  [1,2] |
     gate-1 | [2,0]  | [2,1] |  [2,2] |
     ----------------------------------
    Parameters:
            radarfilepath (str): Path al archivo volumen de radar.
            fields (str o list): Nombres de los campos de radar a extraer el valor.
            coords_lst (list): Lista con nombre de referencia de los puntos.
            lat_lst (list): Lista con la latitud en grados decimales.
            lon_lst (list): Lista con la longitud en grados decimales.
            rhohv_field (str): Nombre del campo RHOHV. Default 'RHOHV'.
            rhohv_threshold (float): Valor de RHOHV (0 a 1) para aplicar mascara. Default 0.8.
            mask (bool): True para aplicar mascara. Default False.

    Returns:
            result (dict): Dict con los resultados.
                           Contiene una clave para cada par de lat,lon ingresado. La clave
                           tiene el correspondiente nombre de la lista coords_lst.
                           Ademas tiene una clave datetime con la fecha y hora.
    """
    
    # Creamos el objeto "radar" con los campos definidos en include_fields
    radar = pyart.io.read(radarfilepath)

    # Extraemos la fecha
    fecha = radar.time['units'][14:]

    # La convertimos en un objeto datetime con el módulo datetime 
    fecha = datetime.strptime(fecha, '%Y-%m-%dT%H:%M:%SZ')
    
    field_dict = {}
    
    if mask:
        for field in fields:
            try:
                # Enmascaro los valores de los campos indicados
                field_dict[field] = np.ma.masked_where(radar.fields[rhohv_field]['data'] < rhohv_threshold,
                                                       radar.fields[field]['data'])
            except KeyError:
                print('Error. Campo no encontrado. Los campos disponibles en el volumen son',radar.fields.keys())

    else:
        for field in fields:
            try:
                field_dict[field] = radar.fields[field]['data']
            except KeyError:
                print('Error. Campo no encontrado. Los campos disponibles son:',radar.fields.keys())
            
    result = {}
    result['datetime'] = fecha
    
    for element in zip(coords_lst,lat_lst,lon_lst):
        
        # Funcion para calcular el gate y azimuth sobre las coordenadas
        #                                            Longitud    Latitud
        gate,alfa = get_nearest_gate_azimuth(radar, element[2], element[1])
        result_temp = []
        
        for field in fields:
            # np.full([shape], fillvalue)
            # shape = [number rows, number cols]
            result_array = np.full([3, 3], None)
            # Extraigo los valores de los campos enmascarados
            result_array[0,0] = field_dict[field][alfa-1,gate+1]
            result_array[0,1] = field_dict[field][alfa  ,gate+1]
            result_array[0,2] = field_dict[field][alfa+1,gate+1]
            result_array[1,0] = field_dict[field][alfa-1,gate]
            result_array[1,1] = field_dict[field][alfa  ,gate]
            result_array[1,2] = field_dict[field][alfa+1,gate]
            result_array[2,0] = field_dict[field][alfa-1,gate-1]
            result_array[2,1] = field_dict[field][alfa  ,gate-1]
            result_array[2,2] = field_dict[field][alfa+1,gate-1]

            result_temp.append(result_array)
            
        result[element[0]] = result_temp
        
    del radar,field_dict
    return result



def create_df_from_radar_windows_vars(datetime,var,fields):
    """
    Funcion para convertir el array con los valores de radar en una ventana de 3x3
    en un DataFrame cuya fila es la fecha/hora y en las columnas los valores de
    cada celda. La notación de las columnas es la misma que la función que la generó
    (radar_variable_window_lat_lon y radar_variable_window_lat_lon_list). Es decir,
            | alfa-1 | alfa  | alfa+1 |
     ----------------------------------
     gate+1 | [0,0]  | [0,1] |  [0,2] |
     gate   | [1,0]  | [1,1] |  [1,2] |
     gate-1 | [2,0]  | [2,1] |  [2,2] |
     ----------------------------------
    Parameters:
            datetime (datetime obj): Objeto datetime con la fecha y hora del volumen.
            var (tuple): Tupla con la salida de la función (radar_variable_window_lat_lon
                         o radar_variable_window_lat_lon_list).
            fields (lst o str): Lista o str de los campos extraidos en var. Debe coincidir
                                con los definidos en la función (radar_variable_window_lat_lon
                                o radar_variable_window_lat_lon_list).

    Returns:
            df (DataFrame): DataFrame de Pandas con los resultados.
    """

    # Bloque Try de control de la variable var
    for x in range(len(var)):
        try:
            dim = var[x].shape
            if dim != (3,3):
                print('Las dimensiones de las variables son distintas de (3,3)')
                sys.exit()
        except:
            raise ValueError("Las dimensiones de las variables son incorrectas")
        
    # Se crea dict vacio
    data = {}
    # Indice k para recorrer los campos de var
    k=0
    # Recorre los campos dentro de fields
    for field in fields:
        # Doble loop para recorrer la malla de celdas
        for i in [0,1,2]:
            for j in [0,1,2]:
                data[field+' ['+str(i)+','+str(j)+']'] = [var[k][i,j]]
        k=k+1

    # Se agrega la Fecha y Hora
    data['t radar[ART]'] = datetime


    # Se guarda la info en un DataFrame
    df = pd.DataFrame(data)
    # Se setea el indice como la Fecha/Hora
    df.set_index('t radar[ART]',inplace=True)
    
    return df

def db_to_linear(in_df):
    """
    Funcion para transformar los valores de las columnas de dB
    a unidades lineales: La operacion que se realiza es la siguiente:
    
                   linear_value = 10**(db_value / 10)
    Parameters: 
              in_df (DataFrame): DataFrame de entrada con los valores en db
    
    Return:
              out_df (DataFrame): DataFrame de salida en unidades lineales
    """
    
    return 10**(in_df/10)



def linear_to_db(in_df):
    """
    Funcion para transformar los valores de las columnas de unidades lineales
    a dB: La operacion que se realiza es la siguiente:
    
                      db_value = 10*Log(linear_value)
    Parameters: 
              in_df (DataFrame): DataFrame de entrada en unidades lineales
    
    Return:
              out_df (DataFrame): DataFrame de salida en db
    """
    
    return 10*np.log10(in_df)



def lat_lon_to_range_azimuth(lat,lon):
    """
    Funcion para calcular el rango y el azimuth de un punto con coordenadas
    lat y lon teniendo como referencia la ubicacion del radar RMA1.
    Parameters:
               lat (float): Latitud en grados decimales.
               lon (float): Longitud en grados decimales.
    Return:
            distance (float): Distancia al radar en km (rango).
            azimuth (float): Azimuth del punto ingresado en grados. Medido en sentido
                             horario con el 0 en el norte geográfico.
    """

    import pyproj

    # Latitud y Longitud del radar meteorológico
    radar_lat = -31.4412824015
    radar_lon = -64.1919061484
    
    # Calculo de la distancia y azimuth del segmento radar-punto
    geodesic = pyproj.Geod(ellps='WGS84')
    fwd_azimuth, back_azimuth, distance = geodesic.inv(radar_lon,
                                                       radar_lat,
                                                       lon,
                                                       lat)

    # Distancia radar-punto en kilometros
    distance = distance/1000.
    
    # Azimuth del segmento radar-punto. 0° en dirección
    # norte y aumentando en sentido horario
    if 0 <= fwd_azimuth <= 180:
        azimuth = fwd_azimuth
    elif -180 <= fwd_azimuth < 0:
        azimuth = 180 + back_azimuth

    return distance,azimuth



def calculate_R_from_Z_R(a,b,Z):
    """
    Funcion para calcular la tasa de precipitacion R a partir de la relacion Z-R.
    Parameters:
               a (float): Parametro a de la relacion Z-R.
               b (float): Parametro b de la relacion Z-R.
               Z (float): Factor de reflectividad horizontal Zh en unidades lineales.
    Return:
            R (float): Tasa de precipitacion en mm/h.
    """
    return (Z/a)**(1/b)



def calculate_R_from_Z_ZDR_R(a,b,c,Z,ZDR):
    """
    Funcion para calcular la tasa de precipitacion R a partir de la relacion Z-ZDR-R.
    Parameters:
               a (float): Parametro a de la relacion Z-R.
               b (float): Parametro b de la relacion Z-R.
               c (float): Parametro c de la relacion Z-R.
               Z (float): Factor de reflectividad horizontal Zh en unidades lineales.
               ZDR (float): Reflectividad diferencial ZDR en unidades lineales.
    Return:
            R (float): Tasa de precipitacion en mm/h.
    """
    return a*(Z**b)*(ZDR**c)

