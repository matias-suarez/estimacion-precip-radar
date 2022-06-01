# estimacion-precip-radar

Funciones para tratamiento de datos de radar elaboradas para el trabajo final de Sensoramiento Remoto 2 (UBA)

Este modulo contiene las siguientes funciones:

1) get_nearest_gate_azimuth: Calcula gate y azimuth dadas las coordenadas geográficas.
2) radar_variable_lat_lon: Función para obtener los valores de variables de radar sobre un punto con coordenadas lat y lon.
3) radar_variable_window_lat_lon: Función para obtener los valores de variables de radar en una ventana de 3x3 celdas centradas sobre un punto con coordenadas lat y lon.
4) radar_variable_window_lat_lon_list: Función para obtener los valores de variables de radar en una ventana de 3x3 celdas centradas sobre una serie de puntos con coordenadas lat y lon.
5) create_df_from_radar_windows_vars: Funcion para convertir el array con los valores de radar en una ventana de 3x3 en un DataFrame cuya fila es la fecha/hora y en las columnas los valores de cada celda.
6) db_to_linear: Funcion para transformar los valores de las columnas de dB a unidades lineales.
7) linear_to_db: Funcion para transformar los valores de las columnas de unidades lineales a dB.
8) lat_lon_to_range_azimuth: Funcion para calcular el rango y el azimuth de un punto con coordenadas lat y lon teniendo como referencia la ubicacion del radar.
9) calculate_R_from_Z_R: Funcion para calcular la tasa de precipitacion R a partir de la relacion Z-R.
10) calculate_R_from_Z_ZDR_R: Funcion para calcular la tasa de precipitacion R a partir de la relacion Z-ZDR-R.
