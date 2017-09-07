# coding: utf8
import glob
from collections import OrderedDict
from datetime import datetime, timedelta

from netCDF4 import Dataset

# Давление на высоте 80км
PRESSURE = 101325 * 2.7 ** (80000 / -7000)
# Универсальная газовая постоянная
GAS_CONSTANT = 8.3144598
# стартовое время от которого в HAMMONIA (GSMTIP) производится отсчет времени.
# в файле с данными HAMMONIA (GSMTIP) указана разница во времени
# между ЭТИМ моментом и моментом проведения измерений (моделирования?)
START_DATE_HAMMONIA = datetime.strptime(
    '2007-12-31 00:00:00', "%Y-%m-%d %H:%M:%S"
)
START_DATE_GSMTIP = datetime.strptime(
    '2009-01-01 00:00:00', "%Y-%m-%d %H:%M:%S"
)

LONGITUDE_VARIABLE_PARAMETERS = OrderedDict()
LONGITUDE_VARIABLE_PARAMETERS['standard_name'] = 'longitude'
LONGITUDE_VARIABLE_PARAMETERS['long_name'] = 'longitude'
LONGITUDE_VARIABLE_PARAMETERS['units'] = 'degrees_east'
LONGITUDE_VARIABLE_PARAMETERS['axis'] = 'X'
LATITUDE_VARIABLE_PARAMETERS = OrderedDict()
LATITUDE_VARIABLE_PARAMETERS['standard_name'] = 'latitude'
LATITUDE_VARIABLE_PARAMETERS['long_name'] = 'latitude'
LATITUDE_VARIABLE_PARAMETERS['units'] = 'degrees_north'
LATITUDE_VARIABLE_PARAMETERS['axis'] = 'Y'
HEIGHT_VARIABLE_PARAMETERS = OrderedDict()
HEIGHT_VARIABLE_PARAMETERS['standard_name'] = 'height'
HEIGHT_VARIABLE_PARAMETERS['long_name'] = 'height'
HEIGHT_VARIABLE_PARAMETERS['units'] = 'm'
HEIGHT_VARIABLE_PARAMETERS['positive'] = 'up'
HEIGHT_VARIABLE_PARAMETERS['axis'] = 'Z'
TIME_VARIABLE_PARAMETERS = OrderedDict()
TIME_VARIABLE_PARAMETERS['standard_name'] = 'time'
TIME_VARIABLE_PARAMETERS['long_name'] = 'time'
TIME_VARIABLE_PARAMETERS['units'] = 'days since 2007-12-31 00:00:00'
TIME_VARIABLE_PARAMETERS['calendar'] = 'standard'
MM_AIR_VARIABLE_PARAMETERS = OrderedDict()
MM_AIR_VARIABLE_PARAMETERS['long_name'] = 'mol. air mass'
MM_AIR_VARIABLE_PARAMETERS['units'] = 'g/mol'
MM_AIR_VARIABLE_PARAMETERS['code'] = '80'
MM_AIR_VARIABLE_PARAMETERS['table'] = '128'
MM_AIR_VARIABLE_PARAMETERS['grid_type'] = 'gaussian'
ST_VARIABLE_PARAMETERS = OrderedDict()
ST_VARIABLE_PARAMETERS['long_name'] = 'temperature'
ST_VARIABLE_PARAMETERS['units'] = 'K'
ST_VARIABLE_PARAMETERS['code'] = '130'
ST_VARIABLE_PARAMETERS['table'] = '128'
ST_VARIABLE_PARAMETERS['grid_type'] = 'gaussian'
V_VARIABLE_PARAMETERS = OrderedDict()
V_VARIABLE_PARAMETERS['long_name'] = 'v-velocity'
V_VARIABLE_PARAMETERS['units'] = 'm/s'
V_VARIABLE_PARAMETERS['code'] = '132'
V_VARIABLE_PARAMETERS['table'] = '128'
V_VARIABLE_PARAMETERS['grid_type'] = 'gaussian'
U_VARIABLE_PARAMETERS = OrderedDict()
U_VARIABLE_PARAMETERS['long_name'] = 'u-velocity'
U_VARIABLE_PARAMETERS['units'] = 'm/s'
U_VARIABLE_PARAMETERS['code'] = '131'
U_VARIABLE_PARAMETERS['table'] = '128'
U_VARIABLE_PARAMETERS['grid_type'] = 'gaussian'

DESCRIPTION_VARIABLES_PARAMETERS = {
    'lon': LONGITUDE_VARIABLE_PARAMETERS,
    'lat': LATITUDE_VARIABLE_PARAMETERS,
    'height': HEIGHT_VARIABLE_PARAMETERS,
    'time': TIME_VARIABLE_PARAMETERS,
    'mm_air': MM_AIR_VARIABLE_PARAMETERS,
    'st': ST_VARIABLE_PARAMETERS,
    'v': V_VARIABLE_PARAMETERS,
    'u': U_VARIABLE_PARAMETERS,
}


def create_hammonia_from_gsmtip(gsmtip_folder):
    gsm_tip_data_dict = create_gsmtip_data_dict_from_gsmtip_files(
        gsmtip_folder)
    write_gsmtip_data_dict_to_hammonia_netcdf_file(gsm_tip_data_dict)


def write_gsmtip_data_dict_to_hammonia_netcdf_file(gsm_tip_data_dict):
    filename = 'tmp.nc'
    converted_keys = create_dimensions_in_netcdf(gsm_tip_data_dict, filename)

    file_netCDF = Dataset(filename, 'a', format='NETCDF4')
    mm_air = file_netCDF.createVariable('mm_air', 'f4',
                                        ('time', 'height', 'lat', 'lon',))
    st = file_netCDF.createVariable('st', 'f4',
                                    ('time', 'height', 'lat', 'lon',))
    u = file_netCDF.createVariable('u', 'f4',
                                   ('time', 'height', 'lat', 'lon',))
    v = file_netCDF.createVariable('v', 'f4',
                                   ('time', 'height', 'lat', 'lon',))

    update_variable_parameters('mm_air', file_netCDF)
    update_variable_parameters('st', file_netCDF)
    update_variable_parameters('u', file_netCDF)
    update_variable_parameters('v', file_netCDF)

    for key, gsmtip_values in gsm_tip_data_dict.items():
        TMP_n_dimensional_array = [[] for i in range(0, len(
            converted_keys['time']['all']))]
        for time_index, time_value in enumerate(converted_keys['time']['all']):
            TMP_n_dimensional_array[time_index] = [
                [] for i in range(0, len(converted_keys['height']['all']))
                ]
            converted_time = converted_keys['time'][time_value]
            for height_index, height_value in enumerate(
                    converted_keys['height']['all']):
                TMP_n_dimensional_array[time_index][height_index] = [
                    [] for i in range(0, len(converted_keys['lat']['all']))
                    ]
                converted_height = converted_keys['height'][height_value]
                for latitude_index, latitude_value in enumerate(
                        converted_keys['lat']['all']):
                    TMP_n_dimensional_array[time_index][height_index][
                        latitude_index] = [
                        [] for i in range(0, len(converted_keys['lon']['all']))
                        ]
                    converted_latitude = converted_keys['lat'][latitude_value]
                    for longitude_index, longitude_value in enumerate(
                            converted_keys['lon']['all']):
                        converted_longitude = converted_keys['lon'][
                            longitude_value]
                        if key == 'Temperature':
                            value = \
                            gsmtip_values[converted_time][converted_height][
                                converted_latitude][converted_longitude]
                            converted_value = value
                        elif key == 'Density':
                            value = \
                            gsmtip_values[converted_time][converted_height][
                                converted_latitude][converted_longitude]
                            temperature = \
                            gsm_tip_data_dict['Temperature'][converted_time][
                                converted_height][converted_latitude][
                                converted_longitude]
                            converted_value = calculate_mm_air_from_air_density(
                                value, temperature)
                        elif key == 'Vwind':
                            value = \
                            gsmtip_values[converted_time][converted_height][
                                converted_latitude][converted_longitude]
                            converted_value = -value
                        elif key == 'Zwind':
                            value = \
                            gsmtip_values[converted_time][converted_height][
                                converted_latitude][converted_longitude]
                            converted_value = value

                        TMP_n_dimensional_array[time_index][height_index][
                            latitude_index][longitude_index] = converted_value

        if key == 'Temperature':
            st[:] = TMP_n_dimensional_array
        elif key == 'Density':
            mm_air[:] = TMP_n_dimensional_array
        elif key == 'Vwind':
            v[:] = TMP_n_dimensional_array
        elif key == 'Zwind':
            u[:] = TMP_n_dimensional_array

    file_netCDF.close()


def create_dimensions_in_netcdf(gsm_tip_data_dict, filename):
    for parameter_name in gsm_tip_data_dict:
        time_keys = list(gsm_tip_data_dict[parameter_name].keys())
        altitude_keys = list(
            gsm_tip_data_dict[parameter_name][time_keys[0]].keys())
        latitude_keys = list(gsm_tip_data_dict[parameter_name][time_keys[0]][
                                 altitude_keys[0]].keys())
        longitude_keys = list(
            gsm_tip_data_dict[parameter_name][time_keys[0]][altitude_keys[0]][
                latitude_keys[0]].keys())
        break

    converted_keys = convert_keys(time_keys, altitude_keys, latitude_keys,
                                  longitude_keys)
    file_netCDF = Dataset(filename, 'w', format='NETCDF4')
    all_keys = ['lon', 'lat', 'height', 'time']
    for key in all_keys:
        data = list(converted_keys[key]['all'])
        fill_dimension(file_netCDF, key, data, typ='f8')
        update_variable_parameters(key, file_netCDF)

    file_netCDF.close()

    return converted_keys


def update_variable_parameters(variable_name, file_netCDF):
    description_parameters = DESCRIPTION_VARIABLES_PARAMETERS[variable_name]
    for parameter_name in description_parameters:
        setattr(file_netCDF.variables[variable_name], parameter_name,
                description_parameters[parameter_name])


def convert_keys(time_keys, altitude_keys, latitude_keys, longitude_keys):
    converted_keys = dict()
    converted_keys['time'] = {}
    converted_keys['time']['all'] = []
    converted_keys['height'] = {}
    converted_keys['height']['all'] = []
    converted_keys['lat'] = {}
    converted_keys['lat']['all'] = []
    converted_keys['lon'] = {}
    converted_keys['lon']['all'] = []

    for key in time_keys:
        date_delta_hammonia = key - START_DATE_HAMMONIA
        hours_delta = date_delta_hammonia.seconds / 86400
        value = date_delta_hammonia.days + hours_delta
        converted_keys['time'][value] = key
        converted_keys['time']['all'].append(value)

    for key in altitude_keys:
        converted_keys['height'][80] = key
        converted_keys['height']['all'].append(80)

    for key in latitude_keys:
        value = key * 5 - 90
        converted_keys['lat'][value] = key
        converted_keys['lat']['all'].append(value)

    for key in longitude_keys:
        value = key * 5
        converted_keys['lon'][value] = key
        converted_keys['lon']['all'].append(value)

    converted_keys['time']['all'].sort()
    converted_keys['height']['all'].sort()
    converted_keys['lon']['all'].sort()
    converted_keys['lat']['all'].sort()

    return converted_keys


def fill_dimension(file_netCDF, name, data, typ='f4'):
    file_netCDF.createDimension(name, len(data))
    file_netCDF.createVariable(name, typ, (name,))
    file_netCDF.variables[name][:] = data


def create_gsmtip_data_dict_from_gsmtip_files(gsmtip_folder):
    parameters_names = ['Density', 'Temperature', 'Vwind', 'Zwind']
    gsmtip_data_dict = {}
    for parameter_name in parameters_names:
        parameter_folder = '{}{}'.format(gsmtip_folder, parameter_name)
        for file_name in glob.glob('{}/*'.format(parameter_folder)):
            file_to_read = open(file_name, 'r')
            lines = file_to_read.readlines()
            file_to_read.close()
            current_file_day = int(file_name.split('_')[1])
            base_days_delta = timedelta(days=current_file_day)
            longitude_i = 0
            for line in lines:
                if 'UT' in line:
                    splitted_line = line.split()
                    hour_UT = int(splitted_line[1])
                    latitude_i = int(splitted_line[3])
                    longitude_i = 0
                    days_delta = base_days_delta + timedelta(hours=hour_UT)
                    gsmtip_datetime = START_DATE_GSMTIP + days_delta
                    continue
                else:
                    splitted_line = line.split()
                    for parameter_value in splitted_line:
                        longitude_i += 1
                        write_value_to_gsmtip_dict(
                            gsmtip_data_dict, parameter_name, gsmtip_datetime,
                            latitude_i, longitude_i, parameter_value
                        )

    return gsmtip_data_dict


def write_value_to_gsmtip_dict(gsmtip_data_dict, parameter_name,
                               gsmtip_datetime, latitude_i, longitude_i,
                               parameter_value):
    altitude_index = 0
    if parameter_name not in gsmtip_data_dict:
        gsmtip_data_dict[parameter_name] = {}

    if gsmtip_datetime not in gsmtip_data_dict[parameter_name]:
        gsmtip_data_dict[parameter_name][gsmtip_datetime] = {}

    if altitude_index not in gsmtip_data_dict[parameter_name][gsmtip_datetime]:
        gsmtip_data_dict[parameter_name][gsmtip_datetime][altitude_index] = {}

    if latitude_i not in gsmtip_data_dict[parameter_name][gsmtip_datetime][
        altitude_index]:
        gsmtip_data_dict[parameter_name][gsmtip_datetime][altitude_index][
            latitude_i] = {}

    if longitude_i not in \
            gsmtip_data_dict[parameter_name][gsmtip_datetime][altitude_index][
                latitude_i]:
        gsmtip_data_dict[parameter_name][gsmtip_datetime][altitude_index][
            latitude_i][longitude_i] = float(parameter_value)


def calculate_mm_air_from_air_density(air_density, temperature):
    return air_density * 1000000 * GAS_CONSTANT * temperature / PRESSURE


# если запускаем этот файл.py а не импортируем из него методы
if __name__ == '__main__':
    # адрес до папки с данными GSMTIP
    gsmtip_folder = '/home/nick/code/python/netCDF/GSMTIP/'

    # создадим на основе данных GSMTIP netcdf-файлы, такого формата,
    # чтобы использовать их в модели HAMMONIA
    create_hammonia_from_gsmtip(gsmtip_folder)
