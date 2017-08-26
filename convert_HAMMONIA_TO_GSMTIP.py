import glob
from datetime import datetime, timedelta
from collections import OrderedDict
from netCDF4 import Dataset

# Давление на высоте 80км
PRESSURE = 101325 * 2.7 ** (80000 / -7000)
# Универсальная газовая постоянная
GAS_CONSTANT = 8.3144598
# стартовое время от которого в HAMMONIA (GSMTIP) производится отсчет времени.
# в файле с данными HAMMONIA (GSMTIP) указана разница во времени
# между ЭТИМ моментом и моментом проведения измерений (моделирования?)
START_DATE_HAMMONIA = datetime.strptime('2007-12-31 00:00:00', "%Y-%m-%d %H:%M:%S")
START_DATE_GSMTIP = datetime.strptime('2009-01-01 00:00:00', "%Y-%m-%d %H:%M:%S")

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


def create_gsmtip_files_from_hammonia(file_to_read):
    """
    Метод для преобразования данных HAMMONIA (netcdf файлы определенного
    формата. не забыть про возможные изменения названий переменных и изменения
    сеток) в текстовые файлы для использования в GSMTIP модели. Создает новые
    файлы температуры, зонального, меридионального ветров и плотностью газа
    на высоте 80км. Пути по которым создаются файлы, см. ниже
    (сделать задаваемую папку для выходных GSMTIP файлов и
    передавать её как аргумент?).

    :param file_to_read: natcdf файл с данными HAMMONIA
    """

    # считаем значения сетки по времени
    times = read_from_netCDF(file_to_read, 'time')
    # считаем значения сетки по долготе
    longitude_old_cells = read_from_netCDF(file_to_read, 'lon')
    # считаем значения сетки по широте и приведем к удобному виду:
    # север-90 градусов индекс-0, юг -90 градусов индекс-i_max. т.е.
    # упорядочим индексы с севера на юг
    latitude_old_cells = list(read_from_netCDF(file_to_read, 'lat'))
    latitude_old_cells.reverse()
    # считаем данные скорости зонального и меридионального ветра,
    # температуры и молярного давления.
    zonal_speed = read_from_netCDF(file_to_read, 'u')
    meridional_speed = read_from_netCDF(file_to_read, 'v')
    temperatures = read_from_netCDF(file_to_read, 'st')
    mm_air = read_from_netCDF(file_to_read, 'mm_air')

    # итерируемся по всем доступным индексам времени измерения (моделирования?)
    for time_index in range(0, len(times)):
        # получим разницу между начальным моментом времени
        # и временем проведения измерения (моделирования?)
        time_delta = timedelta(days=times[time_index])
        # рассчитаем в какое время и в какой день было проведенено измерение
        # (моделирование?)
        corrected_date = START_DATE_HAMMONIA + time_delta
        # получим номер дня в году, для формирования имени GSMTIP файла
        day_of_year = corrected_date.timetuple().tm_yday
        # получим значение часа для записи в заголовок GSMTIP файла
        current_hour = corrected_date.hour
        # дадим имена суточным GSMTIP файлам в зависимости от типа данных
        # скорость зонального ветра
        filename_zwind = 'GSMTIP/Zwind/Zwind_{}'.format(day_of_year)
        # скорость меридионального ветра
        filename_vwind = 'GSMTIP/Vwind/Vwind_{}'.format(day_of_year)
        # температура
        filename_temperature = 'GSMTIP/Temperature/T_{}'.format(day_of_year)
        # плотность
        filename_density = 'GSMTIP/Density/PL_{}'.format(day_of_year)

        # создадим GSMTIP файлы для записи
        file_to_write_zwind = open(filename_zwind, 'a')
        file_to_write_vwind = open(filename_vwind, 'a')
        file_to_write_temperature = open(filename_temperature, 'a')
        file_to_write_density = open(filename_density, 'a')
        # создадим Surfer файлы для записи
        # нужно для проверки преобразований и считывания данных
        file_to_write_zwind_surfer = open('{}_surfer'.format(filename_zwind),'a')
        file_to_write_vwind_surfer = open('{}_surfer'.format(filename_vwind),'a')
        file_to_write_temperature_surfer = open('{}_surfer'.format(filename_temperature), 'a')
        file_to_write_density_surfer = open('{}_surfer'.format(filename_density), 'a')
        # запишем человекопонятный заголовок в Surfer файлы
        surfer_header = 'longitude   latitude   value\n'
        file_to_write_zwind_surfer.write(surfer_header)
        file_to_write_vwind_surfer.write(surfer_header)
        file_to_write_temperature_surfer.write(surfer_header)
        file_to_write_density_surfer.write(surfer_header)

        # !!!5 нужно заменить на шаг сетки!!!
        for latitude in range(-85, 95, 5):
            # узнаем индекс текущей кошироты
            # !!!5 нужно заменить на шаг сетки!!!
            latitude_i = int((90 - latitude) / 5)

            # запишем время и индекс текущей кошироты
            # в заголовок блока GSMTIP файла
            header_message = '  UT={:10d} lat_i={:10d}\n   '.format(
                current_hour, latitude_i
            )
            file_to_write_temperature.write(header_message)
            file_to_write_density.write(header_message)
            file_to_write_zwind.write(header_message)
            file_to_write_vwind.write(header_message)

            # координаты центра сетки, используется для интерполяции
            # !!!2.5 нужно заменить на половину шага сетки!!!
            latitude -= 2.5

            # !!!5 нужно заменить на шаг сетки!!!
            for longitude_index, longitude in enumerate(range(5, 365, 5)):
                # !!!2.5 нужно заменить на половину шага сетки!!!
                longitude -= 2.5
                # сформируем пару координат для интерполяции
                current_coordinates = [longitude, latitude]

                # интерполируем значения параметров для новой ячейки сетки
                # для зонального ветра
                zonal_wind = bilinear_interpolation(
                    zonal_speed[time_index][0],
                    longitude_old_cells, latitude_old_cells,
                    current_coordinates
                )
                # для меридионального ветра
                meridional_wind = bilinear_interpolation(
                    meridional_speed[time_index][0],
                    longitude_old_cells, latitude_old_cells,
                    current_coordinates
                )
                # для температуры
                temperature = bilinear_interpolation(
                    temperatures[time_index][0],
                    longitude_old_cells, latitude_old_cells,
                    current_coordinates
                )
                # для плоотности
                molar_density = bilinear_interpolation(
                    mm_air[time_index][0],
                    longitude_old_cells, latitude_old_cells,
                    current_coordinates
                )
                # преобразуем молярную плотность в обычную, на высоте 80км
                density = calculate_air_density_from_mm_air(molar_density, temperature)

                # отформатируем данные как в GSMTIP
                # 2 пробела между значениями в строке, 4 знака после запятой,
                # экспоненциальный (научный) формат цифр
                temperature_data = '  {:.4E}'.format(temperature)
                density_data = '  {:.4E}'.format(density)
                zwind_data = '  {:.4E}'.format(zonal_wind)
                vwind_data = '  {:.4E}'.format(-meridional_wind)

                # узнаем не закончилась ли строка в блоке.
                # в GSMTIP пишется 6 значений в одну строку,
                # потом переход на новую строку
                new_block_or_string = None
                longitude_index += 1
                if longitude_index % 6 == 0:
                    new_block_or_string = '\n'
                    # когда заканчивается блок делаются дополнительные пробелы
                    # в начале строки с заголовком нового блока
                    # !!!72 зависит от количества ячеек долготной сетки!!!
                    if longitude_index != 72:
                        new_block_or_string += '   '

                # запишем отформатированные данные в GSMTIP файл
                file_to_write_temperature.write(temperature_data)
                file_to_write_density.write(density_data)
                file_to_write_zwind.write(zwind_data)
                file_to_write_vwind.write(vwind_data)

                # применим необходимое форматирование к файлу
                # если закончился блок или строка с данными
                if new_block_or_string:
                    file_to_write_temperature.write(new_block_or_string)
                    file_to_write_density.write(new_block_or_string)
                    file_to_write_zwind.write(new_block_or_string)
                    file_to_write_vwind.write(new_block_or_string)

                # запишем текущие координаты в Surfer файлы
                surfer_coordinates = '\n{}   {} '.format(longitude, latitude)
                file_to_write_zwind_surfer.write(surfer_coordinates)
                file_to_write_vwind_surfer.write(surfer_coordinates)
                file_to_write_temperature_surfer.write(surfer_coordinates)
                file_to_write_density_surfer.write(surfer_coordinates)
                # запишем данные в Surfer файлы
                file_to_write_zwind_surfer.write(zwind_data)
                file_to_write_vwind_surfer.write(vwind_data)
                file_to_write_temperature_surfer.write(temperature_data)
                file_to_write_density_surfer.write(density_data)


def create_hammonia_from_gsmtip(gsmtip_folder):
    gsm_tip_data_dict = create_gsmtip_data_dict_from_gsmtip_files(gsmtip_folder)
    write_gsmtip_data_dict_to_hammonia_netcdf_file(gsm_tip_data_dict)


def write_gsmtip_data_dict_to_hammonia_netcdf_file(gsm_tip_data_dict):
    filename = 'tmp.nc'
    converted_keys = create_dimensions_in_netcdf(gsm_tip_data_dict, filename)

    file_netCDF = Dataset(filename, 'a', format='NETCDF4')
    mm_air = file_netCDF.createVariable('mm_air', 'f4', ('time', 'height', 'lat', 'lon', ))
    st = file_netCDF.createVariable('st', 'f4', ('time', 'height', 'lat', 'lon',))
    u = file_netCDF.createVariable('u', 'f4', ('time', 'height', 'lat', 'lon',))
    v = file_netCDF.createVariable('v', 'f4', ('time', 'height', 'lat', 'lon',))

    update_variable_parameters('mm_air', file_netCDF)
    update_variable_parameters('st', file_netCDF)
    update_variable_parameters('u', file_netCDF)
    update_variable_parameters('v', file_netCDF)

    for key, gsmtip_values in gsm_tip_data_dict.items():
        TMP_n_dimensional_array = [[] for i in range(0, len(converted_keys['time']['all']))]
        for time_index, time_value in enumerate(converted_keys['time']['all']):
            TMP_n_dimensional_array[time_index] = [
                [] for i in range(0, len(converted_keys['height']['all']))
            ]
            converted_time = converted_keys['time'][time_value]
            for height_index, height_value in enumerate(converted_keys['height']['all']):
                TMP_n_dimensional_array[time_index][height_index] = [
                    [] for i in range(0, len(converted_keys['lat']['all']))
                ]
                converted_height = converted_keys['height'][height_value]
                for latitude_index, latitude_value in enumerate(converted_keys['lat']['all']):
                    TMP_n_dimensional_array[time_index][height_index][latitude_index] = [
                        [] for i in range(0, len(converted_keys['lon']['all']))
                    ]
                    converted_latitude = converted_keys['lat'][latitude_value]
                    for longitude_index, longitude_value in enumerate(converted_keys['lon']['all']):
                        converted_longitude = converted_keys['lon'][longitude_value]
                        if key == 'Temperature':
                            value = gsmtip_values[converted_time][converted_height][converted_latitude][converted_longitude]
                            converted_value = value
                        elif key == 'Density':
                            value = gsmtip_values[converted_time][converted_height][converted_latitude][converted_longitude]
                            temperature = gsm_tip_data_dict['Temperature'][converted_time][converted_height][converted_latitude][converted_longitude]
                            converted_value = calculate_mm_air_from_air_density(value, temperature)
                        elif key == 'Vwind':
                            value = gsmtip_values[converted_time][converted_height][converted_latitude][converted_longitude]
                            converted_value = -value
                        elif key == 'Zwind':
                            value = gsmtip_values[converted_time][converted_height][converted_latitude][converted_longitude]
                            converted_value = value

                        TMP_n_dimensional_array[time_index][height_index][latitude_index][longitude_index] = converted_value

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
        altitude_keys = list(gsm_tip_data_dict[parameter_name][time_keys[0]].keys())
        latitude_keys = list(gsm_tip_data_dict[parameter_name][time_keys[0]][altitude_keys[0]].keys())
        longitude_keys = list(gsm_tip_data_dict[parameter_name][time_keys[0]][altitude_keys[0]][latitude_keys[0]].keys())
        break

    converted_keys = convert_keys(time_keys, altitude_keys, latitude_keys, longitude_keys)
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
        setattr(file_netCDF.variables[variable_name], parameter_name, description_parameters[parameter_name])


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
        hours_delta = date_delta_hammonia.seconds/86400
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
    file_netCDF.createVariable(name, typ, (name, ))
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


def write_value_to_gsmtip_dict(gsmtip_data_dict, parameter_name, gsmtip_datetime, latitude_i, longitude_i, parameter_value):
    altitude_index = 0
    if parameter_name not in gsmtip_data_dict:
        gsmtip_data_dict[parameter_name] = {}

    if gsmtip_datetime not in gsmtip_data_dict[parameter_name]:
        gsmtip_data_dict[parameter_name][gsmtip_datetime] = {}

    if altitude_index not in gsmtip_data_dict[parameter_name][gsmtip_datetime]:
        gsmtip_data_dict[parameter_name][gsmtip_datetime][altitude_index] = {}

    if latitude_i not in gsmtip_data_dict[parameter_name][gsmtip_datetime][altitude_index]:
        gsmtip_data_dict[parameter_name][gsmtip_datetime][altitude_index][latitude_i] = {}

    if longitude_i not in gsmtip_data_dict[parameter_name][gsmtip_datetime][altitude_index][latitude_i]:
        gsmtip_data_dict[parameter_name][gsmtip_datetime][altitude_index][latitude_i][longitude_i] = float(parameter_value)


def calculate_air_density_from_mm_air(mm_air, temperature):
    """
    Метод для преобразования молярной плотности в "обычную" плотность,
    на высоте 80км

    :param mm_air: молярная плотность, г/моль
    :param temperature: температура, кельвинах
    :return: плотность газа на высоте 80км, г/см3
    """

    return mm_air * PRESSURE / temperature / GAS_CONSTANT / 1000000


def calculate_mm_air_from_air_density(air_density, temperature):
    return air_density * 1000000 * GAS_CONSTANT * temperature / PRESSURE


def read_from_netCDF(filename, parameter_name):
    """
    Метод для чтения переменных из netCDF файла по имени переменных

    :param filename: имя или полный путь netCDF файла из которого будем читать
    :param parameter_name: название параметра который хотим считать
    :return: значения считываемого параметра в виде списка
    """

    file_to_read = Dataset(filename, 'r', format='NETCDF3')
    values = file_to_read.variables[parameter_name][:]
    file_to_read.close()

    return values


def bilinear_interpolation(
        values, longitude_cells, latitude_cells, current_coordinates
):
    """
    Метод для проведения билинейной интерполяции. Нужен для перевода одной
    размерности сетки в другую. Например 10х10 в 5х5.

    :param values: двумерный список, первый индекс широта, второй долгота
    :param longitude_cells: исходная сетка по долготе,
        содержит список долгот в которых есть значения каких-либо параметров
    :param latitude_cells: исходная сетка по широте
        содержит список широт в которых есть значения каких-либо параметров
    :param current_coordinates: координаты текущей ячейки сетки,
                    для которой хотим получить интерполированное значение
    :return: интерполированное значение в центре квадрата,
        образованного ближайшими значениями долготы и широты
        к текущим координатам
    """

    # получим значения текущих долготы и широты
    longitude = current_coordinates[0]
    latitude = current_coordinates[1]

    # найдем индексы, ближайших к текущим коррдинатам,
    # широты и долготы в старой сетке,
    # для построения "квадрата" вокруг текущих координат
    latitude_indexes = find_near_indexes(latitude_cells, latitude)
    longitude_indexes = find_near_indexes(longitude_cells, longitude)

    # определим значения в углах квадрата
    top_left = values[latitude_indexes[1]][longitude_indexes[0]]
    top_right = values[latitude_indexes[1]][longitude_indexes[1]]
    bottom_left = values[latitude_indexes[0]][longitude_indexes[0]]
    bottom_right = values[latitude_indexes[0]][longitude_indexes[1]]

    # посчитаем интерполированные верхнее и нижнее значения
    # интерполируем левые/правые значения сверху и снизу
    # получаем значения вверху/внизу квадрата, посередине ребра
    denominator = longitude_cells[longitude_indexes[1]] - longitude_cells[
        longitude_indexes[0]]
    numerator_right_side = longitude_cells[longitude_indexes[1]] - longitude
    numerator_left_side = longitude - longitude_cells[longitude_indexes[0]]
    interpolated_top = numerator_right_side / denominator * top_left
    interpolated_top += numerator_left_side / denominator * top_right
    interpolated_bottom = numerator_right_side / denominator * bottom_left
    interpolated_bottom += numerator_left_side / denominator * bottom_right

    # интерполируем верхнее и нижнее значение в среднее, между верхом и низом
    # получаем искомое значение в точке пространства,
    # которая находится по текущим координатам
    denominator = latitude_cells[latitude_indexes[1]] - latitude_cells[
        latitude_indexes[0]]
    numerator_top_side = latitude_cells[latitude_indexes[1]] - latitude
    numerator_bottom_side = latitude - latitude_cells[latitude_indexes[0]]
    converted_value = numerator_top_side / denominator * interpolated_bottom
    converted_value += numerator_bottom_side / denominator * interpolated_top

    return converted_value


def find_near_indexes(old_cells, coordinate):
    """
    Метод для нахождения индексов ближайших координат для указанной сетки и
    координат. Например есть сетка [0, 5, 10, 15] и координаты 7. Ближайшими
    координатами будут 5, 10 (в сетке), и их индексы 1, 2. Нужно для
    формирования "квадрата" вокруг указанных координат, для проведения
    билинейной интерполяции и преобразования шага сетки в произвольный

    :param old_cells: список координат исходной сетки
    :param coordinate: координаты для которых хоти найти ближайшие индексы
    :return: ближайшие индексы для указанных координат для исходной сетки
    """

    # итерируемся по всем индексам кроме последнего (см. ниже)
    for index in range(0, len(old_cells) - 1):

        # если значение координат находится в интервале
        # между текущим и следующим значением сетки, вернем индексы
        if old_cells[index] <= coordinate <= old_cells[index + 1]:
            return index, index + 1

    # вернем нули (а может лучше None?),
    # если текущие координаты не вписываются в исходную сетку
    return 0, 0


# если запускаем этот файл.py а не импортируем из него методы
if __name__ == '__main__':
    # адрес до netcdf файла с данными HAMMONIA
    path_to_hammonia = '/media/nick/Elements/PanOply/HAMMONIA_test_200901_6h.nc'

    # создадим на основе данных HAMMONIA текстовые файлы, такого формата,
    # чтобы использовать их в модели GSMTIP
    # create_gsmtip_files_from_hammonia(path_to_hammonia)

    # адрес до папки с данными GSMTIP
    gsmtip_folder = '/home/nick/code/python/netCDF/GSMTIP/'

    # создадим на основе данных GSMTIP netcdf-файлы, такого формата,
    # чтобы использовать их в модели HAMMONIA
    create_hammonia_from_gsmtip(gsmtip_folder)
