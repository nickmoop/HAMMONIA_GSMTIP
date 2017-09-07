# coding: utf8
import os
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
    times = read_from_netCDF(file_to_read, 'TIME')
    # считаем значения сетки по долготе
    longitude_old_cells = read_from_netCDF(file_to_read, 'LON')
    # считаем значения сетки по широте, значения сетки должны быть вида:
    # север 90 градусов индекс 0, юг -90 градусов индекс i_max. т.е.
    # отсортированы с севера на юг
    latitude_old_cells = read_from_netCDF(file_to_read, 'LAT')
    # считаем данные скорости зонального и меридионального ветра,
    # температуры и молярного давления.
    zonal_speed = read_from_netCDF(file_to_read, 'u')
    meridional_speed = read_from_netCDF(file_to_read, 'v')
    temperatures = read_from_netCDF(file_to_read, 'st')
    densities = read_from_netCDF(file_to_read, 'd')

    # создадим все необходимые для записи данных папки в текущей папке
    create_directories()

    #
    latitudes_list = list(range(-90, 95, 5))
    latitudes_list.reverse()

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
        file_to_write_zwind_surfer = open(
            '{}_surfer.dat'.format(filename_zwind), 'a'
        )
        file_to_write_vwind_surfer = open(
            '{}_surfer.dat'.format(filename_vwind), 'a'
        )
        file_to_write_temperature_surfer = open(
            '{}_surfer.dat'.format(filename_temperature), 'a'
        )
        file_to_write_density_surfer = open(
            '{}_surfer.dat'.format(filename_density), 'a'
        )
        # запишем человекопонятный заголовок в Surfer файлы
        surfer_header = 'longitude   latitude   value\n'
        file_to_write_zwind_surfer.write(surfer_header)
        file_to_write_vwind_surfer.write(surfer_header)
        file_to_write_temperature_surfer.write(surfer_header)
        file_to_write_density_surfer.write(surfer_header)

        # !!!5 нужно заменить на шаг сетки!!!
        for latitude in latitudes_list:
            # узнаем индекс текущей кошироты
            # !!!5 нужно заменить на шаг сетки!!!
            latitude_i = int((90 - latitude) / 5) + 1

            # запишем время и индекс текущей кошироты
            # в заголовок блока GSMTIP файла
            header_message = '  UT={:10d} lat_i={:10d}\n   '.format(
                current_hour, latitude_i
            )
            file_to_write_temperature.write(header_message)
            file_to_write_density.write(header_message)
            file_to_write_zwind.write(header_message)
            file_to_write_vwind.write(header_message)

            # !!!5 нужно заменить на шаг сетки!!!
            for longitude_index, longitude in enumerate(range(0, 360, 5)):
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
                density = bilinear_interpolation(
                    densities[time_index][0],
                    longitude_old_cells, latitude_old_cells,
                    current_coordinates
                )

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
                surfer_coordinates = '{}   {}  '.format(longitude, latitude)
                file_to_write_zwind_surfer.write(surfer_coordinates)
                file_to_write_vwind_surfer.write(surfer_coordinates)
                file_to_write_temperature_surfer.write(surfer_coordinates)
                file_to_write_density_surfer.write(surfer_coordinates)
                # запишем данные в Surfer файлы
                file_to_write_zwind_surfer.write(zwind_data + '\n')
                file_to_write_vwind_surfer.write(vwind_data + '\n')
                file_to_write_temperature_surfer.write(temperature_data + '\n')
                file_to_write_density_surfer.write(density_data + '\n')


def create_directories():
    """
    Метод для создания всех необходимых для записи папок.
    Создает папки с нужными именами в текущей папке
    """

    # создадим список со всеми именами всех нужных папок
    all_directories = [
        'GSMTIP/Zwind', 'GSMTIP/Vwind',
        'GSMTIP/Temperature', 'GSMTIP/Density'
    ]
    for directory in all_directories:
        # проверим нужно существует ли папка
        if not os.path.exists(directory):
            os.makedirs(directory)


def calculate_air_density_from_mm_air(mm_air, temperature):
    """
    Метод для преобразования молярной плотности в "обычную" плотность,
    на высоте 80км

    :param mm_air: молярная плотность, г/моль
    :param temperature: температура, кельвинах
    :return: плотность газа на высоте 80км, г/см3
    """

    return mm_air * PRESSURE / temperature / GAS_CONSTANT / 1000000


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

    # простая замена, если текущие координаты выходят за границы сетки HAMMONIA
    if not latitude_indexes or not longitude_indexes:
        if latitude <= latitude_cells[-1]:
            latitude_index = -1
        if latitude >= latitude_cells[0]:
            latitude_index = 0
        if longitude == 0:
            longitude_index = 0
        if longitude >= longitude_indexes[-1]:
            longitude_index = -1
        return values[latitude_index][longitude_index]

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
        if (old_cells[index] <= coordinate <= old_cells[index + 1] or
                        old_cells[index] >= coordinate >= old_cells[index + 1]
        ):
            return index, index + 1

    # вернем None, если текущие координаты не вписываются в исходную сетку
    return None


# если запускаем этот файл.py а не импортируем из него методы
if __name__ == '__main__':
    # адрес до netcdf файла с данными HAMMONIA
    path_to_hammonia = 'HAMMONIA_80km_200901.nc'

    # создадим на основе данных HAMMONIA текстовые файлы, такого формата,
    # чтобы использовать их в модели GSMTIP
    create_gsmtip_files_from_hammonia(path_to_hammonia)
