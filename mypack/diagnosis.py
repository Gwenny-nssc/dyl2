from metpy.calc import dewpoint_from_relative_humidity
from metpy.units import units
import numpy as np
import xarray as xr
import os
from datetime import datetime
import metpy.constants as constants
import metpy.calc as mpcalc

##输入np数组,其中t的单位要求为开尔文温度，rh单位为%要求两者的维度保持一致
def dtp(t,rh):
    t = t * units("K")
    rh = rh * units("%")
    a = dewpoint_from_relative_humidity(t, rh)
    return a

##需求输入开尔文温度
def a_index(t850,td850,t700,td700,t500,td500):
    '''计算A指数'''
    a = (t850 - t500) - (t850 - td850) - (t700 - td700) - (t500 - td500)
    return a

#计算单层水汽通量和水汽通量散度
#输入的均为np数组,返回水汽通量和水汽通量散度单位分别是通量kg/(m*s) 散度kg/m2•hPa•s
def vapor(lat,lon,level,u,v,q):
    # # # 计算单层水汽通量和水汽通量散度
    qv_u = u * q / (constants.g * 10 ** -2)  # g的单位为m/s2，换算为N/kg，再换算为10-2hPa·m2/kg,最终单层水汽通量的单位是kg/m•hPa•s
    qv_v = v * q / (constants.g * 10 ** -2)  # 计算q*v/g,单位是kg/m•hPa•s
    dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat) # 将经纬度转换为格点距离

    div_qv = np.zeros((level.shape[0], lat.shape[0], lon.shape[0]))
    for j in range(level.shape[0]):
        div_qv[j] = mpcalc.divergence(u = qv_u[j],v = qv_v[j],dx = dx ,dy = dy)   # 单位是kg/m2•hPa•s

    b = div_qv
    a = np.sqrt(qv_u * qv_u + qv_v * qv_v)
    return a,b

def nc_create(lat,lon,data,varname,level=None,time=None):
    new_ds = xr.Dataset()

    new_ds['lat'] = xr.DataArray(lat, dims=['lat']).astype('float64')
    new_ds['lon'] = xr.DataArray(lon, dims=['lon']).astype('float64')

    if level is not None:
        new_ds['level'] = xr.DataArray(level, dims=['level']).astype('float64')
        new_ds['level'].attrs['long_name'] = 'Level'
        new_ds['level'].attrs['units'] = 'hPa'  # 你可以根据实际情况设置单位

    if time is not None:
        new_ds['time'] = xr.DataArray(time, dims=['time'])
        # 设置 'time' 维度的属性
        new_ds['time'].attrs['long_name'] = 'Time'
        new_ds['time'].attrs['units'] = 'datetime64[ns]'

    # 设置 'lat' 维度的属性
    new_ds['lat'].attrs['long_name'] = 'Latitude'
    new_ds['lat'].attrs['units'] = 'degrees'
    new_ds['lat'].attrs['description'] = 'Array of latitude values'
    # 设置 'lon' 维度的属性
    new_ds['lon'].attrs['long_name'] = 'Longitude'
    new_ds['lon'].attrs['units'] = 'degrees'
    new_ds['lon'].attrs['description'] = 'Array of Longitude values'


    # 添加变量
    for i, data in enumerate(data_list):
        var_name = varname[i] if i < len(varname) else f'var{i + 1}'
        coords = {'lat': new_ds['lat'], 'lon': new_ds['lon']}
        if 'level' in new_ds and level is not None:
            coords['level'] = new_ds['level']
        if 'time' in new_ds and time is not None:
            coords['time'] = new_ds['time']
        dims = ['time', 'level', 'lat', 'lon'] if 'level' in new_ds and 'time' in new_ds and level is not None and time is not None else \
               ['time', 'lat', 'lon'] if 'time' in new_ds and time is not None else \
               ['level', 'lat', 'lon'] if 'level' in new_ds and level is not None else \
               ['lat', 'lon']
        new_ds[var_name] = xr.DataArray(data, coords=coords, dims=dims).astype('float32')

    return new_ds


def time_change(time):
    iso_format_str_list = []
    # 将输入字符串解析为 datetime 对象
    for time_str in time:
        try:
            dt_object = datetime.strptime(time_str, '%Y%m%d%H')
            #转换为标准时间格式
            iso_format_str = dt_object.isoformat()
            iso_format_str_list.append(iso_format_str)
        except ValueError:
            print('请输入YYYYMMDDHH格式数据')
    return iso_format_str_list


