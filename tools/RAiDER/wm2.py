import datetime
import logging
import os
from abc import ABC, abstractmethod

import h5py
import numpy as np
import netCDF4

from osgeo import osr
from pyproj import CRS, Transformer

_standard_projection = {


def writeNETCDF4(
        time,
        varList,
        varAttribsDict,
        outName=None, 
        NoDataValue=-3.4028234e+38, 
        chunkSize=128
    ):
    '''
    Write a set of georeferenced variables to a netcdf4 file

    Parameters
    ----------

    Returns
    -------

    '''
    nc_outfile = netCDF4.Dataset(
            outName,
            'w',
            clobber=True,
            format='NETCDF4'
        )

    nc_outfile.setncattr('Conventions','CF-1.6')
    nc_outfile.setncattr('datetime',datetime.datetime.strftime(self._time, "%Y_%m_%dT%H_%M_%S"))
    nc_outfile.setncattr('date_created',datetime.datetime.now().strftime("%Y_%m_%dT%H_%M_%S"))
    title='Weather model data and delay calculations'
    nc_outfile.setncattr('title',title)
    # chunking
    dimidY, dimidX, dimidZ = self._t.shape
    chunk_lines_Y = np.min([chunk, dimidY])
    chunk_lines_X = np.min([chunk, dimidX])
    ChunkSize = [dimidZ, chunk_lines_Y, chunk_lines_X]
    
    
    nc_outfile.createDimension('x',dimidX)
    nc_outfile.createDimension('y',dimidY)
    nc_outfile.createDimension('z',dimidZ)
    
    varname='x'
    datatype=np.dtype('float64')
    dimensions=('x')
    FillValue=None
    var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue)
    var.setncattr('standard_name','projection_x_coordinate')
    var.setncattr('description','weather model native x')
    var.setncattr('units','degrees_east')
    var[:] = self._xs.astype(np.float64)
    
    varname='y'
    datatype=np.dtype('float64')
    dimensions=('y')
    FillValue=None
    var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue)
    var.setncattr('standard_name','projection_y_coordinate')
    var.setncattr('description','weather model native y')
    var.setncattr('units','degrees_north')
    var[:] = self._ys.astype(np.float64)
    
    varname='z'
    datatype=np.dtype('float32')
    dimensions=('z')
    FillValue=None
    var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue)
    var.setncattr('standard_name','projection_z_coordinate')
    var.setncattr('description','vertical coordinate')
    var.setncattr('units','m')
    var[:] = self._zs.astype(np.float32)
    
    epsg = 4326
    srs=osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    tran = [self._xs[0], self._xs[1]-self._xs[0], 0.0, self._ys[0], 0.0, self._ys[1]-self._ys[0]]
    
    mapping_name='WGS84'
    grid_mapping='WGS84'  # need to set this as an attribute for the image variables
    datatype=np.dtype('S1')
    dimensions=()
    FillValue=None
    
    var = nc_outfile.createVariable(mapping_name,datatype,dimensions, fill_value=FillValue)
    # variable made, now add attributes
    
    var.setncattr('grid_mapping_name',grid_mapping)
    var.setncattr('straight_vertical_longitude_from_pole',srs.GetProjParm('central_meridian'))
    var.setncattr('false_easting',srs.GetProjParm('false_easting'))
    var.setncattr('false_northing',srs.GetProjParm('false_northing'))
    var.setncattr('latitude_of_projection_origin',np.sign(srs.GetProjParm('latitude_of_origin'))*90.0)
    var.setncattr('latitude_of_origin',srs.GetProjParm('latitude_of_origin'))
    var.setncattr('semi_major_axis',float(srs.GetAttrValue('GEOGCS|SPHEROID',1)))
    var.setncattr('scale_factor_at_projection_origin',1)
    var.setncattr('inverse_flattening',float(srs.GetAttrValue('GEOGCS|SPHEROID',2)))
    var.setncattr('spatial_ref',srs.ExportToWkt())
    var.setncattr('spatial_proj4',srs.ExportToProj4())
    var.setncattr('spatial_epsg',epsg)
    var.setncattr('GeoTransform',' '.join(str(x) for x in tran))  # note this has pixel size in it - set  explicitly above
    
    varname='latitude'
    datatype=np.dtype('float64')
    dimensions=('z','y','x')
    FillValue=NoDataValue
    var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue,zlib=True, complevel=2, shuffle=True, chunksizes=ChunkSize)
    var.setncattr('grid_mapping',mapping_name)
    var.setncattr('description','latitude')
    var.setncattr('standard_name','latitude')
    var.setncattr('units','degrees_north')
    self._lats[np.isnan(self._lats)] = NoDataValue
    var[:] = self._lats.swapaxes(0,2).swapaxes(1,2).astype(np.float64)
    
    varname='longitude'
    datatype=np.dtype('float64')
    dimensions=('z','y','x')
    FillValue=NoDataValue
    var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue,zlib=True, complevel=2, shuffle=True, chunksizes=ChunkSize)
    var.setncattr('grid_mapping',mapping_name)
    var.setncattr('description','longitude')
    var.setncattr('standard_name','longitude')
    var.setncattr('units','degrees_east')
    self._lons[np.isnan(self._lons)] = NoDataValue
    var[:] = self._lons.swapaxes(0,2).swapaxes(1,2).astype(np.float64)
    
    varname='t'
    datatype=np.dtype('float32')
    dimensions=('z','y','x')
    FillValue=NoDataValue
    var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue,zlib=True, complevel=2, shuffle=True, chunksizes=ChunkSize)
    var.setncattr('grid_mapping',mapping_name)
    var.setncattr('description','temperature')
    var.setncattr('standard_name','temperature')
    var.setncattr('units','K')
    self._t[np.isnan(self._t)] = NoDataValue
    var[:] = self._t.swapaxes(0,2).swapaxes(1,2).astype(np.float32)
    
    varname='p'
    datatype=np.dtype('float32')
    dimensions=('z','y','x')
    FillValue=NoDataValue
    var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue,zlib=True, complevel=2, shuffle=True, chunksizes=ChunkSize)
    var.setncattr('grid_mapping',mapping_name)
    var.setncattr('description','pressure')
    var.setncattr('standard_name','pressure')
    var.setncattr('units','Pa')
    self._p[np.isnan(self._p)] = NoDataValue
    var[:] = self._p.swapaxes(0,2).swapaxes(1,2).astype(np.float32)
    
    varname='e'
    datatype=np.dtype('float32')
    dimensions=('z','y','x')
    FillValue=NoDataValue
    var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue,zlib=True, complevel=2, shuffle=True, chunksizes=ChunkSize)
    var.setncattr('grid_mapping',mapping_name)
    var.setncattr('description','humidity')
    var.setncattr('standard_name','humidity')
    var.setncattr('units','Pa')
    self._e[np.isnan(self._e)] = NoDataValue
    var[:] = self._e.swapaxes(0,2).swapaxes(1,2).astype(np.float32)
    
    varname='wet'
    datatype=np.dtype('float32')
    dimensions=('z','y','x')
    FillValue=NoDataValue
    var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue,zlib=True, complevel=2, shuffle=True, chunksizes=ChunkSize)
    var.setncattr('grid_mapping',mapping_name)
    var.setncattr('description','wet_refractivity')
    var.setncattr('standard_name','wet_refractivity')
    self._wet_refractivity[np.isnan(self._wet_refractivity)] = NoDataValue
    var[:] = self._wet_refractivity.swapaxes(0,2).swapaxes(1,2).astype(np.float32)
    
    varname='hydro'
    datatype=np.dtype('float32')
    dimensions=('z','y','x')
    FillValue=NoDataValue
    var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue,zlib=True, complevel=2, shuffle=True, chunksizes=ChunkSize)
    var.setncattr('grid_mapping',mapping_name)
    var.setncattr('description','hydrostatic_refractivity')
    var.setncattr('standard_name','hydrostatic_refractivity')
    self._hydrostatic_refractivity[np.isnan(self._hydrostatic_refractivity)] = NoDataValue
    var[:] = self._hydrostatic_refractivity.swapaxes(0,2).swapaxes(1,2).astype(np.float32)
    
    varname='wet_total'
    datatype=np.dtype('float32')
    dimensions=('z','y','x')
    FillValue=NoDataValue
    var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue,zlib=True, complevel=2, shuffle=True, chunksizes=ChunkSize)
    var.setncattr('grid_mapping',mapping_name)
    var.setncattr('description','total_wet_refractivity')
    var.setncattr('standard_name','total_wet_refractivity')
    self._wet_ztd[np.isnan(self._wet_ztd)] = NoDataValue
    var[:] = self._wet_ztd.swapaxes(0,2).swapaxes(1,2).astype(np.float32)
    
    varname='hydro_total'
    datatype=np.dtype('float32')
    dimensions=('z','y','x')
    FillValue=NoDataValue
    var = nc_outfile.createVariable(varname,datatype,dimensions, fill_value=FillValue,zlib=True, complevel=2, shuffle=True, chunksizes=ChunkSize)
    var.setncattr('grid_mapping',mapping_name)
    var.setncattr('description','total_hydrostatic_refractivity')
    var.setncattr('standard_name','total_hydrostatic_refractivity')
    self._hydrostatic_ztd[np.isnan(self._hydrostatic_ztd)] = NoDataValue
    var[:] = self._hydrostatic_ztd.swapaxes(0,2).swapaxes(1,2).astype(np.float32)
    
    nc_outfile.sync() # flush data to disk
    nc_outfile.close()
