import math
import numpy as np
import os
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from shapely.geometry import MultiPolygon
from shapely import wkt
import matplotlib.pyplot as plt
import subprocess
import json
import ast


country_buffers_dataframe= pd.read_csv('country_buffer_data.csv')

country_buffers_dataframe['geometry']=country_buffers_dataframe['geometry'].apply(wkt.loads)
# country_buffers_geodataframe = gpd.GeoDataFrame(country_buffers_dataframe, crs='epsg:4326')

#INPUT_PARAMETERS#
country_name = input('Enter Country name :')
catalog_file_name = input('Write name of catalog(should be stored in \"Catalogs\" folder) for selecting regional earthquakes :')
country_dataFrame=pd.read_csv('Catalogs/{file_name}'.format(file_name=catalog_file_name), sep='|')  #here we have to input country catalog file name as input.
longitude_list=country_dataFrame['longitude'].to_list()
latitude_list=country_dataFrame['latitude'].to_list()
boolean_list_for_contained_points=[]
loop_count=0
print('started selecting earthquake for {0}......'.format(country_name))

for i in range(country_buffers_dataframe.shape[0]):

	polygon_id=country_buffers_dataframe.iloc[i]['id']
	#here we used country_name[:-1] because it is a string of form 'name'+'\n'
	if country_name in polygon_id:
		polygon=country_buffers_dataframe.iloc[i]['geometry']
		
		for j in range(country_dataFrame.shape[0]):
			point=Point(longitude_list[j],latitude_list[j])

			if loop_count==0:
				if polygon.contains(point):
					boolean_list_for_contained_points.append(True)
				else:
					boolean_list_for_contained_points.append(False)

			else:
				if polygon.contains(point):
					boolean_list_for_contained_points[j]=True
		loop_count=loop_count+1

print(len(boolean_list_for_contained_points),sum([1 for l in boolean_list_for_contained_points if l==True]))
if boolean_list_for_contained_points!=[]:
	new_country_dataFrame=country_dataFrame.loc[boolean_list_for_contained_points]

	new_country_dataFrame.to_csv('country_buffered_catalogs/final_merge_{0}.csv'.format(country_name),sep='|')
	print('All earthquakes selected for {0} within approx. 100km boundary.'.format(country_name))
	quit()