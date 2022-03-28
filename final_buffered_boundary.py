import math
import numpy as np
import os
import sys
import json
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from shapely.geometry import MultiPolygon
import matplotlib.pyplot as plt


# #We can use lines below to make it a interactive function.
# country=input('Write the name of country:-')
# buffer_diatance=float(input('Give the width of buffer zone'))

#1)***************************************
#changing coordinates from cartesian to polar co-ordinates.
#list_ is of form [x,y,z].
def cartesian2polar(list_):
	x=list_[0]
	y=list_[1]
	z=list_[2]
	Re=6371
	R=math.sqrt((x**2)+(y**2)+(z**2))
	height=R-Re
	lat=math.asin(z/R)
	arctan_to_find_lon=math.atan(y/x)
	if y>=0:
		if arctan_to_find_lon<0:
			lon=math.pi-abs(arctan_to_find_lon)
		else:
			lon=arctan_to_find_lon

	else:
		if arctan_to_find_lon<0:
			lon=arctan_to_find_lon
		else:
			lon=-(math.pi-arctan_to_find_lon)

	res=[math.degrees(lon),math.degrees(lat),R,height]
	return(res)


#2)*******************************
#changing coordinates from geographic to cartesian co-ordinates.
#list_cord is of form [lon,lat]
def geo2cartesian(list_cord):
	Re=6371
	lon=list_cord[0]
	lat=list_cord[1]

	x=Re*(math.cos(math.radians(lat)))*(math.cos(math.radians(lon)))
	y=Re*(math.cos(math.radians(lat)))*(math.sin(math.radians(lon)))
	z=Re*(math.sin(math.radians(lat)))

	res=[x,y,z]
	return(res)

#3)********************************
#function which takes in a vector and gives back
def unit_vec(vec_array):
	mod_p=math.sqrt(sum(vec_array**2))
	unit_vec=vec_array/mod_p
	return(unit_vec)


#4)********************************
#function to get vector perpendicular to an arc between 2 points.
def perpendicular(a1,a2,polygon):
	perp=np.cross(a1,a2)
	unit_perp=unit_vec(perp)
	# Finding direction perpendicular to the arc which is directed outside the polygon.
	def ext_perp_dir(point_vector,unit_perp):
		Re=6371
		vector_array=np.array(point_vector)
		newpoint1=cartesian2polar(vector_array+0.000001*unit_perp)
		lat1,lon1=newpoint1[1],newpoint1[0]

		#Converting a touple to a shapely point object.
		point1=Point(lon1,lat1)

		if not polygon.contains(point1):
			return([+1*unit_perp])
		else:
			return([-1*unit_perp])

	ext_a2=ext_perp_dir(a2,unit_perp)

	return(ext_a2[0])


#5)********************************
#Here list_ is a list of 3 [lon,lat] co-ordinates.
def buffer_point_list(list_,buffer,polygon):
	Re=6371
	a1=np.array(geo2cartesian(list_[0]))
	a2=np.array(geo2cartesian(list_[1]))
	a3=np.array(geo2cartesian(list_[2]))

	#finding tangent vectors at a2.
	diff1=a1-a2
	diff2=a3-a2

	unit_a2=unit_vec(a2)

	v1=diff1-np.dot(unit_a2,diff1)*unit_a2
	v2=diff2-np.dot(unit_a2,diff2)*unit_a2

	#defining signs to differentiate between 4 cases. 
	sign_cos=np.sign(np.dot(v1,v2))
	v1_v2=cartesian2polar(v1+v2)
	point_v1_v2=Point(v1_v2[0],v1_v2[1])
	if polygon.contains(point_v1_v2):
		sign_sine=-1
	else:
		sign_sine=1

	#Different Cases.
	if sign_sine==1:
		p1=perpendicular(a1,a2,polygon)
		p2=perpendicular(a2,a3,polygon)
		f=unit_vec(p1+p2)
		alpha=buffer/Re
		f0=a2+math.tan(alpha)*Re*p1
		f1=a2+math.tan(alpha)*Re*f
		f2=a2+math.tan(alpha)*Re*p2
		polar_f0=cartesian2polar(f0)
		polar_f1=cartesian2polar(f1)
		polar_f2=cartesian2polar(f2)
		return([[polar_f0[0],polar_f0[1]],[polar_f1[0],polar_f1[1]],[polar_f2[0],polar_f2[1]]])
		# return([[polar_f1[0],polar_f1[1]]])

	else:
		unit_diff1=unit_vec(diff1)
		unit_diff2=unit_vec(diff2)
		cos_theta=np.dot(unit_diff2,unit_diff1)
		sin_half_theta=math.sqrt((1-cos_theta)/2)
		#here we take minus sign because diff1 and diff2 are directed inside polygon.
		f=-unit_vec(unit_diff2+unit_diff1)
		alpha=buffer/Re
		mod_p1=math.tan(alpha)*Re
		f1=a2+(mod_p1/sin_half_theta)*f
		polar_f1=cartesian2polar(f1)
		return([[polar_f1[0],polar_f1[1]]])


#6)********************************
#function to convert plotly geometry objects to list of points.
#In next line border is a geopandas dataframe obtained using shapefile. 
def process_data(border):
	g = [i for i in border.geometry]
	all_coords={'border':[],'type':[]}
	for i in range(len(g)):

		if 'multipolygon' in str(type(g[i])):
			# print(type(g[i]))
			mpolygon=[]
			for b in g[i].boundary: # for first feature/row
			    coords = np.dstack(b.coords.xy).tolist()
			    mpolygon.append(*coords)
			all_coords['border'].append(mpolygon)
			all_coords['type'].append('MP')
		else:
			# print(type(g[i]))
			try:
				coords = np.dstack(g[i].boundary.coords.xy).tolist()
				all_coords['border'].append(*coords)
				all_coords['type'].append('P')
			except:
				multipart_polygon=[]
				for b in g[i].boundary:
					coords = np.dstack(b.coords.xy).tolist()
					multipart_polygon.append(*coords)
				all_coords['border'].append(multipart_polygon)
				all_coords['type'].append('MPP')

	name_border_data=pd.DataFrame({'NAME':border['NAME'],'ISO3':border['ISO3'], 'geometry':border['geometry'],'usable_geo':all_coords['border'],'type':all_coords['type']})
	return(name_border_data)


#7)********************************
#Function which takes in information about border,country,border type and return buffered border, Polygon in case of 'P' and 'MPP', list of polygons in case of 'MP',
#and type of country border which is any one of 'P', 'MPP', 'MP'.
def border_buffer(dataframe,country_name,ISO3,bufferDistance):
	try:
		dataframe1=dataframe.loc[dataframe['NAME'] == country_name]
		data_list=dataframe1['usable_geo'].tolist()[0]
	except:
		dataframe1=dataframe.loc[dataframe['ISO3'] == ISO3]
		data_list=dataframe1['usable_geo'].tolist()[0]
	# print(data_list)

	#I have put if else statements in for loop because first and last point in polygon is same and their difference will give 0.
	if dataframe1['type'].tolist()[0]=='P':
		#Converting a list of points to a shapely polygon object.
		polygon=Polygon(data_list)
		buffer_border=[]
		for i in range(len(data_list)):
			if i==len(data_list)-1:
				data_points=buffer_point_list([data_list[i-1],data_list[i],data_list[1]],bufferDistance,polygon)
			elif i==0:
				data_points=buffer_point_list([data_list[-2],data_list[i],data_list[1]],bufferDistance,polygon)
			else:
				data_points=buffer_point_list([data_list[i-1],data_list[i],data_list[(i+1)]],bufferDistance,polygon)
			for i in range(len(data_points)):
				if not polygon.contains(Point(data_points[i])):
					buffer_border=buffer_border + [data_points[i]]
		buffer_border=buffer_border[:] + [buffer_border[0]]
		return(buffer_border,polygon,'P')

	#In this case when we have a hollow polygon I am taking all the earthquakes in region inside the polygon to be used for calculating completeness.
	if dataframe1['type'].tolist()[0]=='MPP':
		#Converting a list of points to a shapely polygon object.
		polygon=Polygon(data_list[0])
		# print(polygon)
		buffer_border=[]
		for i in range(len(data_list[0])):
			if i==len(data_list[0])-1:
				data_points=buffer_point_list([data_list[0][i-1],data_list[0][i],data_list[0][1]],bufferDistance,polygon)
			elif i==0:
				data_points=buffer_point_list([data_list[0][-2],data_list[0][i],data_list[0][1]],bufferDistance,polygon)

			else:
				data_points=buffer_point_list([data_list[0][i-1],data_list[0][i],data_list[0][(i+1)]],bufferDistance,polygon)
			for i in range(len(data_points)):
				if not polygon.contains(Point(data_points[i])):
					buffer_border=buffer_border + [data_points[i]]
		buffer_border=buffer_border[:] + [buffer_border[0]]
		return(buffer_border,polygon,'MPP')

	if dataframe1['type'].tolist()[0]=='MP':
		#Converting a list of points to a shapely multipolygon object.
		data_list1=[Polygon(i) for i in data_list]
		polygon1=data_list1
		polygon=MultiPolygon(data_list1)
		multi_buffer_border=[]
		j=0
		for fig in data_list:
			buffer_border=[]
			for i in range(len(fig)):
				if i==0:
					data_points=buffer_point_list([fig[i-2],fig[i],fig[1]],bufferDistance,polygon)
				elif i==len(fig)-1:
					data_points=buffer_point_list([fig[-2],fig[i],fig[1]],bufferDistance,polygon)
				else:
					data_points=buffer_point_list([fig[i-1],fig[i],fig[(i+1)]],bufferDistance,polygon)
				
				for i in range(len(data_points)):
					if not polygon.contains(Point(data_points[i])):
						buffer_border=buffer_border + [data_points[i]]
			j=j+1

			#We have put this condition because in case of multipolygons there can some polygons whose buffer
			#will completely lie inside the countries multipolygon.
			if buffer_border==[]:
				True
			else:
				buffer_border=buffer_border[:] + [buffer_border[0]]

			# If statement in the next line is for the case when the polygon of a multipolygon is present inside another polygon or just outside it near the boundry of it.
			# and only 2 points of the boundary of this inner polygon lies outside the bigger polygon. Here I assume that these 2 boundary points lie inside the boundary of bigger polygon or will be close to it.
			if len(buffer_border)>2:
				multi_buffer_border.append(buffer_border)
		return(multi_buffer_border,polygon1,'MP')

                               #END OF DEFINITIONS FOR PART 1#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#PROGRAM TO REMOVE BUFFER POINTS WHICH ARE AT A DISTANCE < 100Km FROM COUNTRY BORDER.

#FUNCTIONS.

#4)********************************
#function to get vector perpendicular to an arc between 2 points.
def unit_perpendicular(a1,a2):
	perp=np.cross(a1,a2)
	unit_perp=unit_vec(perp)
	return(unit_perp)

#list_ is a list containing 2 points of form [lon,lat].
def perp_vecs(list_):
	Re=6371
	a1=np.array(geo2cartesian(list_[0]))
	a2=np.array(geo2cartesian(list_[1]))

	#finding tangent vectors at a2.
	diff1=a1-a2

	unit_a2=unit_vec(a2)

	v1=diff1-np.dot(unit_a2,diff1)*unit_a2
	unit_v1=unit_vec(v1)
	return([unit_v1,unit_perpendicular(a1,a2)])

#list_perp_vector is a list containing 2 perpendicular vectors at point a2. 
def reg_polygon_inorout(polygon,number_of_sides,list_perp_vector,point,buffer):

	two_pi_radians=2*math.pi
	division_angle=two_pi_radians/number_of_sides
	Re=6371
	list_of_edges=[]
	for i in range(number_of_sides):
		rel_dir_edge=math.sin(i*division_angle)*list_perp_vector[0]+math.cos(i*division_angle)*list_perp_vector[1]
		alpha=buffer/Re
		edge=point+math.tan(alpha)*Re*rel_dir_edge
		polar_f0=cartesian2polar(edge)
		list_of_edges.append([polar_f0[0],polar_f0[1]])

	j=0
	list_of_boundary_points=[]
	for i in list_of_edges:
		point=Point(i)
		if polygon.contains(point):
			return('not allowed')
			j=j+1

	return('allowed')

#function below selects the points on initial boundary which are at a distance of 100Km from the boundary.
def boundary_corrector(polygon,boundary,n=10,buffer=100):
	Re=6371
	new_boundary=[]
	for i in range(len(boundary)-1):

		list_points=[boundary[i],boundary[i+1]]
		point=geo2cartesian(boundary[i+1])
		point_geo=boundary[i+1]

		entry_check=reg_polygon_inorout(polygon,n,perp_vecs(list_points),point,buffer)
		if entry_check=='allowed':
			new_boundary.append(point_geo)

	try:
		new_boundary=new_boundary[:]+[new_boundary[0]]
	except:
		print('Some error occoured check your buffer distance value as all boundary points are inside the buffer.','\nCheck boundary_corrector function.')
	return(new_boundary)
								#end of functions part 2#
#====================================================================================================#
					  #EXTRACTING AND PLOTTING BUFFER DATA USING BORDER DATA#

#---INPUT---
country_name_list = sys.argv[3:]
country_buffered_border_polygon_list=[]
country_polygon_id_list=[]

loop_cycle=0
for name in country_name_list:

	#---INPUT PARAMETERS---
	buffer_distance=float(sys.argv[1])
	country=name
	country_ISO3='not specified'
	number_of_polygon_sides=int(sys.argv[2])

	#Reading and extracting border data.
	border_dataFrame = gpd.read_file('TM_WORLD_BORDERS-0.3.shp')

	# Processing the data is necessary to coordinates of boundary in a usable form.
	processed_data=process_data(border_dataFrame)
	buffered_border_DATA=border_buffer(processed_data,country,country_ISO3,buffer_distance)

	#setting values for variables for different countries.
	polygon_type=buffered_border_DATA[2]
	plot_polygon=buffered_border_DATA[1]
	raw_border=buffered_border_DATA[0]

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
	# Correcting Raw Borders.

	print('Started Correcting Raw Border')

	if polygon_type=='MP':
		correct_buffered_border=[]
		for i in range(len(raw_border)):
			part_correct_buffered_border=boundary_corrector(plot_polygon[i],raw_border[i],n=number_of_polygon_sides,buffer=buffer_distance)
			correct_buffered_border.append(part_correct_buffered_border)
			correct_buffered_border_x=[j[0] for j in part_correct_buffered_border]
			correct_buffered_border_y=[j[1] for j in part_correct_buffered_border]

	else:
		correct_buffered_border=boundary_corrector(plot_polygon,raw_border,n=number_of_polygon_sides,buffer=buffer_distance)
		correct_buffered_border_x=[i[0] for i in correct_buffered_border]
		correct_buffered_border_y=[i[1] for i in correct_buffered_border]


	#===============================================================================================#
		#BUFFERED POLYGON NEEDED TO FIND IF THE EARTHQUAKE LYING IN BUFFERED REGION OF COUNTRY.#

	if polygon_type=='MP':
		for i in range(len(correct_buffered_border)):   #Multi_correct_buffered_border=[Polygon(i) 
			polygon_id=name+str(i)
			buffered_polygon=Polygon(correct_buffered_border[i])
			country_polygon_id_list.append(polygon_id)
			country_buffered_border_polygon_list.append(buffered_polygon)
	else:
		i=0
		polygon_id=name+str(i)
		buffered_polygon=Polygon(correct_buffered_border)
		country_polygon_id_list.append(polygon_id)
		country_buffered_border_polygon_list.append(buffered_polygon)
	print(type(buffered_polygon))

	dictionary={'id':country_polygon_id_list,'geometry':country_buffered_border_polygon_list}
	print(str(dictionary))
	required_border_data=gpd.GeoDataFrame(dictionary,crs="EPSG:4326")

	# #writing a dictionary to a json file.
	# with open('country_buffer_border_data','w') as file:
	# 	json.dump(dictionary,file)

	#Writing buffer border data in a csv file.
	#++++++++++++++IMPORTANT+++++++++++
	# We should not forget about specifying crs="EPSG:4326" while reading this csv file to a geoDataFrame.
	required_border_data.to_csv('country_buffer_data.csv',index=False)
	
	loop_cycle=loop_cycle+1
	print('border generation completed for {count} countries'.format(count=loop_cycle))