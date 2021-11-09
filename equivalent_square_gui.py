from tkinter import Tk as tk_Tk
from tkinter import Frame as tk_Frame
from tkinter import Label as tk_Label
from tkinter import Entry as tk_Entry
from tkinter import StringVar as tk_StringVar
from tkinter import Button as tk_Button
from tkinter import ttk
from tkinter import IntVar as tk_IntVar
from tkinter import Radiobutton as tk_Radiobutton
from tkinter import Checkbutton as tk_Checkbutton
from tkinter.filedialog import asksaveasfile
from tkinter.filedialog import askopenfilename

import tkinter as tk

import equivalent_square_lib as eq
import equivalent_square_lib as es

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib import colors
import matplotlib.backends.backend_tkagg as tkagg
import matplotlib.patches as patches

import numpy as np

import pandas as pd
from pandastable import Table, TableModel
from decimal import Decimal

#global variables
diff_plot_ster = None
diff_plot_geo = None

diff_plot_ster_bar = None
diff_plot_geo_bar = None

mlc_plot = None

table_data = None
pt = None
check_pt = 0


#plotting window class
class Plotwindow():
	def __init__(self, masterframe, size):
		(w,h)=size
		inchsize=(w/25.4, h/25.4)
		self.figure = Figure(inchsize)
		self.axes = self.figure.add_subplot(111)

		# create canvas as matplotlib drawing area
		self.canvas = FigureCanvasTkAgg(self.figure, master=masterframe)
		self.canvas.get_tk_widget().pack()

	def plotxy(self, x, y):
		self.axes.plot(x,y)
		self.canvas.draw()

	def clearplot(self):
		self.axes.cla()
		self.canvas.draw()

	def to_plot(self,d_min=None,d_max=4,increment=0.1,version=0,sheet_name="Sheet 1", c_map='viridis',pre_calc_a=None,z=10,tpr2010=0.775,repres=0,absolute=0,show_limit=False,limit=0.5,cmap='viridis',check_cmin=False,check_cmax=False,cmin_entry=None,cmax_entry=False,es_mode='WFF-WFF',es_method='Axis dose',energy='6mv'):

		if d_min==None:
			d_min=increment

		field_length = np.linspace(d_min,d_max,num=int((d_max-d_min+increment)/increment))
		rect_fields = []

		if pre_calc_a == None:
			progress['value'] = 0
			progress.update()
			pre_calc_a = []
			progress_count = 0
			progress_max = ((len(field_length)**2+len(field_length))/2)
			for i in range(len(field_length)):
				pre_calc_a_part = [None]*i
				for j in range(i,len(field_length)):
					progress_count += 1
					if es_method == 'Axis dose':
						if es_mode == 'WFF-WFF':
							pre_calc_a_part.append(es.EquivalentSquare(field_length[i],field_length[j],z=z,tpr2010=tpr2010))
						elif es_mode == 'FFF-FFF':
							pre_calc_a_part.append(es.EquivalentSquareFFF(field_length[i],field_length[j],z=z,tpr2010=tpr2010,energy=energy))
						else:
							pre_calc_a_part.append(es.EquivalentSquareFFFWFF(field_length[i],field_length[j],z=z,tpr2010=tpr2010,energy=energy))
					else:
						if es_mode == 'WFF-WFF':
							pre_calc_a_part.append(es.EquivalentSquareFFFTPR(field_length[i],field_length[j],z=z,tpr2010=tpr2010))
						elif es_mode == 'FFF-FFF':
							pre_calc_a_part.append(es.EquivalentSquareFFFTPR(field_length[i],field_length[j],z=z,tpr2010=tpr2010,energy=energy))
						else:
							pre_calc_a_part.append(es.EquivalentSquareFFFWFFTPR(field_length[i],field_length[j],z=z,tpr2010=tpr2010,energy=energy))

					progress['value'] = progress_count*100/progress_max
					progress.update()
				pre_calc_a.append(pre_calc_a_part)



		if repres==0:
			if version == 0:
				self.axes.set_title('Differences Geometric Mean')
				column_list = []
				for i in range(len(field_length)):
					column_list.append((field_length[i],'d_geo'))
					rect_fields_part = [np.nan]*i
					for j in range(i,len(field_length)):
						#es_obj = es.EquivalentSquare(field_length[i],field_length[j])
						es_obj = pre_calc_a[i][j]
						field_result = es_obj.geo_dif
						rect_fields_part.append(field_result)
						
					rect_fields.append(rect_fields_part)
			else:
				self.axes.set_title('Differences Sterling')
				column_list = []
				for i in range(len(field_length)):
					column_list.append((field_length[i],'d_sterling'))
					rect_fields_part = [np.nan]*i
					for j in range(i,len(field_length)):
						#es_obj = es.EquivalentSquare(field_length[i],field_length[j])
						es_obj = pre_calc_a[i][j]
						field_result = es_obj.sterling_dif
						rect_fields_part.append(field_result)
						
					rect_fields.append(rect_fields_part)
		else:

			if version == 0:
				self.axes.set_title('Differences Geometric Mean')
				column_list = []
				for i in range(len(field_length)):
					column_list.append((field_length[i],'d_geo'))
					rect_fields_part = [np.nan]*i
					for j in range(i,len(field_length)):
						#es_obj = es.EquivalentSquare(field_length[i],field_length[j])
						es_obj = pre_calc_a[i][j]
						field_result = es_obj.geo_dif_rel*100
						rect_fields_part.append(field_result)
						
					rect_fields.append(rect_fields_part)
			else:
				self.axes.set_title('Differences Sterling')
				column_list = []
				for i in range(len(field_length)):
					column_list.append((field_length[i],'d_sterling'))
					rect_fields_part = [np.nan]*i
					for j in range(i,len(field_length)):
						#es_obj = es.EquivalentSquare(field_length[i],field_length[j])
						es_obj = pre_calc_a[i][j]
						field_result = es_obj.sterling_dif_rel*100
						rect_fields_part.append(field_result)
						
					rect_fields.append(rect_fields_part)

		
		self.axes.set_xlabel('x [cm]')
		self.axes.set_ylabel('y [cm]')

		if absolute:
			rect_fields = np.abs(rect_fields)

		if show_limit:


			b_min = np.nanmin(rect_fields)

			if b_min > limit:
				b_min = limit

			b_max = np.nanmax(rect_fields)

			if b_max < limit:
				b_max = limit

			cmap_two = colors.ListedColormap(['green', 'red'])
			bounds_two = [b_min,limit,b_max]
			norm_two = colors.BoundaryNorm(bounds_two, cmap_two.N)

			pic = self.axes.imshow(rect_fields,extent=[d_min,d_max,d_max,d_min],cmap=cmap_two,norm=norm_two)
		else:
			if check_cmin == False:
				cmin_entry = None
			if check_cmax == False:
				cmax_entry = None

			pic = self.axes.imshow(rect_fields,extent=[d_min,d_max,d_max,d_min],cmap=cmap,vmin=cmin_entry,vmax=cmax_entry)


		cbar = self.figure.colorbar(pic,ax=self.axes)

		if repres==0:
			cbar.set_label('Differences [cm]')
		else:
			cbar.set_label('Differences [%]')

		self.figure.tight_layout()
		self.canvas.draw()

		return pre_calc_a

	def to_plot_ges(self,d_min=None,d_max=4,increment=0.1,version=0,sheet_name="Sheet 1", c_map='viridis',pre_calc_a=None,z=10,tpr2010=0.775,repres=0,absolute=0,show_limit=False,limit=0.5,cmap='viridis',check_cmin=False,check_cmax=False,cmin_entry=None,cmax_entry=False,es_mode='WFF-WFF',es_method='Axis dose'):

		if d_min==None:
			d_min=increment

		field_length = np.linspace(d_min,d_max,num=int((d_max-d_min+increment)/increment))
		rect_fields = []

		
		progress['value'] = 0
		progress.update()
		progress_count = 0
		progress_max = ((len(field_length)**2+len(field_length))/2)
		field_a = []
		field_b = []
		for i in range(len(field_length)):
			field_a_part = [0]*i
			field_b_part = [0]*i
			for j in range(i,len(field_length)):
				progress_count += 1
				if es_method == 'Axis dose':
					if es_mode == 'WFF-WFF':
						equi_a = es.EquivalentSquare(field_length[i],field_length[j],z=z,tpr2010=tpr2010)
					elif es_mode == 'FFF-FFF':
						equi_a = es.EquivalentSquareFFF(field_length[i],field_length[j],z=z,tpr2010=tpr2010)
					else:
						equi_a = es.EquivalentSquareFFFWFF(field_length[i],field_length[j],z=z,tpr2010=tpr2010)
				else:
					if es_mode == 'WFF-WFF':
						equi_a = es.EquivalentSquareFFFTPR(field_length[i],field_length[j],tpr2010=tpr2010)
					elif es_mode == 'FFF-FFF':
						equi_a = es.EquivalentSquareFFFTPR(field_length[i],field_length[j],tpr2010=tpr2010)
					else:
						equi_a = es.EquivalentSquareFFFWFFTPR(field_length[i],field_length[j],tpr2010=tpr2010)

				field_a_part.append(equi_a.sterling_dif*10)
				field_b_part.append(equi_a.geo_dif*10)

				progress['value'] = progress_count*100/progress_max
				progress.update()
			field_a.append(field_a_part)
			field_b.append(field_b_part)



		if repres==0:
			if version == 0:
				self.axes.set_title('Differences Geometric Mean')
				column_list = []
				for i in range(len(field_length)):
					column_list.append((field_length[i],'d_geo'))
					rect_fields_part = [np.nan]*i
					for j in range(i,len(field_length)):
						#es_obj = es.EquivalentSquare(field_length[i],field_length[j])
						es_obj = pre_calc_a[i][j]
						field_result = es_obj.geo_dif
						rect_fields_part.append(field_result)
						
					rect_fields.append(rect_fields_part)
			else:
				self.axes.set_title('Differences Sterling')
				column_list = []
				for i in range(len(field_length)):
					column_list.append((field_length[i],'d_sterling'))
					rect_fields_part = [np.nan]*i
					for j in range(i,len(field_length)):
						#es_obj = es.EquivalentSquare(field_length[i],field_length[j])
						es_obj = pre_calc_a[i][j]
						field_result = es_obj.sterling_dif
						rect_fields_part.append(field_result)
						
					rect_fields.append(rect_fields_part)
		else:

			if version == 0:
				self.axes.set_title('Differences Geometric Mean')
				column_list = []
				for i in range(len(field_length)):
					column_list.append((field_length[i],'d_geo'))
					rect_fields_part = [np.nan]*i
					for j in range(i,len(field_length)):
						#es_obj = es.EquivalentSquare(field_length[i],field_length[j])
						es_obj = pre_calc_a[i][j]
						field_result = es_obj.geo_dif_rel*100
						rect_fields_part.append(field_result)
						
					rect_fields.append(rect_fields_part)
			else:
				self.axes.set_title('Differences Sterling')
				column_list = []
				for i in range(len(field_length)):
					column_list.append((field_length[i],'d_sterling'))
					rect_fields_part = [np.nan]*i
					for j in range(i,len(field_length)):
						#es_obj = es.EquivalentSquare(field_length[i],field_length[j])
						es_obj = pre_calc_a[i][j]
						field_result = es_obj.sterling_dif_rel*100
						rect_fields_part.append(field_result)
						
					rect_fields.append(rect_fields_part)

		
		self.axes.set_xlabel('x [cm]')
		self.axes.set_ylabel('y [cm]')

		if absolute:
			rect_fields = np.abs(rect_fields)

		if show_limit:


			b_min = np.nanmin(rect_fields)

			if b_min > limit:
				b_min = limit

			b_max = np.nanmax(rect_fields)

			if b_max < limit:
				b_max = limit

			cmap_two = colors.ListedColormap(['green', 'red'])
			bounds_two = [b_min,limit,b_max]
			norm_two = colors.BoundaryNorm(bounds_two, cmap_two.N)

			pic = self.axes.imshow(rect_fields,extent=[d_min,d_max,d_max,d_min],cmap=cmap_two,norm=norm_two)
		else:
			if check_cmin == False:
				cmin_entry = None
			if check_cmax == False:
				cmax_entry = None

			pic = self.axes.imshow(rect_fields,extent=[d_min,d_max,d_max,d_min],cmap=cmap,vmin=cmin_entry,vmax=cmax_entry)


		cbar = self.figure.colorbar(pic,ax=self.axes)

		if repres==0:
			cbar.set_label('Differences [cm]')
		else:
			cbar.set_label('Differences [%]')

		self.figure.tight_layout()
		self.canvas.draw()

		return pre_calc_a

	def draw_mlc(self,field):
		self.axes.set_xlim(-200,200)
		self.axes.set_ylim(-200,200)
		for leaf in field:
			rect_l = patches.Rectangle((-200,leaf['y']*10), 200+leaf['xl']*10,leaf['dy']*10, linewidth=1, edgecolor='black', facecolor='blue')
			rect_r = patches.Rectangle((leaf['xr']*10,leaf['y']*10), 200-leaf['xr']*10,leaf['dy']*10, linewidth=1, edgecolor='black', facecolor='blue')
			self.axes.add_patch(rect_l)
			self.axes.add_patch(rect_r)

		self.axes.set_aspect('equal', adjustable='box')
		self.canvas.draw()



#Funktion to calculate Equivalent Square and write equi, geo and sterling to window
def calculate(event=0):
	x_entry = entry_x.get()
	y_entry = entry_y.get()
	z_entry = entry_z.get()
	tpr_entry = entry_tpr.get()
	modus_entry = entry_modus.get()
	profile_entry = entry_profile.get()
	#method_entry = entry_method.get()

	x_entry = float(x_entry.replace(',','.'))
	y_entry = float(y_entry.replace(',','.'))
	z_entry = float(z_entry.replace(',','.'))
	tpr_entry = float(tpr_entry.replace(',','.'))

	if modus_entry == 'WFF-WFF':
		eq_square = eq.EquivalentSquare(x_entry,y_entry,z=z_entry,tpr2010=tpr_entry)
		eq_square_tpr = eq.EquivalentSquareTPR(x_entry,y_entry,tpr2010=tpr_entry)
	elif modus_entry == 'FFF-FFF':
		eq_square = eq.EquivalentSquareFFF(x_entry,y_entry,z=z_entry,tpr2010=tpr_entry,energy=profile_entry)
		eq_square_tpr = eq.EquivalentSquareFFFTPR(x_entry,y_entry,tpr2010=tpr_entry,energy=profile_entry)
	else:
		eq_square = eq.EquivalentSquareFFFWFF(x_entry,y_entry,z=z_entry,tpr2010=tpr_entry,energy=profile_entry)
		eq_square_tpr = eq.EquivalentSquareFFFWFFTPR(x_entry,y_entry,tpr2010=tpr_entry,energy=profile_entry)
	

	equi_str = str(eq_square.equi_sq) + ' cm'
	equi_output.configure(text=equi_str,font='Verdana 10 bold',fg='green')

	equi_tpr_str = str(eq_square_tpr.equi_sq) + ' cm'
	equi_tpr_output.configure(text=equi_tpr_str,font='Verdana 10 bold',fg='green')

	sterling_str = str(eq_square.sterling) + ' cm'
	sterling_output.configure(text=sterling_str,font='Verdana 9 bold',fg='orange')

	geo_str = str(eq_square.geometric_mean) + ' cm'
	geo_output.configure(text=geo_str,font='Verdana 9 bold',fg='orange')


def calculate_i(event=0):
	# plot_size = (100,100)
	path_entry = entry_path.get()
	tpr_entry_i = entry_tpr_i.get()
	modus_entry_i = entry_modus_i.get()
	profile_entry_i = entry_profile_i.get()
	select_beam_entry = entry_select_beam.get()
	#method_entry = entry_method.get()

	tpr_entry_i = float(tpr_entry_i.replace(',','.'))

	if modus_entry_i == 'WFF-WFF':
		eq_square = eq.EquivalentSquareIrr(path_entry,tpr2010=tpr_entry_i,beam=select_beam_entry)
		eq_square_tpr = eq.EquivalentSquareTPRIrr(path_entry,tpr2010=tpr_entry_i,beam=select_beam_entry)
	elif modus_entry_i == 'FFF-FFF':
		eq_square = eq.EquivalentSquareFFFIrr(path_entry,tpr2010=tpr_entry_i,energy=profile_entry_i,beam=select_beam_entry)
		eq_square_tpr = eq.EquivalentSquareFFFTPRIrr(path_entry,tpr2010=tpr_entry_i,energy=profile_entry_i,beam=select_beam_entry)
	else:
		eq_square = eq.EquivalentSquareFFFWFFIrr(path_entry,tpr2010=tpr_entry_i,energy=profile_entry_i,beam=select_beam_entry)
		eq_square_tpr = eq.EquivalentSquareFFFWFFTPRIrr(path_entry,tpr2010=tpr_entry_i,energy=profile_entry_i,beam=select_beam_entry)

	# global mlc_plot
	# if mlc_plot != None:
	# 	mlc_plot.canvas.get_tk_widget().pack_forget()
	# mlc_plot = Plotwindow(frame_mlc_plot,plot_size)
	# mlc_plot.draw_mlc(eq_square.field)

	equi_str = str(eq_square.equi_sq) + ' cm'
	equi_output_i.configure(text=equi_str,font='Verdana 10 bold',fg='green')

	equi_tpr_str = str(eq_square_tpr.equi_sq) + ' cm'
	equi_tpr_output_i.configure(text=equi_tpr_str,font='Verdana 10 bold',fg='green')

	sterling_str = str(eq_square.sterling) + ' cm'
	sterling_output_i.configure(text=sterling_str,font='Verdana 9 bold',fg='orange')

	geo_str = str(eq_square.geometric_mean) + ' cm'
	geo_output_i.configure(text=geo_str,font='Verdana 9 bold',fg='orange')


def show_diff(event=0):
	plot_size = (120,120)

	min_entry = entry_min.get()
	max_entry = entry_max.get()
	inc_entry = entry_inc.get()
	z_d_entry = entry_z_d.get()
	tpr_d_entry = entry_tpr_d.get()
	limit_entry = entry_limit.get()
	cmin_entry = entry_cmin.get()
	cmax_entry = entry_cmax.get()


	min_entry = float(min_entry.replace(',','.'))
	max_entry = float(max_entry.replace(',','.'))
	inc_entry = float(inc_entry.replace(',','.'))
	z_d_entry = float(z_d_entry.replace(',','.'))
	tpr_d_entry = float(tpr_d_entry.replace(',','.'))
	limit_entry = float(limit_entry.replace(',','.'))

	modus_d_entry = entry_modus_d.get()
	method_d_entry = entry_method_d.get()
	profile_d_entry = entry_profile_d.get()

	if len(cmin_entry) == 0:
		cmin_entry = None
	else:
		cmin_entry = float(cmin_entry.replace(',','.'))
	if len(cmax_entry) == 0:
		cmax_entry = None
	else:
		cmax_entry = float(cmax_entry.replace(',','.'))


	global diff_plot_ster
	global diff_plot_ster_bar
	if diff_plot_ster != None:
		diff_plot_ster.canvas.get_tk_widget().pack_forget()
	diff_plot_ster = Plotwindow(frame_dif_ster,plot_size)
	pre_calc = diff_plot_ster.to_plot(version=1,d_min=min_entry,d_max=max_entry,increment=inc_entry,z=z_d_entry,tpr2010=tpr_d_entry,repres=repres_rel.get(),absolute=repres_abs.get(),show_limit=check_limit_var.get(),limit=limit_entry,cmap=select_cmap.get(),check_cmin=check_cmin_var.get(),cmin_entry=cmin_entry,check_cmax=check_cmax_var.get(),cmax_entry=cmax_entry,es_mode=modus_d_entry,es_method=method_d_entry,energy=profile_d_entry)
	if diff_plot_ster_bar != None:
		diff_plot_ster_bar.destroy()
	diff_plot_ster_bar = tkagg.NavigationToolbar2Tk(diff_plot_ster.canvas, frame_dif_ster)

	global diff_plot_geo
	global diff_plot_geo_bar
	if diff_plot_geo != None:
		diff_plot_geo.canvas.get_tk_widget().pack_forget()
	diff_plot_geo = Plotwindow(frame_dif_geo,plot_size)
	diff_plot_geo.to_plot(pre_calc_a=pre_calc,version=0,d_min=min_entry,d_max=max_entry,increment=inc_entry,z=z_d_entry,tpr2010=tpr_d_entry,repres=repres_rel.get(),absolute=repres_abs.get(),show_limit=check_limit_var.get(),limit=limit_entry,cmap=select_cmap.get(),check_cmin=check_cmin_var.get(),cmin_entry=cmin_entry,check_cmax=check_cmax_var.get(),cmax_entry=cmax_entry,es_mode=modus_d_entry,es_method=method_d_entry,energy=profile_d_entry)
	if diff_plot_geo_bar != None:
		diff_plot_geo_bar.destroy()
	diff_plot_geo_bar = tkagg.NavigationToolbar2Tk(diff_plot_geo.canvas, frame_dif_geo)

def change_status_more_opt():
	if frame_more_opt.winfo_manager():
		frame_more_opt.grid_forget()
	else:
		frame_more_opt.grid(row=2,pady=10)

def check_tab(event):
	if parent_tab.index(parent_tab.select()) == 0:
		frame_differences.pack_forget()
		frame_tables.pack_forget()
		frame_master_i.pack_forget()
		frame_mlc_plot.pack_forget()
	elif parent_tab.index(parent_tab.select()) == 1:
		frame_tables.pack_forget()
		frame_master_i.pack_forget()
		frame_mlc_plot.pack_forget()
		frame_differences.pack(side='top', padx='20', pady='5')
	elif parent_tab.index(parent_tab.select()) == 2:
		frame_differences.pack_forget()
		frame_master_i.pack_forget()
		frame_mlc_plot.pack_forget()
		frame_tables.pack(side='top', padx='20', pady='5')
	else:
		frame_tables.pack_forget()
		frame_differences.pack_forget()
		frame_master_i.pack(side='left', padx='20', pady='5')
		frame_mlc_plot.pack(side='left')

def save_excel():
	if table_data != None:
		files = [('Excel Document', '*.xlsx')]
		file = asksaveasfile(filetypes = files, defaultextension = files)
		xlsx = pd.ExcelWriter(file.name)
		table_data.to_excel(xlsx,'Equivalent Square',index=True,header=True)
		xlsx.save()

def change(x):
	if x=='orange':
		x='red'
	else:
		x='orange'
	return x

def color_test(data):
	num_values = int(np.shape(data)[1]/np.shape(data)[0])
	a = np.full(data.shape, '', dtype='<U24')
	for i in range(len(a)):
		color = 'red'
		for j in range(len(a[i])):
			if j%num_values == 0:
				color = change(color)
			a[i][j] = 'background-color: ' + color
			if i==int(j/num_values):
				a[i][j] = 'background-color: green'
			if i>int(j/num_values):
				a[i][j] = ''


	return pd.DataFrame(a, index=data.index, columns=data.columns)
		

def show_table(event=0):
	
	
	d_min = float(entry_min_t.get())
	d_max = float(entry_max_t.get())
	inc = float(entry_inc_t.get())

	

	z = float(entry_z_d_t.get())
	tpr2010 = float(entry_tpr_d_t.get())
	modus = entry_modus_d_t.get()
	method_d_t_entry = entry_method_d_t.get()
	profile_d_t_entry = entry_profile_d_t.get()

	num = int(round((d_max - d_min+inc)/inc))

	field_length = np.arange(d_min,d_max+inc,inc)

	pre_calc_a = []
	progress_tab['value'] = 0
	progress_tab.update()
	progress_tab_count = 0
	progress_tab_max = ((len(field_length)**2+len(field_length))/2)
	for i in range(len(field_length)):
		pre_calc_a_part = [None]*i
		for j in range(i,len(field_length)):
			progress_tab_count += 1
			if method_d_t_entry == 'Axis dose':
				if modus == 'WFF-WFF':
					pre_calc_a_part.append(es.EquivalentSquare(field_length[i],field_length[j],z=z,tpr2010=tpr2010))
				elif modus == 'FFF-FFF':
					pre_calc_a_part.append(es.EquivalentSquareFFF(field_length[i],field_length[j],z=z,tpr2010=tpr2010,energy=profile_d_t_entry))
				else:
					pre_calc_a_part.append(es.EquivalentSquareFFFWFF(field_length[i],field_length[j],z=z,tpr2010=tpr2010,energy=profile_d_t_entry))
					
			else:
				if modus == 'WFF-WFF':
					pre_calc_a_part.append(es.EquivalentSquareTPR(field_length[i],field_length[j],tpr2010=tpr2010))
				elif modus == 'FFF-FFF':
					pre_calc_a_part.append(es.EquivalentSquareFFFTPR(field_length[i],field_length[j],tpr2010=tpr2010,energy=profile_d_t_entry))
				else:
					pre_calc_a_part.append(es.EquivalentSquareFFFWFFTPR(field_length[i],field_length[j],tpr2010=tpr2010,energy=profile_d_t_entry))

			progress_tab['value'] = progress_tab_count*100/progress_tab_max
			progress_tab.update()
		pre_calc_a.append(pre_calc_a_part)

	#Only equivalent square kernel Model
	
	column_list = []
	rect_fields = []
	square_fields =[]
	not_square_fields = []
	for i in range(len(field_length)):
		column_list.append(round(field_length[i],1))
		rect_fields_part = [None]*i
		for j in range(i,len(field_length)):
			es_obj = pre_calc_a[i][j]
			field_result = es_obj.equi_sq
			rect_fields_part.append(field_result)
			
			
		rect_fields.append(rect_fields_part)

	df = pd.DataFrame(rect_fields,columns=column_list,index=column_list)
	styled = df.style.apply(color_test,axis=None)

	global table_data
	table_data = styled

	global pt
	global check_pt

	if False:

		if pt != None:
			pt.updateModel(TableModel(df))
		else:
			pt = Table(frame_table_output, dataframe=df,showtoolbar=False, showstatusbar=False)


		for i in range(len(field_length)):
			pt.setRowColors(cols=[i],clr='#bfff00',rows=[i])
			for j in range(0,i):
				if i%2 == 0:
					pt.setRowColors(cols=[i],clr='#ffbf00',rows=[j])
				else:
					pt.setRowColors(cols=[i],clr='#ffff00',rows=[j])

		
		pt.showIndex()
		if check_pt == 0:
			pt.show()
			check_pt = 1
		else:
			pt.show()

	else:

		if pt == None:
			pt = Table(frame_table_output,showtoolbar=False, showstatusbar=False)
			pt.showIndex()
			pt.show()

		model = TableModel(dataframe=df)

		pt.updateModel(model)

def select_dicom():
	filetypes = (('DICOM', '*.dcm'),('All files','*.*'))
	filename = askopenfilename(title = 'Open DICOM file', initialdir = './', filetypes=filetypes)
	#print(filename)
	dicom_path.set(filename)

	fields = es.read_mlc(filename)
	keys = list(fields.keys())

	key_selected = keys[0]


	entry_select_beam.config(values=keys)

	entry_select_beam.set(key_selected)


	field = fields[key_selected][0]

	plot_size = (100,100)
	global mlc_plot
	if mlc_plot != None:
		mlc_plot.canvas.get_tk_widget().pack_forget()
	mlc_plot = Plotwindow(frame_mlc_plot,plot_size)
	mlc_plot.draw_mlc(field)

def update_mlc_plot(event=0):
	filename = entry_path.get()


	fields = es.read_mlc(filename)

	key_selected = entry_select_beam.get()

	field = fields[key_selected][0]

	plot_size = (100,100)
	global mlc_plot
	if mlc_plot != None:
		mlc_plot.canvas.get_tk_widget().pack_forget()
	mlc_plot = Plotwindow(frame_mlc_plot,plot_size)
	mlc_plot.draw_mlc(field)

	



#open tk root
root = tk_Tk()
root.title('Equivalent Square Calculator') #add title to window


#add tabs menu
parent_tab = ttk.Notebook(root)

#add tabs
tab_calculator = tk.Frame(parent_tab)

tab_irreg = tk.Frame(parent_tab)

tab_differences = tk.Frame(parent_tab)

tab_tables = tk.Frame(parent_tab)


parent_tab.bind("<<NotebookTabChanged>>", check_tab)

parent_tab.add(tab_calculator,text='Calculator')
parent_tab.add(tab_differences,text='Differences Plots')
parent_tab.add(tab_tables,text='Create Tables')
parent_tab.add(tab_irreg,text='Irregular Fields')


parent_tab.pack(expand=1,fill='both')

### TAB: CALCULATOR ###

#add and pack master frame
frame_master = tk_Frame(master=tab_calculator)
frame_master.pack(side='top', padx='20', pady='5',fill='both')


#add header and text to masterframe
header1 = tk_Label(frame_master,text = 'Usage:',font = "Verdana 12 bold",anchor='w')
header1.pack(side='top',fill='both')
text1 = tk_Label(frame_master,text = 'Input of x and y dimensions of a rectangular field in cm, the depth in cm and the quality factor TPR2010.', font = 'Verdana 9',anchor='w')
text1.pack(side='top',fill='both')
text2 = tk_Label(frame_master,text = 'As mode equivalent square calculation between WFF and WFF, FFF and FFF, FFF and WFF can be selected.', font = 'Verdana 9',anchor='w')
text2.pack(side='top',fill='both')


##INPUT##

#add and pack frame_x_y to master frame: using for x, y, z and tpr2010 inputs
frame_x_y = tk_Frame(master=frame_master)
frame_x_y.pack(side='top',fill='both',pady=10)


#x input
label_x = tk_Label(frame_x_y, text='x:',anchor='w',font='Verdana 10 bold')
label_x.grid(row=0, column=0,padx=10)

entry_x = tk_Entry(frame_x_y)
entry_x.bind("<Return>",calculate)
entry_x.grid(row=0, column=1)

unit_x = tk_Label(frame_x_y, text='cm',anchor='w',font='Verdana 10')
unit_x.grid(row=0, column=2)


#y input
label_y = tk_Label(frame_x_y, text='y:',anchor='w',font='Verdana 10 bold')
label_y.grid(row=0, column=3,padx=10)

entry_y = tk_Entry(frame_x_y)
entry_y.bind("<Return>",calculate)
entry_y.grid(row=0, column=4)

unit_y = tk_Label(frame_x_y, text='cm',anchor='w',font='Verdana 10')
unit_y.grid(row=0, column=5)


#z input
label_z = tk_Label(frame_x_y, text='depth:',anchor='w',font='Verdana 10 bold')
label_z.grid(row=1, column=0,pady=5,padx=10)

default_z = tk_StringVar(root, value='10')

entry_z = tk_Entry(frame_x_y,textvariable=default_z)
entry_z.bind("<Return>",calculate)
entry_z.grid(row=1, column=1)

unit_z = tk_Label(frame_x_y, text='cm',anchor='w',font='Verdana 10')
unit_z.grid(row=1, column=2)


#tpr2010 input
label_tpr = tk_Label(frame_x_y, text='TPR2010:',anchor='w',font='Verdana 10 bold')
label_tpr.grid(row=1, column=3,padx=10)

default_tpr = tk_StringVar(root, value='0.671')

entry_tpr = tk_Entry(frame_x_y,textvariable=default_tpr)
entry_tpr.bind("<Return>",calculate)
entry_tpr.grid(row=1, column=4)


#Mode input
label_modus = tk_Label(frame_x_y, text='mode:',anchor='w',font='Verdana 10 bold')
label_modus.grid(row=2, column=0,padx=10)

default_modus = tk_StringVar(root, value='WFF-WFF')

entry_modus = ttk.Combobox(frame_x_y,values=['WFF-WFF','FFF-FFF','FFF-WFF'],textvariable=default_modus)
entry_modus.bind("<Return>",calculate)
entry_modus.grid(row=2, column=1)

#FFF Profile input
label_profile = tk_Label(frame_x_y, text='FFF profile*:',anchor='w',font='Verdana 10 bold')
label_profile.grid(row=2, column=3,padx=10)

default_profile = tk_StringVar(root, value='6mv')

entry_profile = ttk.Combobox(frame_x_y,values=['6mv','10mv'],textvariable=default_profile)
entry_profile.bind("<Return>",calculate)
entry_profile.grid(row=2, column=4)

frame_note = tk_Frame(master=frame_master)
frame_note.pack(side='top',fill='both',pady=10)

note_1 = tk_Label(frame_note,text='*Only relevant for modes including FFF-fields.',anchor='w',font='Verdana 9')
note_1.pack(side='top',fill='both')
note_2 = tk_Label(frame_note,text='Selection of the implemented FFF-profiles, used as weighting function in the calculation.',anchor='w',font='Verdana 9')
note_2.pack(side='top',fill='both')
note_3 = tk_Label(frame_note,text='Obtained for ELEKTA Versa HD LINACS at a depth of 1cm.',anchor='w',font='Verdana 9')
note_3.pack(side='top',fill='both')


##OUTPUT##

#title for equivalent square output
equi_title = tk_Label(frame_master,text='Equivalent Square:',anchor='w',font='Verdana 12 bold')
equi_title.pack(side='top',fill='both',pady=10)

#add and pack frame frame_equi in master frame: using for output equi, geo and sterling
frame_equi = tk_Frame(master = frame_master)
frame_equi.pack(side='top',fill='both',pady=0)


#equi output
equi_label = tk_Label(frame_equi,text='Equal axis dose definition:',anchor='w',font='Verdana 10 bold')
equi_label.grid(row=0,column=0,padx=10,sticky='w')

equi_output = tk_Label(frame_equi,anchor='w',font='Verdana 10 bold')
equi_output.grid(row=0,column=1)

#equi TPR2010 output
equi_tpr_label = tk_Label(frame_equi,text='Equal TPR2010 definition:',anchor='w',font='Verdana 10 bold')
equi_tpr_label.grid(row=0,column=2,padx=10,sticky='w')

equi_tpr_output = tk_Label(frame_equi,anchor='w',font='Verdana 10 bold')
equi_tpr_output.grid(row=0,column=3)



#sterling output
sterling_label = tk_Label(frame_equi,text='Sterling equation:',anchor='w',font='Verdana 9 bold')
sterling_label.grid(row=1,column=0,padx=10,sticky='w',pady=5)

sterling_output = tk_Label(frame_equi,anchor='w',font='Verdana 9 bold')
sterling_output.grid(row=1,column=1)


#geo output
geo_label = tk_Label(frame_equi,text='Geometric mean:',anchor='w',font='Verdana 9 bold')
geo_label.grid(row=1,column=2,padx=10,sticky='w')

geo_output = tk_Label(frame_equi,anchor='w',font='Verdana 9 bold')
geo_output.grid(row=1,column=3)


#add calculation button to start calculation of equivalent square
calc_button = tk_Button(frame_master, text="Calculate", fg="black",command=calculate,font='Verdana 10 bold')
calc_button.pack(side='top',pady=10)





### TAB: IRREGULAR FIELDS ###

#add and pack master frame
# frame_irreg = tk_Frame(master=tab_irreg)
# frame_irreg.pack(side='top')
frame_master_i = tk_Frame(master=tab_irreg)
frame_master_i.pack(side='left', padx='20', pady='5')
#frame_master_i.pack(side='left', padx='20', pady='5',fill='both')


#add header and text to masterframe
header1_i = tk_Label(frame_master_i,text = 'Usage:',font = "Verdana 12 bold",anchor='w')
header1_i.pack(side='top',fill='both')
text1_i = tk_Label(frame_master_i,text = 'Input of x and y dimensions of a rectangular field in cm, the depth in cm and the quality factor TPR2010.', font = 'Verdana 9',anchor='w')
text1_i.pack(side='top',fill='both')
text2_i = tk_Label(frame_master_i,text = 'As mode equivalent square calculation between WFF and WFF, FFF and FFF, FFF and WFF can be selected.', font = 'Verdana 9',anchor='w')
text2_i.pack(side='top',fill='both')


##INPUT##

#add and pack frame_input_i to master frame: using for x, y, z and tpr2010 inputs
frame_input_i = tk_Frame(master=frame_master_i)
frame_input_i.pack(side='top',fill='both',pady=10)


#path input
label_path = tk_Label(frame_input_i, text='DICOM-path:',anchor='w',font='Verdana 10 bold')
label_path.grid(row=0, column=0,padx=10)

dicom_path = tk_StringVar(root, value='')

entry_path = tk_Entry(frame_input_i, textvariable=dicom_path)
entry_path.bind("<Return>",calculate_i)
entry_path.grid(row=0, column=1)

button_path = ttk.Button(frame_input_i,text='browse',command=select_dicom)
button_path.grid(row=0, column=2)


#Beam input
label_select_beam = tk_Label(frame_input_i, text='Beam:',anchor='w',font='Verdana 10 bold')
label_select_beam.grid(row=0, column=3,padx=10)

default_select_beam = tk_StringVar(root, value='')

entry_select_beam = ttk.Combobox(frame_input_i,values=[],textvariable=default_select_beam)
entry_select_beam.bind("<Return>",calculate_i)
entry_select_beam.bind('<<ComboboxSelected>>', update_mlc_plot)
entry_select_beam.grid(row=0, column=4)

# unit_x = tk_Label(frame_input_i, text='cm',anchor='w',font='Verdana 10')
# unit_x.grid(row=0, column=2)


# #y input
# label_y = tk_Label(frame_input_i, text='y:',anchor='w',font='Verdana 10 bold')
# label_y.grid(row=0, column=3,padx=10)

# entry_y = tk_Entry(frame_input_i)
# entry_y.bind("<Return>",calculate)
# entry_y.grid(row=0, column=4)

# unit_y = tk_Label(frame_input_i, text='cm',anchor='w',font='Verdana 10')
# unit_y.grid(row=0, column=5)


# #z input
# label_z = tk_Label(frame_input_i, text='depth:',anchor='w',font='Verdana 10 bold')
# label_z.grid(row=1, column=0,pady=5,padx=10)

# default_z = tk_StringVar(root, value='10')

# entry_z = tk_Entry(frame_input_i,textvariable=default_z)
# entry_z.bind("<Return>",calculate)
# entry_z.grid(row=1, column=1)

# unit_z = tk_Label(frame_input_i, text='cm',anchor='w',font='Verdana 10')
# unit_z.grid(row=1, column=2)


#tpr2010 input
label_tpr_i = tk_Label(frame_input_i, text='TPR2010:',anchor='w',font='Verdana 10 bold')
label_tpr_i.grid(row=1, column=0,padx=10)

default_tpr_i = tk_StringVar(root, value='0.671')

entry_tpr_i = tk_Entry(frame_input_i,textvariable=default_tpr)
entry_tpr_i.bind("<Return>",calculate_i)
entry_tpr_i.grid(row=1, column=1)


#Mode input
label_modus_i = tk_Label(frame_input_i, text='mode:',anchor='w',font='Verdana 10 bold')
label_modus_i.grid(row=1, column=3,padx=10)

default_modus_i = tk_StringVar(root, value='WFF-WFF')

entry_modus_i = ttk.Combobox(frame_input_i,values=['WFF-WFF','FFF-FFF','FFF-WFF'],textvariable=default_modus_i)
entry_modus_i.bind("<Return>",calculate_i)
entry_modus_i.grid(row=1, column=4)

#FFF Profile input
label_profile_i = tk_Label(frame_input_i, text='FFF profile*:',anchor='w',font='Verdana 10 bold')
label_profile_i.grid(row=2, column=0,padx=10)

default_profile_i = tk_StringVar(root, value='6mv')

entry_profile_i = ttk.Combobox(frame_input_i,values=['6mv','10mv'],textvariable=default_profile_i)
entry_profile_i.bind("<Return>",calculate_i)
entry_profile_i.grid(row=2, column=1)

frame_note_i = tk_Frame(master=frame_master_i)
frame_note_i.pack(side='top',fill='both',pady=10)

note_1_i = tk_Label(frame_note_i,text='*Only relevant for modes including FFF-fields.',anchor='w',font='Verdana 9')
note_1_i.pack(side='top',fill='both')
note_2_i = tk_Label(frame_note_i,text='Selection of the implemented FFF-profiles, used as weighting function in the calculation.',anchor='w',font='Verdana 9')
note_2_i.pack(side='top',fill='both')
note_3_i = tk_Label(frame_note_i,text='Obtained for ELEKTA Versa HD LINACS at a depth of 1cm.',anchor='w',font='Verdana 9')
note_3_i.pack(side='top',fill='both')


##OUTPUT##

#title for equivalent square output
equi_title_i = tk_Label(frame_master_i,text='Equivalent Square:',anchor='w',font='Verdana 12 bold')
equi_title_i.pack(side='top',fill='both',pady=10)

#add and pack frame frame_equi in master frame: using for output equi, geo and sterling
frame_equi_i = tk_Frame(master = frame_master_i)
frame_equi_i.pack(side='top',fill='both',pady=0)


#equi output
equi_label_i = tk_Label(frame_equi_i,text='Equal axis dose definition:',anchor='w',font='Verdana 10 bold')
equi_label_i.grid(row=0,column=0,padx=10,sticky='w')

equi_output_i = tk_Label(frame_equi_i,anchor='w',font='Verdana 10 bold')
equi_output_i.grid(row=0,column=1)

#equi TPR2010 output
equi_tpr_label_i = tk_Label(frame_equi_i,text='Equal TPR2010 definition:',anchor='w',font='Verdana 10 bold')
equi_tpr_label_i.grid(row=0,column=2,padx=10,sticky='w')

equi_tpr_output_i = tk_Label(frame_equi_i,anchor='w',font='Verdana 10 bold')
equi_tpr_output_i.grid(row=0,column=3)



#sterling output
sterling_label_i = tk_Label(frame_equi_i,text='Sterling equation:',anchor='w',font='Verdana 9 bold')
sterling_label_i.grid(row=1,column=0,padx=10,sticky='w',pady=5)

sterling_output_i = tk_Label(frame_equi_i,anchor='w',font='Verdana 9 bold')
sterling_output_i.grid(row=1,column=1)


#geo output
geo_label_i = tk_Label(frame_equi_i,text='Geometric mean:',anchor='w',font='Verdana 9 bold')
geo_label_i.grid(row=1,column=2,padx=10,sticky='w')

geo_output_i = tk_Label(frame_equi_i,anchor='w',font='Verdana 9 bold')
geo_output_i.grid(row=1,column=3)


#add calculation button to start calculation of equivalent square
calc_button_i = tk_Button(frame_master_i, text="Calculate", fg="black",command=calculate_i,font='Verdana 10 bold')
calc_button_i.pack(side='top',pady=10)

frame_mlc_plot = tk_Frame(master=tab_irreg)
frame_mlc_plot.pack(side='left')




### TAB: DIFFERENCES ###

##INPUT##

#BASIC#

#add and pack differences frame
frame_differences = tk_Frame(master=tab_differences)
frame_differences.pack(side='top', padx='20', pady='5')

#add and pack frame_x_y to master frame: using for x, y, z and tpr2010 inputs
frame_input = tk_Frame(master=frame_differences)
frame_input.grid(row=0,pady=10)

#min input
label_min = tk_Label(frame_input, text='Smin:',anchor='w',font='Verdana 10 bold')
label_min.grid(row=0, column=0,padx=10)

default_min = tk_StringVar(root, value='0.1')

entry_min = tk_Entry(frame_input,textvariable=default_min)
entry_min.bind("<Return>",show_diff)
entry_min.grid(row=0, column=1)

unit_min = tk_Label(frame_input, text='cm',anchor='w',font='Verdana 10')
unit_min.grid(row=0, column=2)

#max input
label_max = tk_Label(frame_input, text='Smax:',anchor='w',font='Verdana 10 bold')
label_max.grid(row=0, column=3,padx=10)

default_max = tk_StringVar(root, value='4')

entry_max = tk_Entry(frame_input,textvariable=default_max)
entry_max.bind("<Return>",show_diff)
entry_max.grid(row=0, column=4)

unit_max = tk_Label(frame_input, text='cm',anchor='w',font='Verdana 10')
unit_max.grid(row=0, column=5)


#inc input
label_inc = tk_Label(frame_input, text='Increment:',anchor='w',font='Verdana 10 bold')
label_inc.grid(row=0, column=6,padx=10)

default_inc = tk_StringVar(root, value='0.1')

entry_inc = tk_Entry(frame_input,textvariable=default_inc)
entry_inc.bind("<Return>",show_diff)
entry_inc.grid(row=0, column=7)

unit_inc = tk_Label(frame_input, text='cm',anchor='w',font='Verdana 10')
unit_inc.grid(row=0, column=8)


#z input
label_z_d = tk_Label(frame_input, text='depth:',anchor='w',font='Verdana 10 bold')
label_z_d.grid(row=1, column=0,pady=5,padx=10)

default_z_d = tk_StringVar(root, value='10')

entry_z_d = tk_Entry(frame_input,textvariable=default_z_d)
entry_z_d.bind("<Return>",show_diff)
entry_z_d.grid(row=1, column=1)

unit_z_d = tk_Label(frame_input, text='cm',anchor='w',font='Verdana 10')
unit_z_d.grid(row=1, column=2)


#tpr2010 input
label_tpr_d = tk_Label(frame_input, text='TPR2010:',anchor='w',font='Verdana 10 bold')
label_tpr_d.grid(row=1, column=3,padx=10)

default_tpr_d = tk_StringVar(root, value='0.671')

entry_tpr_d = tk_Entry(frame_input,textvariable=default_tpr_d)
entry_tpr_d.bind("<Return>",show_diff)
entry_tpr_d.grid(row=1, column=4)

#Mode input
label_modus_d = tk_Label(frame_input, text='mode:',anchor='w',font='Verdana 10 bold')
label_modus_d.grid(row=1, column=6,padx=10)

default_modus_d = tk_StringVar(root, value='WFF-WFF')

entry_modus_d = ttk.Combobox(frame_input,values=['WFF-WFF','FFF-FFF','FFF-WFF'],textvariable=default_modus_d)
entry_modus_d.bind("<Return>",show_diff)
entry_modus_d.grid(row=1, column=7)

#Method input
label_method_d = tk_Label(frame_input, text='method:',anchor='w',font='Verdana 10 bold')
label_method_d.grid(row=2, column=0,padx=10)

default_method_d = tk_StringVar(root, value='Axis dose')

entry_method_d = ttk.Combobox(frame_input,values=['Axis dose','TPR2010'],textvariable=default_method_d)
entry_method_d.bind("<Return>",show_diff)
entry_method_d.grid(row=2, column=1)

#FFF Profile input
label_profile_d = tk_Label(frame_input, text='FFF profile*:',anchor='w',font='Verdana 10 bold')
label_profile_d.grid(row=2, column=3,padx=10)

default_profile_d = tk_StringVar(root, value='6mv')

entry_profile_d = ttk.Combobox(frame_input,values=['6mv','10mv'],textvariable=default_profile_d)
entry_profile_d.bind("<Return>",show_diff)
entry_profile_d.grid(row=2, column=4)



#MORE OPTIONS#

more_opt_button = tk_Button(frame_differences, text="More Options", fg="black",command=change_status_more_opt,font='Verdana 10 bold')
more_opt_button.grid(row=1,pady=10)

frame_more_opt = tk_Frame(master=frame_differences)

#representation input

label_repres = tk_Label(frame_more_opt, text='Representation:',anchor='w',font='Verdana 10 bold')
label_repres.grid(row=0, column=0,padx=10,pady=5)

repres_rel = tk_IntVar()
repres_abs = tk_IntVar()

check_repres_rel = tk_Checkbutton(frame_more_opt,text='relative',variable=repres_rel, onvalue=1, offvalue=0)
check_repres_rel.grid(row=0, column=1,padx=10)

check_repres_abs = tk_Checkbutton(frame_more_opt,text='absolute values',variable=repres_abs,onvalue=1, offvalue=0)
check_repres_abs.grid(row=0, column=4,padx=10)

#colors to be completed

label_cmap = tk_Label(frame_more_opt, text='Color Map:',anchor='w',font='Verdana 10 bold')
label_cmap.grid(row=1, column=0,padx=10,pady=5)

default_cmap = tk_StringVar(root, value='viridis')

select_cmap = ttk.Combobox(frame_more_opt, values=['viridis'],textvariable=default_cmap)
select_cmap.grid(row=1,column=1,padx=10)

check_cmin_var = tk_IntVar()

check_cmin = tk_Checkbutton(frame_more_opt, text='min color:',variable=check_cmin_var, onvalue=1, offvalue=0,font='Verdana 10')
check_cmin.grid(row=1, column=3,padx=10)

entry_cmin = tk_Entry(frame_more_opt)
entry_cmin.grid(row=1, column=4,padx=10)

unit_cmin = tk_Label(frame_more_opt, text='cm/%',anchor='w',font='Verdana 10')
unit_cmin.grid(row=1, column=5)

check_cmax_var = tk_IntVar()

check_cmax = tk_Checkbutton(frame_more_opt, text='max color:',variable=check_cmax_var, onvalue=1, offvalue=0,font='Verdana 10')
check_cmax.grid(row=1, column=6,padx=10)

entry_cmax = tk_Entry(frame_more_opt)
entry_cmax.grid(row=1, column=7,padx=10)

unit_cmax = tk_Label(frame_more_opt, text='cm/%',anchor='w',font='Verdana 10')
unit_cmax.grid(row=1, column=8)

#upper limit

label_limit_plot = tk_Label(frame_more_opt, text='Limit Plot:',anchor='w',font='Verdana 10 bold')
label_limit_plot.grid(row=2, column=0,padx=10,pady=5)

check_limit_var = tk_IntVar()

check_limit = tk_Checkbutton(frame_more_opt, text='show',variable=check_limit_var, onvalue=1, offvalue=0,font='Verdana 10')
check_limit.grid(row=2, column=1,padx=10)

label_limit = tk_Label(frame_more_opt, text='Limit:',anchor='w',font='Verdana 10')
label_limit.grid(row=2, column=3,padx=10)

default_limit = tk_StringVar(root, value='0.5')

entry_limit = tk_Entry(frame_more_opt,textvariable=default_limit)
entry_limit.bind("<Return>",show_diff)
entry_limit.grid(row=2, column=4)

unit_limit = tk_Label(frame_more_opt, text='cm/%',anchor='w',font='Verdana 10')
unit_limit.grid(row=2, column=5)


##OUPUT##

frame_output = tk_Frame(master=frame_differences)
#frame_output.pack(side='top',pady=10)
frame_output.grid(row=3,pady=10)

frame_dif_ster = tk_Frame(master=frame_output)
frame_dif_ster.pack(side='left')
frame_dif_geo = tk_Frame(master=frame_output)
frame_dif_geo.pack(side='left',padx=10)


#add calculation button to start calculation of equivalent square

show_button = tk_Button(frame_differences, text="Show", fg="black",command=show_diff,font='Verdana 10 bold')
#show_button.pack(side='top',pady=10)
show_button.grid(row=4,pady=10)

s = ttk.Style()
s.theme_use('default')
#s.configure("red.Horizontal.TProgressbar", foreground='green', background='green',lightcolor='green',darkcolor='green')
s.configure("red.Horizontal.TProgressbar", background='#29931d')


progress = ttk.Progressbar(tab_differences, orient = 'horizontal',length = 100, mode = 'determinate',style="red.Horizontal.TProgressbar")
progress.pack(side='bottom',padx=10,fill='both',pady=10)

frame_note_d = tk_Frame(master=tab_differences)
frame_note_d.pack(side='bottom',padx=10,fill='both',pady=10)

note_1_d = tk_Label(frame_note_d,text='*Only relevant for modes including FFF-fields.',anchor='w',font='Verdana 9')
note_1_d.pack(side='top',fill='both')
note_2_d = tk_Label(frame_note_d,text='Selection of the implemented FFF-profiles, used as weighting function in the calculation. Obtained for ELEKTA Versa HD LINACS at a depth of 1cm.',anchor='w',font='Verdana 9')
note_2_d.pack(side='top',fill='both')



### TAB: CREATE TABLES ###

##INPUT##

#BASIC#

#add and pack differences frame
frame_tables = tk_Frame(master=tab_tables)
frame_tables.pack(side='top', padx='20', pady='5')

#add and pack frame_x_y to master frame: using for x, y, z and tpr2010 inputs
frame_input_t = tk_Frame(master=frame_tables)
frame_input_t.grid(row=0,pady=10)

#min input
label_min_t = tk_Label(frame_input_t, text='Smin:',anchor='w',font='Verdana 10 bold')
label_min_t.grid(row=0, column=0,padx=10)

default_min_t = tk_StringVar(root, value='0.1')

entry_min_t = tk_Entry(frame_input_t,textvariable=default_min_t)
entry_min_t.bind("<Return>",show_table)
entry_min_t.grid(row=0, column=1)

unit_min_t = tk_Label(frame_input_t, text='cm',anchor='w',font='Verdana 10')
unit_min_t.grid(row=0, column=2)

#max input
label_max_t = tk_Label(frame_input_t, text='Smax:',anchor='w',font='Verdana 10 bold')
label_max_t.grid(row=0, column=3,padx=10)

default_max_t = tk_StringVar(root, value='4')

entry_max_t = tk_Entry(frame_input_t,textvariable=default_max_t)
entry_max_t.bind("<Return>",show_table)
entry_max_t.grid(row=0, column=4)

unit_max_t = tk_Label(frame_input_t, text='cm',anchor='w',font='Verdana 10')
unit_max_t.grid(row=0, column=5)


#inc input
label_inc_t = tk_Label(frame_input_t, text='Increment:',anchor='w',font='Verdana 10 bold')
label_inc_t.grid(row=0, column=6,padx=10)

default_inc_t = tk_StringVar(root, value='0.1')

entry_inc_t = tk_Entry(frame_input_t,textvariable=default_inc_t)
entry_inc_t.bind("<Return>",show_table)
entry_inc_t.grid(row=0, column=7)

unit_inc_t = tk_Label(frame_input_t, text='cm',anchor='w',font='Verdana 10')
unit_inc_t.grid(row=0, column=8)


#z input
label_z_d_t = tk_Label(frame_input_t, text='depth:',anchor='w',font='Verdana 10 bold')
label_z_d_t.grid(row=1, column=0,pady=5,padx=10)

default_z_d_t = tk_StringVar(root, value='10')

entry_z_d_t = tk_Entry(frame_input_t,textvariable=default_z_d_t)
entry_z_d_t.bind("<Return>",show_table)
entry_z_d_t.grid(row=1, column=1)

unit_z_d_t = tk_Label(frame_input_t, text='cm',anchor='w',font='Verdana 10')
unit_z_d_t.grid(row=1, column=2)


#tpr2010 input
label_tpr_d_t = tk_Label(frame_input_t, text='TPR2010:',anchor='w',font='Verdana 10 bold')
label_tpr_d_t.grid(row=1, column=3,padx=10)

default_tpr_d_t = tk_StringVar(root, value='0.671')

entry_tpr_d_t = tk_Entry(frame_input_t,textvariable=default_tpr_d_t)
entry_tpr_d_t.bind("<Return>",show_table)
entry_tpr_d_t.grid(row=1, column=4)


#Mode input
label_modus_d_t = tk_Label(frame_input_t, text='mode:',anchor='w',font='Verdana 10 bold')
label_modus_d_t.grid(row=1, column=6,padx=10)

default_modus_d_t = tk_StringVar(root, value='WFF-WFF')

entry_modus_d_t = ttk.Combobox(frame_input_t,values=['WFF-WFF','FFF-FFF','FFF-WFF'],textvariable=default_modus_d_t)
entry_modus_d_t.bind("<Return>",show_table)
entry_modus_d_t.grid(row=1, column=7)

#Method input
label_method_d_t = tk_Label(frame_input_t, text='method:',anchor='w',font='Verdana 10 bold')
label_method_d_t.grid(row=2, column=0,padx=10)

default_method_d_t = tk_StringVar(root, value='Axis dose')

entry_method_d_t = ttk.Combobox(frame_input_t,values=['Axis dose','TPR2010'],textvariable=default_method_d_t)
entry_method_d_t.bind("<Return>",show_table)
entry_method_d_t.grid(row=2, column=1)

#FFF Profile input
label_profile_d_t = tk_Label(frame_input_t, text='FFF profile*:',anchor='w',font='Verdana 10 bold')
label_profile_d_t.grid(row=2, column=3,padx=10)

default_profile_d_t = tk_StringVar(root, value='6mv')

entry_profile_d_t = ttk.Combobox(frame_input_t,values=['6mv','10mv'],textvariable=default_profile_d_t)
entry_profile_d_t.bind("<Return>",show_table)
entry_profile_d_t.grid(row=2, column=4)



#Buttons
button_frame = tk_Frame(master=frame_tables)
button_frame.grid(row=2,pady=10)

show_button_t = tk_Button(button_frame, text="Show", fg="black",command=show_table,font='Verdana 10 bold')
show_button_t.grid(row=0,column=0,pady=10)

save_button_t = tk_Button(button_frame, text="Save", fg="black",command=save_excel,font='Verdana 10 bold')
save_button_t.grid(row=0,column=1,pady=10)

##OUPUT##

frame_table_output = tk_Frame(master=frame_tables)
frame_table_output.grid(row=3,pady=10,sticky='ew')


progress_tab = ttk.Progressbar(tab_tables, orient = 'horizontal',length = 100, mode = 'determinate',style="red.Horizontal.TProgressbar")
progress_tab.pack(side='bottom',padx=10,fill='both',pady=10)

frame_note_t = tk_Frame(master=tab_tables)
frame_note_t.pack(side='bottom',padx=10,fill='both',pady=10)

note_1_t = tk_Label(frame_note_t,text='*Only relevant for modes including FFF-fields.',anchor='w',font='Verdana 9')
note_1_t.pack(side='top',fill='both')
note_2_t = tk_Label(frame_note_t,text='Selection of the implemented FFF-profiles, used as weighting function in the calculation. Obtained for ELEKTA Versa HD LINACS at a depth of 1cm.',anchor='w',font='Verdana 9')
note_2_t.pack(side='top',fill='both')


#start mainloop
root.mainloop()
