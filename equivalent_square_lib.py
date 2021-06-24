from scipy import integrate
from scipy.optimize import curve_fit
import numpy as np
from scipy.misc import derivative


#Functions to Calculate the beam kernel parameters Az, az, Bz, bz using the beam-quality factor TPR20,10 (here called t) according to Nyholm et. all. 2006

def A1(t):
	return 0.0128018 - 0.0577391*t + 0.1790839*t**2 - 0.2467955*t**3 + 0.1328192*t**4 - 0.0194684*t**5

def A2(t):
	return 16.7815028 - 279.4672663*t + 839.0016549*t**2 - 978.4915013*t**3 + 470.5317337*t**4 - 69.2485573*t**5

def A3(t):
	return -0.0889669 - 0.2587584*t + 0.7069203*t**2 - 0.3654033*t**3 + 0.0029760*t**4 - 0.0003786*t**5

def A4(t):
	return 0.0017089 - 0.0169150*t + 0.0514650*t**2 - 0.0639530*t**3 + 0.0324490*t**4 - 0.0049121*t**5

def A5(t):
	return 0.1431447 - 0.2134626*t + 0.5825546*t**2 - 0.2969273*t**3 - 0.0011436*t**4 + 0.0002219*t**5

def B1(t):
	return -42.7607523 + 264.3424720*t - 633.4540368*t**2 + 731.5311577*t**3 - 402.5280374*t**4 + 82.4936551*t**5

def B2(t):
	return 0.2428359 - 2.5029336*t + 7.6128101*t**2 - 9.5273454*t**3 + 4.8249840*t**4 - 0.7097852*t**5

def B3(t):
	return -0.0910420 - 0.2621605*t + 0.7157244*t**2 - 0.3664126*t**3 + 0.0000930*t**4 - 0.0000232*t**5

def B4(t):
	return 0.0017284 - 0.0172146*t + 0.0522109*t**2 - 0.0643946*t**3 + 0.0322177*t**4 - 0.0047015*t**5

def B5(t):
	return -30.4609625 + 354.2866078*t - 1073.2952368*t**2 + 1315.2670101*t**3 - 656.3702845*t**4 + 96.5983711*t**5

def a1(t):
	return -0.0065985 + 0.0242136*t - 0.0647001*t**2 + 0.0265272*t**3 + 0.0072169*t**4 - 0.0020479*t**5

def a2(t):
	return -26.3337419 + 435.6865552*t - 1359.8342546*t**2 + 1724.6602381*t**3 - 972.7565415*t**4 + 200.3468023*t**5

def b1(t):
	return -80.7027159 + 668.1710175*t - 2173.2445309*t**2 + 3494.2393490*t**3 - 2784.4670834*t**4 + 881.2276510*t**5

def b2(t):
	return 3.4685991 - 41.2468479*t + 124.9729952*t**2 - 153.2610078*t**3 + 76.5242757*t**4 - 11.2624113*t**5

def b3(t):
	return -39.6550497 + 277.7202038*t - 777.0749505*t**2 + 1081.5724508*t**3 - 747.1056558*t**4 + 204.5432666*t**5

def b4(t):
	return 0.6514859 - 4.7179961*t + 13.6742202*t**2 - 19.7521659*t**3 + 14.1873606*t**4 - 4.0478845*t**5

def b5(t):
	return 0.4695047 - 3.6644336*t + 10.0039321*t**2 - 5.1195905*t**3 - 0.0007387*t**4 + 0.0002360*t**5




def az(z,t):
	return a2(t) + a1(t)*z 

def bz(z,t):
	return (b1(t)*(1 - np.exp(b2(t)*(z**2 + b5(t)**2)**0.5))*np.exp(b3(t)*z + b4(t)*z**2))

def Az(z,t):
	return (A1(t)*(1 - np.exp(A2(t)*(z**2 + A5(t)**2)**0.5))*np.exp(A3(t)*z + A4(t)*z**2))*az(z,t)

def Bz(z,t):
	return (B1(t)*(1 - np.exp(B2(t)*(z**2 + B5(t)**2)**0.5))*np.exp(B3(t)*z + B4(t)*z**2))*bz(z,t)

##FLAT FIELDS##
#Functions to integrate the beam kernals for square and rectangular fields using the parameters Az, az, Bz, bz according to Ahnesjoe et. all. 1992 / Nyholm et. all. 2006


def Pz_tri_part_b(x,A,a,ankat):
	return A/a*np.exp(-a*ankat/(np.cos(x)))

def Pz_tri_part_a_int(A,a,theta):
	return A*theta/a

def Pz_tri_part_b_int(A,a,ankat,theta):
	return integrate.quad(Pz_tri_part_b, 0, theta, args=(A,a,ankat))

def Pz_tri_part_ges_int(A,a,ankat,theta):
	part_b_ges = Pz_tri_part_b_int(A,a,ankat,theta)
	result = Pz_tri_part_a_int(A,a,theta) - part_b_ges[0]
	error = part_b_ges[1]
	return (result, error)

def Pz_tri_ges_int(A,a,B,b,ankat,theta):
	A_ges = Pz_tri_part_ges_int(A,a,ankat,theta)
	B_ges = Pz_tri_part_ges_int(B,b,ankat,theta)
	result = A_ges[0] + B_ges[0]
	error = ((A_ges[1])**2 + (B_ges[1])**2)**(1/2)
	return (result,error)

def Dz_square(A,a,B,b,d):
	tri_ges = Pz_tri_ges_int(A,a,B,b,d/2,np.pi/4)
	result = 8*tri_ges[0]
	error = 8*tri_ges[1]
	return (result, error)

def Dz_rect(A,a,B,b,dx,dy):
	theta_a = np.arctan(dy/dx)
	theta_b = np.arctan(dx/dy)

	tri_a = Pz_tri_ges_int(A,a,B,b,dx/2,theta_a)
	tri_b = Pz_tri_ges_int(A,a,B,b,dy/2,theta_b)

	result = 4*tri_a[0] + 4*tri_b[0]
	error = ((4*tri_a[1])**2 + (4*tri_b[1])**2)**(1/2)

	return (result,error)

def Dz_square_tpr(d,z,t):
	A = Az(z,t)
	a = az(z,t)
	B = Bz(z,t)
	b = bz(z,t)
	return Dz_square(A,a,B,b,d)

def Dz_rect_tpr(dx,dy,z,t):
	A = Az(z,t)
	a = az(z,t)
	B = Bz(z,t)
	b = bz(z,t)
	return Dz_rect(A,a,B,b,dx,dy)

#newton optimation

def newton_f(d,z,t,k_rect):
	return Dz_square_tpr(d,z,t)[0] - k_rect

def newton_df(d,f_d,z,t,k_rect):
	h = 1e-6
	return (newton_f(d + h,z,t,k_rect) - f_d)/h


def newton_df_scipy(d,z,t,k_rect):
	return derivative(newton_f,x0=d,dx=1e-6,args=(z,t,k_rect))

#Sterling equation

def sterling(dx,dy):
	return 4*dx*dy/(2*(dx+dy))



def newton(f,Df,x0,epsilon,max_iter,args):
	xn = x0
	for n in range(0,max_iter):
		fxn = f(xn,*args)
		if abs(fxn) < epsilon:
			return xn
		Dfxn = Df(xn,fxn,*args)
		if Dfxn == 0:
			print('Zero derivative. No solution found.')
			return None
		xn = xn - fxn/Dfxn
	print('Exceeded maximum iterations. No solution found.')
	return None


#Object to calculate the equivalent square by comparing rectangular field kernels with square field kernels

#Equivalent Square WFF/WFF
class EquivalentSquare:
	def __init__(self,dx,dy,z=10,tpr2010=0.671,epsilon=0.000001,max_iter=100,no_kernel=True):
		self.dx = dx
		self.dy = dy
		self.z = z
		self.tpr2010 = tpr2010
		self.epsilon = epsilon
		self.max_iter = max_iter
		self.no_kernel = no_kernel
		self.update()


	#SET FUNCTIONS

	def update(self):
		dx = self.dx
		dy = self.dy
		z = self.z
		t = self.tpr2010
		epsilon = self.epsilon
		max_iter = self.max_iter
		no_kernel = self.no_kernel
		

		re_res, re_err = Dz_rect_tpr(dx,dy,z,t)

		round_digits = 2

		geo_mean = round((dx*dy)**0.5,round_digits)
		sterling =round(4*dx*dy/(2*(dx+dy)),round_digits)

		self.geometric_mean = geo_mean
		self.sterling = sterling
		self.kernel_re = re_res
		self.kernal_re_er = re_err

		if dx==dy:
			self.equi_sq_raw = dx
			self.equi_sq = dx
			self.kernel_sq = re_res
			self.kernel_sq_er = re_err
			self.kernel_dif = 0
		else:

			equi_newton = newton(f=newton_f,Df=newton_df,x0=sterling,epsilon=epsilon,max_iter=max_iter,args=(z,t,re_res))
			self.equi_sq_raw = equi_newton
			self.equi_sq = round(equi_newton,round_digits)

			if no_kernel:
				self.kernel_sq = 0
				self.kernel_sq_er = 0
				self.kernel_dif = 0
			else:
				sq_res, sq_err = Dz_square_tpr(equi_newton,z,t)

				self.kernel_sq = sq_res
				self.kernel_sq_er = sq_err
				self.kernel_dif = re_res - sq_res

		#self.geo_dif = round(abs(self.geometric_mean-self.equi_sq),round_digits)
		self.geo_dif = round(self.geometric_mean-self.equi_sq,round_digits)
		self.geo_dif_rel = round(self.geo_dif/self.equi_sq,4)
		#self.sterling_dif = round(abs(self.sterling-self.equi_sq),round_digits)
		self.sterling_dif = round(self.sterling-self.equi_sq,round_digits)
		self.sterling_dif_rel = round(self.sterling_dif/self.equi_sq,4)


	#GET FUNCTIONS:

	#Properties

	@property
	def geo_all(self):
		return (self.geometric_mean,self.geo_dif,self.geo_dif_rel)

	@property
	def sterling_all(self):
		return (self.sterling,self.sterling_dif,self.sterling_dif_rel)

	@property
	def input(self):
		return ([self.dx,self.dy],self.z,self.tpr2010)

	@property   
	def overview(self):
		return (self.equi_sq,self.geo_all,self.sterling_all,self.input)




##FFF FIELDS

#FFF Kernel
def gaussian(x,a,b):
	return a*np.exp(-x**2/(2*b**2))

def sigmoid(x):
	return 1/(1+np.exp(-np.array(x)))

#Extended sigmoid function: with hight 'a', slope 'b' and x-shift 'c'
def sig_ext(x,a,b,c):
	return a*sigmoid(b*(x-c))

def sig_ext_sym(x,a,b,c,d,e):
	return sig_ext(x,a,b,-c+d) + sig_ext(-x,a,b,-c-d) -a + e

def sig_gaussian2(x,a,b,c,e,f,g,h,i):
	return sig_ext_sym(x,a,b,c,0,e) + gaussian(x,f,g) + gaussian(x,h,i)



def Pz_pre_int(r,A,a):
	return A*np.exp(-a*r)


def W_pre_int(x,energy='6mv'):
	if energy=='6mv':
		return sig_gaussian2(x,5.15229707e-01,4.70403758e+00,2.00309769e+01,1.56120538e-02,5.73646811e-01,9.94490425e+00,1.26526257e-01,-3.44333164e+00)
	else:
		return sig_gaussian2(x,3.71770516e-01,4.22423908e+00,2.00093005e+01,1.25045275e-02,5.27197010e-01,9.24571498e+00,1.76659943e-01,-3.12352475e+00)

def Pz_w_pre_int(r,A,a,energy='6mv'):
	n_fac = 1/W_pre_int(0,energy)
	return Pz_pre_int(r,A,a)*W_pre_int(r,energy)*n_fac

def Pz_w_int_r(phi,A,a,L,energy='6mv'):
	return  integrate.quad(Pz_w_pre_int,0,L/(np.cos(phi)),args=(A,a,energy))[0]

def Pz_w_int_phi(theta,A,a,L,energy='6mv'):
	return integrate.quad(Pz_w_int_r,0,theta,args=(A,a,L,energy))[0]

def Pz_w_int_ges_tri(A,a,B,b,L,theta,energy='6mv'):
	A_ges = Pz_w_int_phi(theta,A,a,L,energy)
	B_ges = Pz_w_int_phi(theta,B,b,L,energy)

	return A_ges + B_ges

def Dz_w_square(A,a,B,b,d,energy='6mv'):
	tri_ges = Pz_w_int_ges_tri(A,a,B,b,d/2,np.pi/4,energy)
	return 8*tri_ges

def Dz_w_rect(A,a,B,b,dx,dy,energy='6mv'):
	theta_a = np.arctan(dy/dx)
	theta_b = np.arctan(dx/dy)

	tri_a = Pz_w_int_ges_tri(A,a,B,b,dx/2,theta_a,energy)
	tri_b = Pz_w_int_ges_tri(A,a,B,b,dy/2,theta_b,energy)

	return 4*tri_a + 4*tri_b

def Dz_w_square_tpr(d,z,t,energy='6mv'):
	A = Az(z,t)
	a = az(z,t)
	B = Bz(z,t)
	b = bz(z,t)
	return (Dz_w_square(A,a,B,b,d,energy),0)

def Dz_w_rect_tpr(dx,dy,z,t,energy='6mv'):
	A = Az(z,t)
	a = az(z,t)
	B = Bz(z,t)
	b = bz(z,t)
	return (Dz_w_rect(A,a,B,b,dx,dy,energy),0)

def newton_f_fff(d,z,t,k_rect,energy='6mv'):
	return Dz_w_square_tpr(d,z,t,energy)[0] - k_rect

def newton_df_fff(d,f_d,z,t,k_rect,energy='6mv'):
	h = 1e-6
	return (newton_f_fff(d + h,z,t,k_rect,energy) - f_d)/h

#Equivalent Square FFF/FFF
class EquivalentSquareFFF:
	def __init__(self,dx,dy,z=10,tpr2010=0.671,epsilon=0.000001,max_iter=100,no_kernel=True,energy='6mv'):
		self.dx = dx
		self.dy = dy
		self.z = z
		self.tpr2010 = tpr2010
		self.epsilon = epsilon
		self.max_iter = max_iter
		self.no_kernel = no_kernel
		self.energy = energy
		self.update()


	#SET FUNCTIONS

	def update(self):
		dx = self.dx
		dy = self.dy
		z = self.z
		t = self.tpr2010
		epsilon = self.epsilon
		max_iter = self.max_iter
		no_kernel = self.no_kernel
		energy = self.energy
		

		re_res, re_err = Dz_w_rect_tpr(dx,dy,z,t,energy)

		round_digits = 2

		geo_mean = round((dx*dy)**0.5,round_digits)
		sterling =round(4*dx*dy/(2*(dx+dy)),round_digits)

		self.geometric_mean = geo_mean
		self.sterling = sterling
		self.kernel_re = re_res
		self.kernal_re_er = re_err

		if dx==dy:
			self.equi_sq_raw = dx
			self.equi_sq = dx
			self.kernel_sq = re_res
			self.kernel_sq_er = re_err
			self.kernel_dif = 0
		else:

			equi_newton = newton(f=newton_f_fff,Df=newton_df_fff,x0=sterling,epsilon=epsilon,max_iter=max_iter,args=(z,t,re_res))
			self.equi_sq_raw = equi_newton
			self.equi_sq = round(equi_newton,round_digits)

			if no_kernel:
				self.kernel_sq = 0
				self.kernel_sq_er = 0
				self.kernel_dif = 0
			else:
				sq_res, sq_err = Dz_w_square_tpr(equi_newton,z,t,energy)

				self.kernel_sq = sq_res
				self.kernel_sq_er = sq_err
				self.kernel_dif = re_res - sq_res

		self.geo_dif = round(self.geometric_mean-self.equi_sq,round_digits)
		self.geo_dif_rel = round(self.geo_dif/self.equi_sq,4)
		self.sterling_dif = round(self.sterling-self.equi_sq,round_digits)
		self.sterling_dif_rel = round(self.sterling_dif/self.equi_sq,4)


	#GET FUNCTIONS:

	#Properties

	@property
	def geo_all(self):
		return (self.geometric_mean,self.geo_dif,self.geo_dif_rel)

	@property
	def sterling_all(self):
		return (self.sterling,self.sterling_dif,self.sterling_dif_rel)

	@property
	def input(self):
		return ([self.dx,self.dy],self.z,self.tpr2010)

	@property	
	def overview(self):
		return (self.equi_sq,self.geo_all,self.sterling_all,self.input)


#Equivalent Square FFF/WFF
class EquivalentSquareFFFWFF:
	def __init__(self,dx,dy,z=10,tpr2010=0.671,epsilon=0.000001,max_iter=1000,no_kernel=True,energy='6mv'):
		self.dx = dx
		self.dy = dy
		self.z = z
		self.tpr2010 = tpr2010
		self.epsilon = epsilon
		self.max_iter = max_iter
		self.no_kernel = no_kernel
		self.energy = energy
		self.update()


	#SET FUNCTIONS

	def update(self):
		dx = self.dx
		dy = self.dy
		z = self.z
		t = self.tpr2010
		epsilon = self.epsilon
		max_iter = self.max_iter
		no_kernel = self.no_kernel
		energy = self.energy
		
		

		re_res, re_err = Dz_w_rect_tpr(dx,dy,z,t,energy)

		round_digits = 2

		geo_mean = round((dx*dy)**0.5,round_digits)
		sterling =round(4*dx*dy/(2*(dx+dy)),round_digits)

		self.geometric_mean = geo_mean
		self.sterling = sterling
		self.kernel_re = re_res
		self.kernal_re_er = re_err

		#if dx==dy:
		if False:
			self.equi_sq_raw = dx
			self.equi_sq = dx
			self.kernel_sq = re_res
			self.kernel_sq_er = re_err
			self.kernel_dif = 0
		else:
			equi_newton = newton(f=newton_f,Df=newton_df,x0=sterling,epsilon=epsilon,max_iter=max_iter,args=(z,t,re_res))
			self.equi_sq_raw = equi_newton
			self.equi_sq = round(equi_newton,round_digits)

			if no_kernel:
				self.kernel_sq = 0
				self.kernel_sq_er = 0
				self.kernel_dif = 0
			else:
				sq_res, sq_err = Dz_square_tpr(equi_newton,z,t)

				self.kernel_sq = sq_res
				self.kernel_sq_er = sq_err
				self.kernel_dif = re_res - sq_res

		self.geo_dif = round(self.geometric_mean-self.equi_sq,round_digits)
		self.geo_dif_rel = round(self.geo_dif/self.equi_sq,4)
		self.sterling_dif = round(self.sterling-self.equi_sq,round_digits)
		self.sterling_dif_rel = round(self.sterling_dif/self.equi_sq,4)


	#GET FUNCTIONS:

	#Properties

	@property
	def geo_all(self):
		return (self.geometric_mean,self.geo_dif,self.geo_dif_rel)

	@property
	def sterling_all(self):
		return (self.sterling,self.sterling_dif,self.sterling_dif_rel)

	@property
	def input(self):
		return ([self.dx,self.dy],self.z,self.tpr2010)

	@property	
	def overview(self):
		return (self.equi_sq,self.geo_all,self.sterling_all,self.input)


##Calculating Equivalent Squares Using same TPR2010

def calc_tpr(f,dx,dy=0,tpr2010=0.775):
	if dy==0:
		d_20 = f(dx,20,tpr2010)[0]
		d_10 = f(dx,10,tpr2010)[0]
	else:
		d_20 = f(dx,dy,20,tpr2010)[0]
		d_10 = f(dx,dy,10,tpr2010)[0]

	return d_20/d_10

def newton_tpr_f(d,f_sq,res_re,tpr2010=0.775):
	return calc_tpr(f_sq,d,tpr2010=tpr2010) - res_re

def newton_tpr_df(d,res_f,f_sq,res_re,tpr2010=0.775):
	h = 1e-6
	return (newton_tpr_f(d + h,f_sq,res_re,tpr2010) - res_f)/h


def calc_tpr_e(f,dx,dy=0,tpr2010=0.775,energy='6mv'):
	if dy==0:
		d_20 = f(dx,20,tpr2010,energy)[0]
		d_10 = f(dx,10,tpr2010,energy)[0]
	else:
		d_20 = f(dx,dy,20,tpr2010,energy)[0]
		d_10 = f(dx,dy,10,tpr2010,energy)[0]

	return d_20/d_10












#Equivalent Square WFF/WFF: Using same TPR2010 as Equivalent Square definition
class EquivalentSquareTPR:
	def __init__(self,dx,dy,z=10,tpr2010=0.671,epsilon=0.000001,max_iter=100,no_kernel=False):
		self.dx = dx
		self.dy = dy
		self.tpr2010 = tpr2010
		self.epsilon = epsilon
		self.max_iter = max_iter
		self.no_kernel = no_kernel
		self.update()


	#SET FUNCTIONS

	def update(self):
		dx = self.dx
		dy = self.dy
		t = self.tpr2010
		epsilon = self.epsilon
		max_iter = self.max_iter
		no_kernel = self.no_kernel
		

		re_res = calc_tpr(Dz_rect_tpr,dx,dy,tpr2010=t)

		round_digits = 2

		geo_mean = round((dx*dy)**0.5,round_digits)
		sterling =round(4*dx*dy/(2*(dx+dy)),round_digits)

		self.geometric_mean = geo_mean
		self.sterling = sterling
		self.tpr = re_res

		if dx==dy:
			self.equi_sq_raw = dx
			self.equi_sq = dx
			self.tpr_sq = re_res
			self.tpr_dif = 0
		else:

			equi_newton = newton(f=newton_tpr_f,Df=newton_tpr_df,x0=sterling,epsilon=epsilon,max_iter=max_iter,args=(Dz_square_tpr,re_res,t))
			self.equi_sq_raw = equi_newton
			self.equi_sq = round(equi_newton,round_digits)

			if no_kernel:
				self.tpr_sq = 0
				self.tpr_dif = 0
			else:
				sq_res = calc_tpr(Dz_square_tpr,equi_newton,tpr2010=t)

				self.tpr_sq = sq_res
				self.tpr_dif = re_res - sq_res

		self.geo_dif = round(self.geometric_mean-self.equi_sq,round_digits)
		self.geo_dif_rel = round(self.geo_dif/self.equi_sq,4)
		self.sterling_dif = round(self.sterling-self.equi_sq,round_digits)
		self.sterling_dif_rel = round(self.sterling_dif/self.equi_sq,4)


	#GET FUNCTIONS:

	#Properties

	@property
	def geo_all(self):
		return (self.geometric_mean,self.geo_dif,self.geo_dif_rel)

	@property
	def sterling_all(self):
		return (self.sterling,self.sterling_dif,self.sterling_dif_rel)

	@property
	def input(self):
		return ([self.dx,self.dy],self.z,self.tpr2010)

	@property   
	def overview(self):
		return (self.equi_sq,self.geo_all,self.sterling_all,self.input)



#Equivalent Square FFF/FFF: Using same TPR2010 as Equivalent Square definition
class EquivalentSquareFFFTPR:
	def __init__(self,dx,dy,z=10,tpr2010=0.671,epsilon=0.000001,max_iter=100,no_kernel=False,energy='6mv'):
		self.dx = dx
		self.dy = dy
		self.tpr2010 = tpr2010
		self.epsilon = epsilon
		self.max_iter = max_iter
		self.no_kernel = no_kernel
		self.energy = energy
		self.update()


	#SET FUNCTIONS

	def update(self):
		dx = self.dx
		dy = self.dy
		t = self.tpr2010
		epsilon = self.epsilon
		max_iter = self.max_iter
		no_kernel = self.no_kernel
		energy = self.energy
		

		re_res = calc_tpr_e(Dz_w_rect_tpr,dx,dy,tpr2010=t,energy=energy)

		round_digits = 2

		geo_mean = round((dx*dy)**0.5,round_digits)
		sterling =round(4*dx*dy/(2*(dx+dy)),round_digits)

		self.geometric_mean = geo_mean
		self.sterling = sterling
		self.tpr = re_res

		if dx==dy:
			self.equi_sq_raw = dx
			self.equi_sq = dx
			self.tpr_sq = re_res
			self.tpr_dif = 0
		else:

			equi_newton = newton(f=newton_tpr_f,Df=newton_tpr_df,x0=sterling,epsilon=epsilon,max_iter=max_iter,args=(Dz_w_square_tpr,re_res,t,energy))
			self.equi_sq_raw = equi_newton
			self.equi_sq = round(equi_newton,round_digits)

			if no_kernel:
				self.tpr_sq = 0
				self.tpr_dif = 0
			else:
				sq_res = calc_tpr(Dz_square_tpr,equi_newton,tpr2010=t)

				self.tpr_sq = sq_res
				self.tpr_dif = re_res - sq_res

		self.geo_dif = round(self.geometric_mean-self.equi_sq,round_digits)
		self.geo_dif_rel = round(self.geo_dif/self.equi_sq,4)
		self.sterling_dif = round(self.sterling-self.equi_sq,round_digits)
		self.sterling_dif_rel = round(self.sterling_dif/self.equi_sq,4)


	#GET FUNCTIONS:

	#Properties

	@property
	def geo_all(self):
		return (self.geometric_mean,self.geo_dif,self.geo_dif_rel)

	@property
	def sterling_all(self):
		return (self.sterling,self.sterling_dif,self.sterling_dif_rel)

	@property
	def input(self):
		return ([self.dx,self.dy],self.z,self.tpr2010)

	@property   
	def overview(self):
		return (self.equi_sq,self.geo_all,self.sterling_all,self.input)


#Equivalent Square FFF/FFF: Using same TPR2010 as Equivalent Square definition
class EquivalentSquareFFFWFFTPR:
	def __init__(self,dx,dy,z=10,tpr2010=0.671,epsilon=0.000001,max_iter=100,no_kernel=False,energy='6mv'):
		self.dx = dx
		self.dy = dy
		self.tpr2010 = tpr2010
		self.epsilon = epsilon
		self.max_iter = max_iter
		self.no_kernel = no_kernel
		self.energy = energy
		self.update()


	#SET FUNCTIONS

	def update(self):
		dx = self.dx
		dy = self.dy
		t = self.tpr2010
		epsilon = self.epsilon
		max_iter = self.max_iter
		no_kernel = self.no_kernel
		energy = self.energy
		

		re_res = calc_tpr_e(Dz_w_rect_tpr,dx,dy,tpr2010=t,energy=energy)

		round_digits = 2

		geo_mean = round((dx*dy)**0.5,round_digits)
		sterling =round(4*dx*dy/(2*(dx+dy)),round_digits)

		self.geometric_mean = geo_mean
		self.sterling = sterling
		self.tpr = re_res

		if False:
			self.equi_sq_raw = dx
			self.equi_sq = dx
			self.tpr_sq = re_res
			self.tpr_dif = 0
		else:

			equi_newton = newton(f=newton_tpr_f,Df=newton_tpr_df,x0=sterling,epsilon=epsilon,max_iter=max_iter,args=(Dz_square_tpr,re_res,t))
			self.equi_sq_raw = equi_newton
			self.equi_sq = round(equi_newton,round_digits)

			if no_kernel:
				self.tpr_sq = 0
				self.tpr_dif = 0
			else:
				sq_res = calc_tpr(Dz_square_tpr,equi_newton,tpr2010=t)

				self.tpr_sq = sq_res
				self.tpr_dif = re_res - sq_res

		self.geo_dif = round(self.geometric_mean-self.equi_sq,round_digits)
		self.geo_dif_rel = round(self.geo_dif/self.equi_sq,4)
		self.sterling_dif = round(self.sterling-self.equi_sq,round_digits)
		self.sterling_dif_rel = round(self.sterling_dif/self.equi_sq,4)


	#GET FUNCTIONS:

	#Properties

	@property
	def geo_all(self):
		return (self.geometric_mean,self.geo_dif,self.geo_dif_rel)

	@property
	def sterling_all(self):
		return (self.sterling,self.sterling_dif,self.sterling_dif_rel)

	@property
	def input(self):
		return ([self.dx,self.dy],self.z,self.tpr2010)

	@property   
	def overview(self):
		return (self.equi_sq,self.geo_all,self.sterling_all,self.input)