"""Atom class"""
import math


class Atom(object):
	
	#Get the symbol corresponding to Z number Znb
	def get_symbol(self):
		return {
		0 : 'Unknown',
		1 : 'H',
		2 : 'He',
		3 : 'Li',
		4 : 'Be',
		6 : 'C',
		7 : 'N',
		8 : 'O',
		9 : 'F',
		10: 'Ne',
		13: 'Al',
		15: 'P',
		16: 'S',
		17: 'Cl',
		18: 'Ar',
		23: 'V',
		24: 'Cr',
		25: 'Mn',
		26: 'Fe',
		27: 'Co',
		28: 'Ni',
		29: 'Cu',
		30: 'Zn',
		36: 'Kr',
		43: 'Tc',
		46: 'Pd',
		64: 'Gd',
		76: 'Os',
		78: 'Pt',
		79: 'Au',
		-1: 'X'
		}[self.Znb]
		
	def get_Znb(self):
		return {
		'Unknown':0,
		'H':1,
		'He':2,
		'Li':3,
		'Be':4,
		'C':6,
		'N':7,
		'O':8,
		'F':9,
		'Ne':10,
		'Al':13,
		'P':15,
		'S':16,
		'Cl':17,
		'Ar':18,
		'V':23,
		'Cr':24,
		'Mn':25,
		'Fe':26,
		'Co':27,
		'Ni':28,
		'Cu':29,
		'Zn':30,
		'Kr':36,
		'Tc':43,
		'Pd':46,
		'Gd':64,
		'Os':76,
		'Pt':78,
		'Au':79,
		'X': -1
		}[self.symbol]

	def get_radius(self):
		return {
		'Unknown':0,
		'H':5,
		'He':5,
		'Be':5,
		'N':10,
		'C':10,
		'O':10,
		'F':10,
		'Ne':10,
		'Al':12,
		'P':12,
		'S':12,
		'Cl':12,
		'V':15,
		'Mn':20,
		'Fe':20,
		'Co':22,
		'Ni':20,
		'Cu':20,
		'Zn':20,
		'Tc':10,
		'Pd':10,
		'Gd':10,
		'Os':10,
		'Pt':18,
		'Au':18,
		'X':0
		}[self.symbol]
	
	def get_color(self):
		return {
		'Unknown':[0,0,0],
		'H':[0.937,0.937,0.937],
		'He':[0.90,0.90,0.90],
		'Be':[0.50,0.50,0.00],
		'N':[0.373, 0.553, 0.827],
		'O':[0.784, 0.216, 0.216],
		'C':[0.302, 0.302, 0.302],
		'F':[0.0, 0.302, 0.784],
		'Ne':[0.0,0.25,0.75],
		'Al':[0.5,0.5,0.5],
		'P':[0.7,0.3,0.0],
		'S':[1.0, 0.8, 0.0],
		'Cl':[0.0,0.784,0.302],
		'V':[0.5,0.0,0.5],
		'Mn':[0.831, 0.667, 0.0],
		'Fe':[0.831, 0.667, 0.0],
		'Co':[0.931, 0.767, 0.0],
		'Ni':[0.831, 0.667, 0.0],
		'Cu':[0.831, 0.667, 0.0],
		'Zn':[0.831, 0.667, 0.0],
		'Au':[0.831, 0.667, 0.0],
		'Pd':[0.8, 0.8, 1.0],
		'Pt':[0.8, 0.8, 1.0],
		'X':[0.0, 0.0, 0.0]
		}[self.symbol]
	
	def get_color_blender(self):
		return {
		'Unknown':[0.0,0.0,0.0],
		'H':[0.925,0.925,0.925],
		'He':[0.90,0.90,0.90],
		'Be':[0.50,0.50,0.00],
		'N':[0.373,0.553,0.827],
		'O':[0.498,0.0,0.0],
		'C':[0.1,0.1,0.1],
		'F':[0.0,0.3,0.7],
		'Ne':[0.0,0.25,0.75],
		'P':[0.7,0.3,0.0],
		'S':[0.0,0.3,0.7],
		'Al':[0.3,0.3,0.3],
		'S':[1.0, 0.94, 0.0],
		'Cl':[0.0,0.7,0.3],
		'V':[0.5, 0.0, 0.5],
		'Mn':[0.831,0.667,0.0],
		'Fe':[0.831,0.667,0.0],
		'Co':[0.2,0.2,0.2],
		'Ni':[0.831,0.667,0.0],
		'Cu':[0.831,0.667,0.0],
		'Zn':[0.831,0.667,0.0],
		'Tc':[0.0,0.4,1.0],
		'Gd':[0.0,0.4,1.0],
		'Os':[0.0,0.4,1.0],
		'Au':[1.0,0.8,0.0],
		'Pd':[0.8,0.8,1.0],
		'Pt':[0.8,0.8,1.0],
		'X':[0.0,0.0,0.0]
		}[self.symbol]

	def min_coord(self, a2):
		xm=0.0
		ym=0.0
		#zm=0.0
		if self.x < a2.x:
			xm=self.x
		else:
			xm=a2.x
			
		if self.y < a2.y:
			ym=self.y
		else:
			ym=a2.y
			
		#if self.x < a2.z:
		#	zm=self.z
		#else:
		#	zm=a2.z
		
		m = [xm, ym]
		return m
	
	def max_coord(self, a2):
		xM=0.0
		yM=0.0
		#zm=0.0
		if self.x < a2.x:
			xM=a2.x
		else:
			xM=self.x
			
		if self.y < a2.y:
			yM=a2.y
		else:
			yM=self.y
			
		#if self.x < a2.z:
		#	zm=self.z
		#else:
		#	zm=a2.z
		
		M = [xM, yM]
		return M

		
			
	#Constructor method
	def __init__(self, index=0, Znb=0, x=0, y=0, z=0):
		self.index = index
		self.Znb=Znb
		self.symbol=self.get_symbol()
		self.x=x
		self.y=y
		self.z=z
	
	def __str__(self):
		return '%s%d %10.5f %10.5f %10.5f\n' % (self.symbol, self.index, self.x, self.y, self.z)
		
	def str_no_index(self):
		return '%s %10.5f %10.5f %10.5f\n' % (self.symbol, self.x, self.y, self.z)

	def str_for_cube(self):
		return '%d %10.5f %10.5f %10.5f %10.5f\n' % (self.Znb, self.Znb, self.x, self.y, self.z)
	
	#Curently only works for 2D molecules
	def str_for_ps(self, step=25, scale=1.0):
		radius = scale*self.get_radius()
		#if self.z<12.5:
		#	radius=radius-5*((12.5-abs(self.z))/12.5)
		#if self.z>12.5:
		#	radius=radius+5*((12.5-abs(self.z))/12.5)
		color = self.get_color()
		com = '%'
		symbNb='%s%d' % (self.symbol, self.index)
		col ='%.3f %.3f %.3f setrgbcolor\n' % (color[0], color[1], color[2])
		posRad = '%10.5f %10.5f\n%d\n0 360 arc fill\n\n' % (step*self.x, step*self.y, radius)
		#colText = ''
		##txt='0 setgray\n%10.5f %10.5f moveto\n/Helvetica findfont\n10 scalefont setfont\n(%s) show\n\n' % (step*self.x, step*self.y, symbNb)
		toReturn = com+symbNb+'\n'+col+posRad#+txt
		return toReturn
		
	def distance(self, a2):
		return math.sqrt((self.x-a2.x)**2 + (self.y-a2.y)**2 + (self.z-a2.z)**2)

	def distToPoint(self,x,y,z):
		dist = math.sqrt((self.x-x)**2 + (self.y-y)**2 + (self.z-z)**2)
		#txt='(%.2f %.2f) to (%.2f %.2f): %.2f A' % (self.x, self.y, x, y, dist)
		#print txt
		return dist

	#Curently only works for 2D
	def bond_for_ps(self, a2, step=25):
		dist = math.sqrt((self.x-a2.x)**2 + (self.y-a2.y)**2 + (self.z-a2.z)**2)
		tx = self.min_coord(a2)[0]+0.5*abs(self.x-a2.x)
		ty = self.min_coord(a2)[1]+0.5*abs(self.y-a2.y)
		com='%Bond between:'
		symbNb='%s%d - %s%d\nnewpath\n10.0 setlinewidth\n0.502 setgray\n' % (self.symbol, self.index, a2.symbol, a2.index)
		line = '%10.5f %10.5f moveto\n%10.5f %10.5f lineto\nstroke\n\n' % (step*self.x, step*self.y, step*a2.x, step*a2.y)
		#text=''#'%10.5f %10.5f moveto\n/Helvetica findfont\n10 scalefont setfont\n(%10.3f) show\n\n' % (25*a2.x, 25*a2.y, dist)
		##text2='0 setgray\n%10.5f %10.5f moveto\n/Helvetica findfont\n10 scalefont setfont\n(%10.3f) show\n\n' % (step*tx-20, step*ty, dist)
		toReturn = com+symbNb+line #+text2
		return toReturn

	def bond(self, a2):
		dist = math.sqrt((self.x-a2.x)**2 + (self.y-a2.y)**2 + (self.z-a2.z)**2)
		tx = self.min_coord(a2)[0]+0.5*abs(self.x-a2.x)
		ty = self.min_coord(a2)[1]+0.5*abs(self.y-a2.y)
		com='%Bond between:'
		symbNb='%s%d - %s%d\nnewpath\n10.0 setlinewidth\n0.5 setgray\n' % (self.symbol, self.index, a2.symbol, a2.index)
		line = '%10.5f %10.5f moveto\n%10.5f %10.5f lineto\nstroke\n\n' % (25*self.x, 25*self.y, 25*a2.x, 25*a2.y)
		#text=''#'%10.5f %10.5f moveto\n/Helvetica findfont\n10 scalefont setfont\n(%10.3f) show\n\n' % (25*a2.x, 25*a2.y, dist)
		#text2='0 setgray\n%10.5f %10.5f moveto\n/Helvetica findfont\n10 scalefont setfont\n(%10.3f) show\n\n' % (25*tx-20, 25*ty, dist)
		toReturn = com+symbNb+line
		return toReturn

	
	def rectangle_bond_for_ps(self, a2, step=5):
		dist = math.sqrt((self.x-a2.x)**2 + (self.y-a2.y)**2 + (self.z-a2.z)**2)
		com='%Bond between:'
		symbNb='%s%d - %s%d\nnewpath\n1.0 setlinewidth\n' % (self.symbol, self.index, a2.symbol, a2.index)
		line1 = '%10.5f %10.5f moveto\n%10.5f %10.5f lineto\n' % (25*self.x+step, 25*self.y+step, 25*self.x-step, 25*self.y-step)
		line2 = '%10.5f %10.5f lineto\n%10.5f %10.5f lineto\nclosepath\nstroke\n\n' % (25*a2.x-step, 25*a2.y-step, 25*a2.x+step, 25*a2.y+step)
		text='%10.5f %10.5f moveto\n/Helvetica findfont\n10 scalefont setfont\n(%10.3f) show\n\n' % (25*a2.x, 25*a2.y, dist)
		toReturn = com+symbNb+line1+line2+text
		return toReturn
	#def angle(self, a2, a3)
	
	def move(self, x, y, z):
		self.x=x
		self.y=y
		self.z=z
		
	def move_along(self, vx, vy, vz):
		newx = self.x+vx
		newy = self.y+vy
		newz = self.z+vz
		self.x=newx
		self.y=newy
		self.z=newz
	
	def reset_all(self,index,Znb,symbol,x,y,z):
		self.index=index
		self.Znb=Znb
		self.symbol=symbol
		self.x=x
		self.y=y
		self.z=z
	
	def reset(self, index, Znb, x, y, z):
		self.index = index
		self.Znb = Znb
		self.symbol = self.get_symbol()
		self.x = x
		self.y = y
		self.z = z
	
	#reset only the coordinates
	def reset_xyz(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z
	
	def reset_symb(self, index, symbol, x, y, z):
		self.index = index
		self.symbol = symbol
		self.Znb = self.get_Znb()
		self.x = x
		self.y = y
		self.z = z
	
	def copy(self, c):
		self.index=c.index
		self.Znb=c.Znb
		self.x=c.x
		self.y=c.y
		self.z=c.z
		self.symbol=c.symbol
		
	def mul(self, const):
		self.x = self.x * const
		self.y = self.y * const
		self.z = self.z * const
	
	def get_x(self):
		return self.x
	
	def get_y(self):
		return self.y
	
	def get_z(self):
		return self.z
