"""Molecular geometry class. It uses Atom objects. It is essentially an array of atoms
characterized by the index, symbol, Z number, and spatial coordinates.
The geometry can be read form a gaussian output file, written to a xyz file;
It is possible to add atoms or delete them"""

import Atoms
import random
import math
import bpy

def compare(numbers, size, num):
	ok = 1
	for i in range(size):
		if numbers[i] == num:
			ok = 0
	return ok

def gen_urand_num(start, stop, howMany):
	random.seed()
	numbers = []
	for i in range(howMany):
		numbers.append(random.randrange(start, stop+1))
		while(compare(numbers, i-1, numbers[i]) == 0):
			numbers[i] = random.randrange(start, stop+1)
	return numbers

def permute(a=[]):
	b=[]
	for i in range(0,len(a)+1):
        	b.append(0)
	for i in range(1, len(b)):
        	b[i]=a[i-1]
	b[0]=b[-1]
	del b[-1]
	return b

class Geometry(object):
	
	#CONSTRUCTOR method
	def __init__(self, NAtoms=0):
		self.list = []
		self.NAtoms = NAtoms
		for i in range(0, NAtoms):
			self.list.append(Atoms.Atom())

	def remove_all(self):
		for i in range(self.NAtoms):
			self.list.pop()
		self.NAtoms = 0
		self.list = []

	def __str__(self):
		text=''
		for a in self.list:
			text=text+a.__str__()+'\n'
		return(text)

	def read_bonds(self,fileName='default.txt'):
		pairs=[]
		input = open(fileName, 'r')
		for line in input:
			parts = line.split()
			if parts:
				A=int(parts[0])
				B=int(parts[1])
				pairs.append([A,B])
		return pairs
	

	#returns a list of atoms of the type element
	def get_atoms(self, element='C'):
		ats = []
		for at in self.list:
			if at.symbol == element:
				ats.append(at.index)
		return ats
	
	#read xyz procedure; No need to know the number of atoms beforehand
	def read_xyz(self, fileName):
		input = open(fileName, 'r')
		someLine=input.readline()
		parts=someLine.split()
		N = int(parts[0], base = 10)
		someLine = input.readline()
		
		for i in range(N):
			someLine = input.readline()
			parts = someLine.split()
			index = i+1
			symbol = parts[0]
			x = float(parts[1])
			y = float(parts[2])
			z = float(parts[3])
			self.add_atom_symb(symbol,x,y,z)
		input.close()
	
	def read_VASP(self, fileName, Nxyz=[1,1,1]):
		input = open(fileName, 'r')
		#skip first two lines
		someLine = input.readline()
		someLine = input.readline()
		#read lattice vectors
		
		forA1 = input.readline()
		forA2 = input.readline()
		forA3 = input.readline()
		#Currently only for lattice vectors having two zero components
		parts = forA1.split()
		if parts:
			a1 = [float(parts[0]), float(parts[1]), float(parts[2])]
			for i in range(0,3):
				if abs(a1[i])<0.00000001:
					a1[i]=0
			#print(a1)
			
		parts = forA2.split()
		if parts:
			a2 = [float(parts[0]), float(parts[1]), float(parts[2])]
			for i in range(0,3):
				if abs(a2[i])<0.00000001:
					a2[i]=0
			#print a2
			
		parts = forA3.split()
		if parts:
			a3 = [float(parts[0]), float(parts[1]), float(parts[2])]
			for i in range(0,3):
				if abs(a3[i])<0.00000001:
					a3[i]=0
			#print a3
	
		#Read line with atomic species 
		someLine = input.readline()
		atomicSp = someLine.split()
		#print atomicSp
		#Read line with number of atoms
		someLine = input.readline()
		number = someLine.split()
		nb_atoms = []
		tot = 0
		for i in range(0, len(number)):
			sp = int(number[i])
			nb_atoms.append(sp)
			tot = tot+sp
		print(nb_atoms)
		print(tot)
		
		#Skip 1 line
		someLine = input.readline()
		ktot = nb_atoms[0]
		#print tot
		k = 0
		for i in range(0, tot):
			someLine = input.readline()
			parts = someLine.split()
			if parts:
				#print(i+1)
				ax = float(parts[0])
				ay = float(parts[1])
				az = float(parts[2])
				
				x = ax*a1[0]+ay*a2[0]+az*a3[0]
				y = ax*a1[1]+ay*a2[1]+az*a3[1]
				z = ax*a1[2]+ay*a2[2]+az*a3[2]
				if i >= ktot:
					ktot = ktot+nb_atoms[k+1]
					k = k+1
				symbol = atomicSp[k]
				for j in range(0, Nxyz[0]):
					for l in range(0, Nxyz[1]):
						for m in range(0, Nxyz[2]):
							#print(i+1)
							self.add_atom_symb(atomicSp[k], x, y, z)
							self.list[self.NAtoms-1].move_along(j*a1[0]+l*a2[0]+m*a3[0],j*a1[1]+l*a2[1]+m*a3[1],j*a1[2]+l*a2[2]+m*a3[2])


					#self.list[i].reset_symb(i+1,symbol,x,y,z)
				#print symbol
		#self.list.pop()
		return [a1, a2, a3]
	
	def write_to_file(self, fileName):
		output = open(fileName, 'w')
		nb='%d' % self.NAtoms
		nb=nb+'\n\n'
		output.write(nb)
		for a in self.list:
			#output.write(a.__str__())
			output.write(a.str_no_index())
		output.close()

	#add a specific atom, defined by the Z number, to your molecule/system
	def add_atom(self, Znb, x, y, z, const):
		self.list.append(Atoms.Atom())
		self.NAtoms = self.NAtoms + 1 
		i = self.NAtoms
		self.list[i-1].reset(i,Znb,x,y,z)
		self.list[i-1].mul(const)

	#add a specific atom, defined by its atomic symbol, to your molecule/system
	def add_atom_symb(self, symbol, x, y, z):
		self.list.append(Atoms.Atom())
		self.NAtoms = self.NAtoms + 1 
		i = self.NAtoms
		self.list[i-1].reset_symb(i,symbol,x,y,z)
		#self.list[i-1].mul(const)

	#Sorts according to Znb, from small to large	
	def sort_atoms(self):
		a = Atoms.Atom()
		b = Atoms.Atom()
		for i in range(self.NAtoms-1):
			for j in range(i+1, self.NAtoms):
				if (self.list[i].Znb > self.list[j].Znb):
					a.copy(self.list[i])
					b.copy(self.list[j])
					self.list[i].copy(b)
					self.list[j].copy(a)
		for i in range(self.NAtoms):
			#print self.list[i].index
			self.list[i].reset(i+1, self.list[i].Znb, self.list[i].x, self.list[i].y, self.list[i].z)

	#sorts according to Znb, from large to small			
	def reverse_atoms(self):
		a = Atoms.Atom()
		b = Atoms.Atom()
		for i in range(self.NAtoms-1):
			for j in range(i+1, self.NAtoms):
				if (self.list[i].Znb < self.list[j].Znb):
					a.copy(self.list[i])
					b.copy(self.list[j])
					self.list[i].copy(b)
					self.list[j].copy(a)
		for i in range(self.NAtoms):
			self.list[i].index=i+1
			#self.list[i].reset_all(i+1, self.list[i].Znb, self.list[i].x, self.list[i].y, self.list[i].z)

	def sortToFirst(self):
		a = Atoms.Atom()
		b = Atoms.Atom()
		for i in range(1, self.NAtoms-1):
			for j in range(i+1, self.NAtoms):
				if (self.list[i].distance(self.list[0]) > self.list[j].distance(self.list[0])):
					a.copy(self.list[i])
					b.copy(self.list[j])
					self.list[i].copy(b)
					self.list[j].copy(a)
		#Reset atom index according to the new order of the atoms
		for i in range(self.NAtoms):
			self.list[i].index=i+1


	def sort_x(self):
		a = Atoms.Atom()
		b = Atoms.Atom()
		for i in range(self.NAtoms-1):
			for j in range(i+1, self.NAtoms):
				if (self.list[i].x > self.list[j].x):
					a.copy(self.list[i])
					b.copy(self.list[j])
					self.list[i].copy(b)
					self.list[j].copy(a)
		for i in range(self.NAtoms):
			#print self.list[i].index
			self.list[i].reset(i+1, self.list[i].Znb, self.list[i].x, self.list[i].y, self.list[i].z)

	def sort_z(self):
		a = Atoms.Atom()
		b = Atoms.Atom()
		for i in range(self.NAtoms-1):
			for j in range(i+1, self.NAtoms):
				if (self.list[i].z > self.list[j].z):
					a.copy(self.list[i])
					b.copy(self.list[j])
					self.list[i].copy(b)
					self.list[j].copy(a)
		#for i in range(self.NAtoms):
			#print self.list[i].index
		#	self.list[i].reset(i+1, self.list[i].Znb, self.list[i].x, self.list[i].y, self.list[i].z)

			
	def translate_along(self, vx, vy, vz):
		for a in self.list:
			newx = a.x+vx
			newy = a.y+vy
			newz = a.z+vz
			a.reset_xyz(newx,newy,newz)

	def translate_list(self, atoms=[], vx=0, vy=0, vz=0):
		for index in atoms:
			i=index-1
			newx = self.list[i].x+vx
			newy = self.list[i].y+vy
			newz = self.list[i].z+vz
			self.list[i].x=newx
			self.list[i].y=newy
			self.list[i].z=newz
			#self.list[i].reset_xyz(newx,newy,newz)

	#calculate the x,y,z coordinates of the center of tje molecule molecule is centered
	def calculate_CM(self,limit=0):
		Xcm=0.0
		Ycm=0.0
		Zcm=0.0
		if limit==0:
			limit=self.NAtoms
		T=limit
		for a in self.list:
			if a.index<=limit:
				Xcm=Xcm+a.x
				Ycm=Ycm+a.y
				Zcm=Zcm+a.z
		Xcm=Xcm/T
		Ycm=Ycm/T
		Zcm=Zcm/T
		return [Xcm,Ycm,Zcm]



	def multiply_by(self, vx=1.0, vy=1.0, vz=1.0):
		for a in self.list:
			newx = a.x*vx
			newy = a.y*vy
			newz = a.z*vz
			a.reset_xyz(newx,newy,newz)

	def multiply_by_list(self, atoms=[], vx=1.0, vy=1.0, vz=1.0):
		for index in atoms:
			i=index-1
			newx = self.list[i].x*vx
			newy = self.list[i].y*vy
			newz = self.list[i].z*vz
			self.list[i].reset_xyz(newx,newy,newz)


			
	def switch_xz(self):
		for a in self.list:
			temp = a.z
			a.z = a.x
			a.x = temp
	
	def rotate_by(self, angle=0.0, axis='z'):
		if axis=='x':
			self.rotate_aroundX(angle=angle)
		elif axis=='y':
			self.rotate_aroundY(angle=angle)
		else:
			self.rotate_aroundZ(angle=angle)
 
	#rotate around the Z axis
	def rotate_aroundZ(self, angle):
		for a in self.list:
			xrot = a.x * math.cos(angle) + a.y * math.sin(angle) 
			yrot = a.y * math.cos (angle) - a.x * math.sin(angle)
			zrot = a.z
			a.reset_xyz(xrot, yrot, zrot)
	
	#rotate around the X axis
	def rotate_aroundX(self, angle):
		for a in self.list:
			yrot = a.y * math.cos(angle) + a.z * math.sin(angle) 
			zrot = a.z * math.cos (angle) - a.y * math.sin(angle)
			xrot = a.x
			a.reset_xyz(xrot, yrot, zrot)
	
	#rotate around the Y axis
	def rotate_aroundY(self, angle):
		for a in self.list:
			zrot = a.z * math.cos(angle) + a.x * math.sin(angle) 
			xrot = a.x * math.cos (angle) - a.z * math.sin(angle)
			yrot = a.y
			a.reset_xyz(xrot, yrot, zrot)


	#returns an array of the different atomic species in the molecule
	def species(self):
		atom_species=[]
		for a in self.list:
			if a.symbol not in atom_species:
				atom_species.append(a.symbol)
		return atom_species

	#returns the number of atoms of the same species in the molecule
	def how_many_symbol(self, s='Unkown'):
		tot=0
		for a in self.list:
			if a.symbol==s:
				tot=tot+1
		return tot
  
	
	#Curently works only for 2D molecules
	def write_ps_noLengths(self, postScriptName='default.ps',molecule='unknown',pairs=[],step=25):
		psFile = open(postScriptName,'w')
		comment ='%!PS-Adobe-2.0 EPSF-2.0\n%%For: '+molecule
		comment = comment+'\n%%By: Iulia\n'
		psFile.write(comment)
		psFile.write('\n0.2 0.2 0.2 setrgbcolor\n')
		
		for p in pairs:
			A=self.list[p[0]-1]
			B=self.list[p[1]-1]
			#print A
			#print B
			#print A.bond_for_ps(B)
			psFile.write(A.bond(B))
			
		for a in self.list:
			psFile.write(a.str_for_ps(step))

		Xaxis = '\n0 setgray\n1.0 setlinewidth\n10.0 10.0 moveto\n312.5 10.0 lineto\nstroke\n\n'
		Yaxis = '0.5 setgray\n1.0 setlinewidth\n10.0 10.0 moveto\n10.0 312.5 lineto\nstroke\n'
		psFile.write(Xaxis)
		psFile.write(Yaxis)
		showPage='showpage\n'
		psFile.write(showPage)
		psFile.close()


	def write_ps(self, postScriptName='default.ps',molecule='unknown',pairs=[],step=25, scale=1.0):
		psFile = open(postScriptName,'w')
		comment ='%!PS-Adobe-2.0 EPSF-2.0\n%%For: '+molecule
		comment = comment+'\n%%By: Iulia\n'
		psFile.write(comment)
		psFile.write('\n0.2 0.2 0.2 setrgbcolor\n')
		psFile.write('/Helvetica findfont\n8 scalefont\nsetfont\n')

		
		for p in pairs:
			A=self.list[p[0]-1]
			B=self.list[p[1]-1]
			#print A
			#print B
			#print A.bond_for_ps(B)
			#psFile.write(A.bond(B)) #use bond_for_ps(B) for bondlengths
			psFile.write(A.bond_for_ps(B, step))
			
		for a in self.list:
			#if (abs(a.z - 12.5) < 2):
			#print(a)
			psFile.write(a.str_for_ps(step,scale))

		Xaxis = '\n0 setgray\n1.0 setlinewidth\n10.0 10.0 moveto\n312.5 10.0 lineto\nstroke\n\n'
		Yaxis = '0.5 setgray\n1.0 setlinewidth\n10.0 10.0 moveto\n10.0 312.5 lineto\nstroke\n'
		psFile.write(Xaxis)
		psFile.write(Yaxis)
		showPage='showpage\n'
		psFile.write(showPage)
		psFile.close()

	def write_ps_smart(self, psName='default.ps', molecule='Unknown', step=5, threshold=1.97):
		Bonds=[]
		#threshold=1.5
		N=self.NAtoms
		for i in range(N-1):
			for j in range(i+1, N):
				A=self.list[i]
				B=self.list[j]
				magnitude=math.sqrt((A.x-B.x)**2 + (A.y-B.y)**2 + (A.z-B.z)**2)
				if magnitude<threshold:
					Bonds.append([i+1,j+1])
		#print(Bonds)
		self.write_ps(psName, molecule, Bonds, step)
		

	def blender_draw(self, scaling_factor=0.05, pairs=[], bond_size=0.1, diffCol=[], Axes=[]):	
		atom_types=[]

		#Delete any previously existing objects
		for obj in bpy.data.objects:
			obj.select=True
		bpy.ops.object.delete()

	
		#Create group called Molecule and add all atoms and bonds to it
		#bpy.ops.group.create(name='Molecule')
		
		MatName='special'
		bpy.data.materials.new(name=MatName)
		bpy.data.materials[MatName].diffuse_color = [0.6,0.0,0.0]#[0.184,0.553,0.6]
		bpy.data.materials[MatName].specular_intensity = 0.2

		#Create atoms and add materials to them
		for a in self.list:
			x=a.x
			y=a.y
			z=a.z
			s=a.get_radius()*scaling_factor
			bpy.ops.mesh.primitive_uv_sphere_add(size=s, location=(x,y,z))
			bpy.ops.object.shade_smooth()
			if a.index in diffCol:
				bpy.context.active_object.data.materials.append(bpy.data.materials[MatName])
			else: 	
				if a.symbol not in atom_types:
					atom_types.append(a.symbol)
					bpy.data.materials.new(name=a.symbol)
					bpy.data.materials[a.symbol].diffuse_color = a.get_color_blender()
					bpy.data.materials[a.symbol].specular_intensity = 0.2
				bpy.context.active_object.data.materials.append(bpy.data.materials[a.symbol])

		#define bond color (grey)	
		bond_color=[0.302,0.302,0.302]
		bpy.data.materials.new(name='bond')
		bpy.data.materials['bond'].diffuse_color = bond_color
		bpy.data.materials['bond'].specular_intensity = 0.2
		
		#Create bonds 	
		for p in pairs:
			A=self.list[p[0]-1]
			B=self.list[p[1]-1]
			center=[0.5*(A.x+B.x), 0.5*(A.y+B.y), 0.5*(A.z+B.z)]
			magnitude=math.sqrt((A.x-B.x)**2 + (A.y-B.y)**2 + (A.z-B.z)**2)
			dx=B.x-A.x
			dy=B.y-A.y
			dz=B.z-A.z
			angle = math.acos(dz/magnitude)
			bpy.ops.mesh.primitive_cylinder_add(radius=bond_size, depth=magnitude, location=center)
			bpy.ops.object.shade_smooth()
			bpy.ops.transform.rotate(value=(angle), axis=((-1)*dy, dx, 0.0))
			bpy.context.active_object.data.materials.append(bpy.data.materials['bond'])
		
		a_colors=[[0.75,0.1,0.1],[0.1,0.75,0.1],[0.1,0.1,0.75]]
		index_color=0
		for ax in Axes:
			center=[0.5*ax[0], 0.5*ax[1], 0.5*ax[2]]
			magnitude=math.sqrt((ax[0])**2 + (ax[1])**2 + (ax[2])**2)
			#dx=B.x-A.x
			#dy=B.y-A.y
			#dz=B.z-A.z
			dx=ax[0]
			dy=ax[1]
			dz=ax[2]
			angle = math.acos(dz/magnitude)
			bpy.ops.mesh.primitive_cylinder_add(radius=bond_size, depth=magnitude, location=center)
			bpy.ops.object.shade_smooth()
			bpy.ops.transform.rotate(value=(angle), axis=((-1)*dy, dx, 0.0))
			mat_name='axis_%d' % index_color
			bpy.data.materials.new(name=mat_name)
			bpy.data.materials[mat_name].diffuse_color = a_colors[index_color]
			bpy.data.materials[mat_name].specular_intensity = 0.2
			bpy.context.active_object.data.materials.append(bpy.data.materials[mat_name])
			index_color=index_color+1
			
		#add all atoms and bonds to the Molecule group
		#for obj in bpy.data.objects:
		#	bpy.context.scene.objects[obj.name].select = True
		
		#bpy.ops.object.join()

	#scale allows you to change units. Use scale=1.0 for Angstrom	
	def blender_draw_smart(self, scaling_factor=0.05, bond_size=0.1, diffCol=[], threshold=1.9, scale=1.0):	
		atom_types=[]
		#threshold=1.9
		#Delete any previously existing objects
		for obj in bpy.data.objects:
			obj.select=True
		bpy.ops.object.delete()

	
		#Create group called Molecule and add all atoms and bonds to it
		#bpy.ops.group.create(name='Molecule')
		
		MatName='special'
		bpy.data.materials.new(name=MatName)
		bpy.data.materials[MatName].diffuse_color = [0.6,0.0,0.0]#[0.184,0.553,0.6]
		bpy.data.materials[MatName].specular_intensity = 0.2

		#Create atoms and add materials to them
		for a in self.list:
			x=scale*a.x
			y=scale*a.y
			z=scale*a.z
			s=a.get_radius()*scaling_factor
			bpy.ops.mesh.primitive_uv_sphere_add(size=s, location=(x,y,z))
			bpy.ops.object.shade_smooth()
			if a.index in diffCol:
				bpy.context.active_object.data.materials.append(bpy.data.materials[MatName])
			else: 	
				if a.symbol not in atom_types:
					atom_types.append(a.symbol)
					bpy.data.materials.new(name=a.symbol)
					bpy.data.materials[a.symbol].diffuse_color = a.get_color_blender()
					bpy.data.materials[a.symbol].specular_intensity = 0.2
				bpy.context.active_object.data.materials.append(bpy.data.materials[a.symbol])

		#define bond color (grey)	
		bond_color=[0.302,0.302,0.302]
		bpy.data.materials.new(name='bond')
		bpy.data.materials['bond'].diffuse_color = bond_color
		bpy.data.materials[a.symbol].specular_intensity = 0.2
		print(self.NAtoms)
		#Create bonds 	
		for i in range(self.NAtoms-1):
			for j in range(i+1, self.NAtoms):
				A=self.list[i]
				B=self.list[j]
				magnitude=math.sqrt((A.x-B.x)**2 + (A.y-B.y)**2 + (A.z-B.z)**2)
				if magnitude<threshold:
					center=[scale*0.5*(A.x+B.x), scale*0.5*(A.y+B.y), scale*0.5*(A.z+B.z)]
					txt=[0.5*(A.x+B.x), 0.5*(A.y+B.y), 0.5*(A.z+B.z)+1.0]
					dx=B.x-A.x
					dy=B.y-A.y
					dz=B.z-A.z
					angle = math.acos(dz/magnitude)
					bpy.ops.mesh.primitive_cylinder_add(radius=bond_size, depth=scale*magnitude, location=center)
					bpy.ops.object.shade_smooth()
					bpy.ops.transform.rotate(value=(angle), axis=((-1)*dy, dx, 0.0))
					bpy.context.active_object.data.materials.append(bpy.data.materials['bond'])

					#add bond length as text
					#text='%6.3f' % magnitude
						
					#bpy.ops.object.text_add(location=txt)
					#ob=bpy.context.object
					#ob.data.body = text
					#ob.scale*=0.5
					#bpy.context.active_object.data.materials.append(bpy.data.materials['special'])


if __name__ == '__main__':
	Molec = Geometry()
