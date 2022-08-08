import os
import sys
from collections import defaultdict
from Bio.PDB import *
import numpy
import math

def filter(wd):

	pdbs = os.listdir(wd) #Lists the files you want to use (PDBs)

	os.system("mkdir /home/douglas_lima/pdb/unzipped/RMN_rmk/")
	os.system("mkdir /home/douglas_lima/pdb/unzipped/Carbohydrates_rmk/")

	os.chdir(wd)#Changes to that directory


	#Data structures used to saved selected atoms and molecules:

	save_resolution = [] #Saves the different values of resolution in the pdbs
	tags = [] #Saves the tags of the selected sugars
	hetatms = [] #Saves the Heteroatoms that are present in the sugar pdb file
	wholetags = [] #Used to count number of occurences of sugars
	names = [] # Saves the names of the selected sugars
	whole_resolution = [] #Used to count number of occurences of a given resolution
	pdb_tag = defaultdict(list) # Dictionary {Name of the pdb : List of tags of sugars present in the pdb}
	owab = {} # Dictionary to store owab values
	b_factors = {}
	total_atoms = {}
	output = open("OWAB.txt","w")


	content = open("/home/douglas_lima/pdb/unzipped/Sugar_in_pdbs_rmk.txt", "w")



	for name in pdbs:

		#print name

		arq = open(name, "r")

		hetnams = []
		sugars = []
		lines = []
		res = []
		owab[name] = 0
		total_atoms[name] = 0
		b_factors[name] = 0

		for line in arq.readlines():
			lines.append(line)
		
		# Calculo OWAB
		for linha in lines:
			if linha[0:4] == "ATOM" or linha[0:6] == "HETATM":
				owab[name] = owab[name] + (float(linha[54:59])*float(linha[60:65]))
				total_atoms[name] = total_atoms[name] + 1 # Itera sobre todos atomos somando ocupancia
				output.write(name + "      " + str(round(owab[name]/total_atoms[name],2)) +"\n")	#59-61 dentro do for?
				b_factors[name] = round(owab[name]/total_atoms[name],2)
				print(b_factors[name])

		arq.close()

		
		for i in lines:
			x = i.split()
			if x == []:
				pass
			if x[0] == "REMARK" and x[1] == "2":
				for b in range(len(x)):
					if "RESOLUTION" in x[b]:
						if x[3] == "NOT" or x[3] == "NULL":
							print(name)
							os.system("mv " + name + " /home/douglas_lima/pdb/unzipped/RMN_rmk/")
						else:
							if float(x[3]) <= 2.0 and b_factors[name] <= 60.0:
								res.append(x[3])
								resolution = x[3]
							else: 
								resolution = x[3] 
			if x[0] == "HETNAM":
				if float(resolution) <= 2.0:
					
					for b in range(len(x)):
						if "SIALIC" in x[b]:
							print(name, resolution) 
							sugars.append(x)
						if "NEURAMINIC" in x[b]: 
							print(name, resolution)
							sugars.append(x)




				hetnams.append(x)

		for x in range(len(sugars)):
			y = sugars[x]
			if y[1] == "2" or y[1] == "3" or y[1] == "4" or y[1] == "5" or y[1] == "6":
				tag = y[2]
				if tag not in tags:
					tags.append(tag)
					names.append(" ".join(y[3:]))

			else:
				tag = y[1]
				if tag not in tags:
					tags.append(tag)
					names.append(" ".join(y[2:]))

			y = "-".join(sugars[x][1:])
			content.write(name + "	" + y)
			content.write("\n")
			wholetags.append(tag)
			pdb_tag[name].append(tag)

		for y in hetnams:
			if y not in hetatms:
				hetatms.append(y)

		for x in range(len(res)):
			whole_resolution.append(res[x])
			if res[x] not in save_resolution:
				save_resolution.append(res[x])


		if len(sugars) != 0:
			os.chdir("/home/douglas_lima/pdb/unzipped/Carbohydrates_rmk/")
			file = open(name + ".carbo" + ".pdb", "w")

			for i in lines:
				x = i.split()
				for j in tags:
					if x == []:
						pass
					if x[0][:6] == "HETATM":
						if i[17:20] == j:
							file.write(i)

			file.close()
		os.chdir(wd)

	content.close()	



	os.chdir("/home/douglas_lima/pdb/unzipped/")

	file = open("Resolutions_rmk.txt", "w")
	for i in range (len(save_resolution)):
		file.write(save_resolution[i])
		file.write(" " + str(whole_resolution.count(save_resolution[i])))
		file.write("\n")
	file.close()

	#Use if you want to remove sugar with occurence lower than x
	'''
	remove = []
	for i in range(len(tags)):
		if wholetags.count(tags[i]) <= x:
			remove.append(tags[i])
			print tags[i], wholetags.count(tags[i])	
	'''


	file = open("Heteroatoms_in_Selected_rmk.txt", "w")
	for i in range (len(hetatms)):
		file.write(str(hetatms[i]))
		file.write("\n")
	file.close()


	file = open("Selected_Sugars_rmk.txt", "w")
	for i in range (len(tags)):
		file.write(tags[i])
		file.write("	" + str(names[i]))
		file.write("	" + str(wholetags.count(tags[i])))
		file.write("\n")
	file.close()



##########################################################################################################

def Ring_sep(name):

	ring_atoms = ["O6", "C6", "C2", "C3", "C4", "C5"]	
	
	arq = open(name, "r")

	lines = []

	os.system("mkdir /home/douglas_lima/pdb/unzipped/Rings_rmk/")

	for line in arq:
		if line[0:4] == "ATOM" or line[0:6] == "HETATM":
			lines.append(line)


	res = {}
	resid = []
	chainid = []
	ident = []

	for j in lines:
		aux = j[21:26]
		if aux not in ident:
			ident.append(aux)
			key = j[21:26]
			res[key] = []				
		if j[22:26] not in resid:
			resid.append(j[22:26])
		if j[21:22] not in chainid: 
			chainid.append(j[21:22])

		if j[22:26] == ident[-1][-4:] and j[21:22] == ident[-1][0]:
			res[key].append(j)

	for x in res:
		tag = x.replace(" ", "")
		d = open(name[:7] + ".Ring" + str(tag) + ".pdb", "w")
		for y in range (len(res[x])):
			d.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(res[x][y][0:6], int(res[x][y][6:11]), res[x][y][12:16], res[x][y][16:17], res[x][y][17:20], res[x][y][21:22], int(res[x][y][22:26]), res[x][y][26:27], float(res[x][y][30:38]), float(res[x][y][38:46]), float(res[x][y][46:54]), float(res[x][y][54:60]), float(res[x][y][60:66]), res[x][y][76:78], res[x][y][78:80]))
		d.close()
		os.system("mv " + name[:7] + ".Ring" + str(tag) + ".pdb /home/douglas_lima/pdb/unzipped/Rings_rmk/")

			
		
	return

def order_ring(file1):

	name = file1
	out_name = file1[:-4] + "_ord.pdb"

	order1 = ["O5", "C1", "C2", "C3", "C4", "C5"]
	order2 = ["O", "C1", "C2", "C3", "C4", "C5"]

	arq = open(name, "r")

	lines = []

	for line in arq:
		if line[0:4] == "ATOM" or line[0:6] == "HETATM":
			lines.append(line)
	arq.close()

	arq = open(out_name, "w")

	for i in range(len(lines)):
		if order1[0] == lines[i][12:16].replace(" ", "").replace("'", "").replace("B", "") or order2[0] == lines[i][12:16].replace(" ", "").replace("'", "").replace("B", "").replace("'", "").replace("B", "").replace("'", "").replace("B", ""):
			arq.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(lines[i][0:6], int(lines[i][6:11]), lines[i][12:16], lines[i][16:17], lines[i][17:20], lines[i][21:22], int(lines[i][22:26]), lines[i][26:27], float(lines[i][30:38]), float(lines[i][38:46]), float(lines[i][46:54]), float(lines[i][54:60]), float(lines[i][60:66]), lines[i][76:78], lines[i][78:80]))

	for i in range(len(lines)):
		if order1[1] == lines[i][12:16].replace(" ", "").replace("'", "").replace("B", "") or order2[1] == lines[i][12:16].replace(" ", "").replace("'", "").replace("B", ""):
			arq.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(lines[i][0:6], int(lines[i][6:11]), lines[i][12:16], lines[i][16:17], lines[i][17:20], lines[i][21:22], int(lines[i][22:26]), lines[i][26:27], float(lines[i][30:38]), float(lines[i][38:46]), float(lines[i][46:54]), float(lines[i][54:60]), float(lines[i][60:66]), lines[i][76:78], lines[i][78:80]))

	for i in range(len(lines)):
		if order1[2] == lines[i][12:16].replace(" ", "").replace("'", "").replace("B", "") or order2[2] == lines[i][12:16].replace(" ", "").replace("'", "").replace("B", ""):
			arq.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(lines[i][0:6], int(lines[i][6:11]), lines[i][12:16], lines[i][16:17], lines[i][17:20], lines[i][21:22], int(lines[i][22:26]), lines[i][26:27], float(lines[i][30:38]), float(lines[i][38:46]), float(lines[i][46:54]), float(lines[i][54:60]), float(lines[i][60:66]), lines[i][76:78], lines[i][78:80]))

	for i in range(len(lines)):
		if order1[3] == lines[i][12:16].replace(" ", "").replace("'", "").replace("B", "") or order2[3] == lines[i][12:16].replace(" ", "").replace("'", "").replace("B", ""):
			arq.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(lines[i][0:6], int(lines[i][6:11]), lines[i][12:16], lines[i][16:17], lines[i][17:20], lines[i][21:22], int(lines[i][22:26]), lines[i][26:27], float(lines[i][30:38]), float(lines[i][38:46]), float(lines[i][46:54]), float(lines[i][54:60]), float(lines[i][60:66]), lines[i][76:78], lines[i][78:80]))

	for i in range(len(lines)):
		if order1[4] == lines[i][12:16].replace(" ", "").replace("'", "").replace("B", "") or order2[4] == lines[i][12:16].replace(" ", "").replace("'", "").replace("B", ""):
			arq.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(lines[i][0:6], int(lines[i][6:11]), lines[i][12:16], lines[i][16:17], lines[i][17:20], lines[i][21:22], int(lines[i][22:26]), lines[i][26:27], float(lines[i][30:38]), float(lines[i][38:46]), float(lines[i][46:54]), float(lines[i][54:60]), float(lines[i][60:66]), lines[i][76:78], lines[i][78:80]))

	for i in range(len(lines)):
		if order1[5] == lines[i][12:16].replace(" ", "").replace("'", "").replace("B", "") or order2[5] == lines[i][12:16].replace(" ", "").replace("'", "").replace("B", ""):
			arq.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(lines[i][0:6], int(lines[i][6:11]), lines[i][12:16], lines[i][16:17], lines[i][17:20], lines[i][21:22], int(lines[i][22:26]), lines[i][26:27], float(lines[i][30:38]), float(lines[i][38:46]), float(lines[i][46:54]), float(lines[i][54:60]), float(lines[i][60:66]), lines[i][76:78], lines[i][78:80]))

	arq.close()

	return out_name

def mkpdbdat(file1):

	name = order_ring(file1)
	out_name = name[:-4] + ".dat"

	arq = open(name, "r")

	lines = []

	for line in arq:
		if line[0:4] == "ATOM" or line[0:6] == "HETATM":
			lines.append(line)

	arq.close()

	arq = open(out_name, "w")

	arq.write(out_name + " ")

	for i in range(0,6):
		arq.write(lines[i][30:38] + " " + lines[i][38:46] + " " + lines[i][46:54] + " ")

	arq.close()

	return out_name

##########################################################################

def norm(a):
	return math.sqrt(numpy.sum(a*a))

def tofloat(a):
	b=[]
	for i in a:
		b.append(float(i))
	return b

def calculate_pucker(pdbinput):

	datname = mkpdbdat(pdbinput)

	inputfile = file(datname, 'r')

	atoms=numpy.zeros((6,3),dtype='float64')

	for line in inputfile:
		sline=line.split()
		header=sline[0]
		for i in range(1, 19, 3):
			atoms[(i-1)/3]=tofloat(sline[i:i+3])

		center = numpy.add.reduce(atoms)/6.
		atoms = atoms - center

		r1a = numpy.zeros((3),dtype='float64')
		r2a = numpy.zeros((3),dtype='float64')
		for j,i in enumerate(atoms[0:6]):
			r1a += i * math.sin(2.*math.pi*j/6.)
			r2a += i * math.cos(2.*math.pi*j/6.)

		n = numpy.cross(r1a,r2a)
		n = n / norm(n)

		z = numpy.dot(atoms,n)
		q2cosphi = 0.
		q2sinphi = 0.
		q1cosphi = 0.
		q1sinphi = 0.
		q3 = 0.
		bigQ = 0.
		sqrt_2 = math.sqrt(2.)
		inv_sqrt_6 = math.sqrt(1./6.)
		for j,i in enumerate(z):
			q2cosphi += i*math.cos(2.*math.pi*2.*j/6.)
			q2sinphi -= i*math.sin(2.*math.pi*2.*j/6.)
			q1cosphi += i*math.cos(2.*math.pi*j/6.)
			q1sinphi -= i*math.sin(2.*math.pi*j/6.)
			q3 += i*math.cos(j*math.pi)
			bigQ += i*i
		q2cosphi = sqrt_2 * inv_sqrt_6 * q2cosphi
		q2sinphi = sqrt_2 * inv_sqrt_6 * q2sinphi
		q3 = inv_sqrt_6 * q3
		q2 = math.sqrt(q2cosphi*q2cosphi + q2sinphi*q2sinphi)
		q1 = math.sqrt(q1cosphi*q1cosphi + q1sinphi*q1sinphi)
		bigQ = math.sqrt(bigQ)

		if (q2cosphi > 0.):
			if (q2sinphi > 0.):
				phi = math.degrees(math.atan(q2sinphi/q2cosphi))
			else:
				phi = 360. - abs(math.degrees(math.atan(q2sinphi/q2cosphi)))
		else:
			if (q2sinphi > 0.):
				phi = 180. - abs(math.degrees(math.atan(q2sinphi/q2cosphi)))
			else:
				phi = 180. + abs(math.degrees(math.atan(q2sinphi/q2cosphi)))
				theta = math.degrees(math.atan(q2/q3))
		if (q3 > 0.):
			if (q2 > 0.):
				theta = math.degrees(math.atan(q2/q3))
			else:
				theta = 360. - abs(math.degrees(math.atan(q2/q3)))
		else:
			if (q2 > 0.):
				theta = 180. - abs(math.degrees(math.atan(q2/q3)))
			else:
				theta = 180. + abs(math.degrees(math.atan(q2/q3)))

		#bigQ2 = numpy.array([q1,q2,q3],dtype='float64')
		#bigQ2 = math.sqrt((bigQ2*bigQ2).sum())

		string1 = str('%7.3f' % (phi))
		string2 = str('%7.3f' % (theta))
		string3 = str('%7.3f' % (bigQ))

	inputfile.close()

	return string1, string2, string3

#################################################################################################

def Separate():
	pdbs = os.listdir("/home/douglas_lima/pdb/unzipped/Carbohydrates/remake/")

	os.chdir("/home/douglas_lima/pdb/unzipped/")

	os.system("mkdir Separated_rmk/")

	o_atoms = ["O1A", "O2", "O3", "O4", "O5", "O8", "O9", "O10", "O7", "O1B"]

	os.chdir("/home/douglas_lima/pdb/unzipped/")

	report = open("Sugars_Links_rmk.txt", "w")

	os.chdir("/home/douglas_lima/pdb/unzipped/Carbohydrates/remake/")

	for i in pdbs:
		print(i)
		connected = []
		linked = {}
		keys = []
		disordered = defaultdict(list)
		io = PDBIO()
		parser = PDBParser()
		structure = parser.get_structure(i,i)
		atom_list = Selection.unfold_entities(structure, 'A') # A for atoms
		ns = NeighborSearch(atom_list)
		for x in atom_list:
			for y in o_atoms:
				if y in x.get_name():
					if x.get_parent().is_disordered():
						center = x.get_coord()
						neighbors = ns.search(center, 2) # 2.0 for distance in angstrom
						atoms_list = Selection.unfold_entities(neighbors, 'A') # A for atoms
						for j in atoms_list:
							if j.get_name()[-1] != x.get_name()[-1] and j.get_altloc() == x.get_altloc() and j.get_name()[:-1] == x.get_name()[:-1]: 
								value = (j.get_parent(), x.get_parent())
								value = set(value)
								key = str(j.get_parent().get_id()[1]) + str(x.get_parent().get_id()[1]) + str(j.get_altloc()) + str(j.get_name()[-1])
								key = key.replace(" ", "")
								if str(j.get_name()[-1]) not in disordered[str(x.get_parent().get_id()[1])]: 
									disordered[str(x.get_parent().get_id()[1])].append(str(j.get_name()[-1])) 
								if str(x.get_name()[-1]) not in disordered[str(x.get_parent().get_id()[1])]:
									disordered[str(x.get_parent().get_id()[1])].append(str(x.get_name()[-1]))
								if keys != []:
									for q in range(len(keys)):
										if value.intersection(linked[keys[q]]):
											aux = value.union(linked[keys[q]])
											del linked[keys[q]]
											new_key = key + keys[q]
											del keys[q]
											linked[new_key] = aux
											keys.append(new_key)
										else:
											if keys[q] == keys[-1]:
												keys.append(key)
												linked[key] = value
								else:
									keys.append(key)
									linked[key] = value
								connected.append(j.get_parent())
								connected.append(x.get_parent())

							if j.get_parent() != x.get_parent() and j.get_altloc() == x.get_altloc(): 
								value = (j.get_parent(), x.get_parent())
								value = set(value)
								key = str(j.get_parent().get_id()[1]) + str(x.get_parent().get_id()[1]) + str(j.get_altloc()) + str(j.get_parent().get_parent().get_id())
								key = key.replace(" ", "")
								if str(j.get_name()[-1]) not in disordered[str(x.get_parent().get_id()[1])]: 
									disordered[str(x.get_parent().get_id()[1])].append(str(j.get_name()[-1])) 
								if str(x.get_name()[-1]) not in disordered[str(x.get_parent().get_id()[1])]:
									disordered[str(x.get_parent().get_id()[1])].append(str(x.get_name()[-1]))
								if keys != []:
									for q in range(len(keys)):
										if value.intersection(linked[keys[q]]):
											aux = value.union(linked[keys[q]])
											del linked[keys[q]]
											new_key = key + keys[q]
											del keys[q]
											linked[new_key] = aux
											keys.append(new_key)
										else:
											if keys[q] == keys[-1]:
												keys.append(key)
												linked[key] = value
								else:
									keys.append(key)
									linked[key] = value
								connected.append(j.get_parent())
								connected.append(x.get_parent())

					else:
						center = x.get_coord()
						neighbors = ns.search(center, 2) # 2.0 for distance in angstrom
						atoms_list = Selection.unfold_entities(neighbors, 'A') # A for atoms
						for j in atoms_list:
							if j.get_parent() != x.get_parent() and j.get_parent().get_parent().get_id() == x.get_parent().get_parent().get_id(): 
								value = (j.get_parent(), x.get_parent())
								value = set(value)
								key = str(j.get_parent().get_id()[1]) + str(x.get_parent().get_id()[1]) + str(j.get_parent().get_parent().get_id())
								key = key.replace(" ", "")
								if keys != []:
									for q in range(len(keys)):
										if value.intersection(linked[keys[q]]):
											aux = value.union(linked[keys[q]])
											del linked[keys[q]]
											new_key = key + keys[q]
											del keys[q]
											linked[new_key] = aux
											keys.append(new_key)
										else:
											if keys[q] == keys[-1]:
												keys.append(key)
												linked[key] = value
								else:
									keys.append(key)
									linked[key] = value			
								connected.append(j.get_parent())
								connected.append(x.get_parent())

							if ("'" in j.get_name() and "'" not in x.get_name()) or  ("'" in x.get_name() and "'" not in j.get_name()):
								io.set_structure(j.get_parent())
								os.chdir("/home/douglas_lima/pdb/unzipped/Separated_rmk/")
								io.save(i[:7] + "_2_" + str(j.get_parent().get_parent().get_id()) + ".pdb")
								connected.append(j.get_parent())
								connected.append(x.get_parent())
								os.chdir("/home/douglas_lima/pdb/unzipped/Carbohydrates/remake/")
							

		for p in atom_list:
			if p.get_parent() not in connected:
				os.chdir("/home/douglas_lima/pdb/unzipped/Separated_rmk/")
				if not os.path.exists(i[:7] + "_1_" + str(p.get_parent().get_id()[1]) + str(p.get_parent().get_parent().get_id()) + ".pdb"):
					io.set_structure(p.get_parent())
					io.save(i[:7] + "_1_" + str(p.get_parent().get_id()[1]) + str(p.get_parent().get_parent().get_id()) + ".pdb")
					report.write(i[:7] + "	" + str(1) + "	" + str(p.get_parent().get_resname()) + "\n")
		os.chdir("/home/douglas_lima/pdb/unzipped/Carbohydrates/remake/")



		if len(linked) != 0:
			
			bonded = []

			for r in linked:
				bonded.append(linked[r])

			ids = defaultdict(list)
			


			for o in range(len(bonded)):
				w = list(bonded[o])
				monomers = ""
				for c in w:
					ids[o].append(str(list(c.get_id())[1])) 
					monomers = monomers + "-" + str(c.get_resname())
				if c.is_disordered():
					report.write(i[:7] + "	" + str(len(disordered[str(c.get_id()[1])])) + "	" + monomers + "\n")
				else:
					report.write(i[:7] + "	" + str(len(ids[o])) + "	" + monomers + "\n")
				ids[o].append(str(c.get_parent().get_id()))

			arq = open(i,"r")

			lines = []

			for n in arq.readlines():
				lines.append(n)

			os.chdir("/home/douglas_lima/pdb/unzipped/Separated_rmk/")

			for p in ids:
				if isinstance(ids[p][-1], int):
					file = open(i[:7] + "_" + str(len(ids[p])) + "_" + str(p+1) + ".pdb", "a")
				else:
					file = open(i[:7] + "_" + str(len(ids[p])-1) + "_" + str(ids[p][-1]) + "_" + str(p+1) + ".pdb", "a")

					
				for l in lines:
					if str(l[22:26].replace(" ", "")) in ids[p] and l[21:22] == ids[p][-1]:
						file.write(l)
				file.close()

			os.chdir("/home/douglas_lima/pdb/unzipped/Carbohydrates/remake/")

			arq.close()

	report.close()


##########################################################################################################

wd = "/home/douglas_lima/pdb/unzipped/" #Location of your PDB files to be analysed

filter(wd)

os.chdir(wd)#Changes to that directory

output = open("OWAB.txt", "r")

first_interval = 0
second_interval = 0
third_interval = 0
last_interval = 0

distribution_OWAB = open("distribution_OWAB.txt", "w")

for line in output:
        x = line.split()
        if float(x[1]) <= 19.9:
                first_interval += 1
        elif 20.0 <= float(x[1]) <= 39.9:
                second_interval += 1
        elif 40.0 <= float(x[1]) <= 59.9:
                third_interval += 1
        elif 60.0 <= float(x[1]):
                last_interval += 1

output.close()

distribution_OWAB.write("1.0 to 19.9" + "       " + str(first_interval) + "\n")
distribution_OWAB.write("20.0 to 39.9" + "      " + str(second_interval) + "\n")
distribution_OWAB.write("40.0 to 59.9" + "      " + str(third_interval) + "\n")
distribution_OWAB.write("60.0 or more" + "      " + str(last_interval) + "\n")

distribution_OWAB.close()

'''
Separate()

pdbs = os.listdir("/home/douglas_lima/pdb/unzipped/Separated_rmk/")

os.chdir("/home/douglas_lima/pdb/unzipped/Separated_rmk/")

for i in pdbs:
	Ring_sep(i)

pdbs = os.listdir("/home/douglas_lima/pdb/unzipped/Rings_rmk/")

os.chdir("/home/douglas_lima/pdb/unzipped/Rings_rmk/")

os.system("mkdir /home/douglas_lima/pdb/unzipped/Not_Rings_rmk")
os.system("mkdir /home/douglas_lima/pdb/unzipped/Rings2_rmk")
os.system("mkdir /home/douglas_lima/pdb/unzipped/Puckering_Files_rmk")

output1 = "phi_rmk.txt"
output2 = "theta_rmk.txt"
output3 = "bigQ_rmk.txt"

out1 = open("/home/douglas_lima/pdb/unzipped/Puckering_Files_rmk/" + output1, 'w')
out1.write("""# This file was created by Glados!\n@    title "Puckering"\n@    xaxis  label "Time (ns)"\n@    yaxis  label "Phi angle"\n@TYPE xy\n""")	
out1.close()

out2 = open("/home/douglas_lima/pdb/unzipped/Puckering_Files_rmk/" + output2, 'w')
out2.write("""# This file was created by Glados!\n@    title "Puckering"\n@    xaxis  label "Time (ns)"\n@    yaxis  label "Theta angle"\n@TYPE xy\n""")
out2.close()

out3 = open("/home/douglas_lima/pdb/unzipped/Puckering_Files_rmk/" + output3, 'w')
out3.write("""# This file was created by Glados!\n@    title "Puckering"\n@    xaxis  label "Time (ns)"\n@    yaxis  label "Big Q"\n@TYPE xy\n""")
out3.close()


for name in pdbs:

	arq = open(name, "r")

	for line in arq:
		if line[0:6] == "HETATM" or line[0:4] == "ATOM":
			tag = line[17:20]

	arq.close()
	try:
		string1, string2, string3 = calculate_pucker(name)

		out1 = open("/home/douglas_lima/pdb/unzipped/" + output1, 'a')
		out2 = open("/home/douglas_lima/pdb/unzipped/" + output2, 'a')
		out3 = open("/home/douglas_lima/pdb/unzipped/" + output3, 'a')

		out1.write("." + name[3:-4] + "	" + tag + "	" + string1 + "\n")
		out2.write("." + name[3:-4] + "	" + tag + "	" + string2 + "\n")
		out3.write("." + name[3:-4] + "	" + tag + "	" + string3 + "\n")

		out1.close()
		out2.close()
		out3.close()

	except IndexError:
		os.system("mv " + name + " /home/douglas_lima/pdb/unzipped/Not_Rings_rmk/")




	os.system("mv " + name + " /home/douglas_lima/pdb/unzipped/Rings2_rmk/")
	os.system("mv  *ord* /home/douglas_lima/pdb/unzipped/Puckering_Files_rmk/")
'''
#################################################################################################################
