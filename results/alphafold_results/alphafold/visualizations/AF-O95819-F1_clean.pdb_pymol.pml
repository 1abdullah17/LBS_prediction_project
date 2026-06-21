from pymol import cmd,stored

set depth_cue, 1
set fog_start, 0.4

set_color b_col, [36,36,85]
set_color t_col, [10,10,10]
set bg_rgb_bottom, b_col
set bg_rgb_top, t_col      
set bg_gradient

set  spec_power  =  200
set  spec_refl   =  0

load "data/AF-O95819-F1_clean.pdb", protein
create ligands, protein and organic
select xlig, protein and organic
delete xlig

hide everything, all

color white, elem c
color bluewhite, protein
#show_as cartoon, protein
show surface, protein
#set transparency, 0.15

show sticks, ligands
set stick_color, magenta




# SAS points

load "data/AF-O95819-F1_clean.pdb_points.pdb.gz", points
hide nonbonded, points
show nb_spheres, points
set sphere_scale, 0.2, points
cmd.spectrum("b", "green_red", selection="points", minimum=0, maximum=0.7)


stored.list=[]
cmd.iterate("(resn STP)","stored.list.append(resi)")    # read info about residues STP
lastSTP=stored.list[-1] # get the index of the last residue
hide lines, resn STP

cmd.select("rest", "resn STP and resi 0")

for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.4","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))



set_color pcol1 = [0.361,0.576,0.902]
select surf_pocket1, protein and id [207,208,153,1273,750,324,1265,307,755,771,1188,767,150,769,770,152,775,779,780,1263,1174,1272,1274,151,155,790,793,787,786,818,1266,1162,1189,166,167,160,162,156,158] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.310,0.702]
select surf_pocket2, protein and id [4383,4384,4385,4834,4378,4390,4380,4382,5231,4825,6279,6281,6283,6293,6295,6302,6303,4091,6288,4859,4836,4837,4842,4851,4856,4857,4858,5257,5247,4399,4400,4401,5966,5954,5960,5962,5965,5608,5614,5950,5599,5606,5610,5240,5241,5602,5613,6277,4082,4084,4077,4372,5945,5947,5941] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.498,0.361,0.902]
select surf_pocket3, protein and id [3968,3969,6371,6184,6185,6186,6195,6197,6083,6214,6216,6202,6207,6210,3987,6451,6454,6449] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.525,0.278,0.702]
select surf_pocket4, protein and id [4129,4136,4214,4648,4118,4649,4655,4657,4290,4292,4113,4115,4117,4624,4620,4619,4621,4109] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.851,0.361,0.902]
select surf_pocket5, protein and id [592,593,594,595,288,290,596,250,251,253,254,104,106,618,620,94,105] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.600]
select surf_pocket6, protein and id [6293,6295,6302,4091,6288,4390,4093,4095,4097,4100,4101,4399,4400,4401,4402,4406] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.596]
select surf_pocket7, protein and id [1376,1685,1649,4921,1654,4907,4905,4906,4919,1374,5064,5061,5063] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.325]
select surf_pocket8, protein and id [4605,4641,4313,4314,4315,4342,4345,4604,4562,4753,4757,4759] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.478,0.361]
select surf_pocket9, protein and id [4552,4766,4577,4772,5052,5097,5104,5106,5109,4553,4781,5116] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.510,0.278]
select surf_pocket10, protein and id [6489,6492,6513,6515,5996,5998,6014,6015,6167,6483,6482,6484,3949,6488,3951,6100,6191,6192,6172,6516,3952] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.831,0.361]
select surf_pocket11, protein and id [4871,5262,5273,4889,5638,5632,5633,5277,5672,5673,5272,5643] 
set surface_color,  pcol11, surf_pocket11 
   

deselect

orient
