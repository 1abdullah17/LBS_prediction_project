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

load "data/AF-P07949-F1_clean.pdb", protein
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

load "data/AF-P07949-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [5706,5178,5199,5331,5346,5707,5176,5138,5196,5198,5588,6161,6179,6180,6254,6256,5480,5514,5348,5174,5347,6260,6262,6263,5479,5759,5762,5764,6150,6158,5137,5139,5144,5179,5175,5177,5732,5756,5733,5141,5142,5712,5729] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.302,0.278,0.702]
select surf_pocket2, protein and id [60,1196,1198,1318,1154,1156,1334,1896,1319,1320,1164,1165,89,1304,1321,1325,206,1329,1345,209,233,80,82,232,201,210,86,90,108,1359,1880,1886] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.631,0.361,0.902]
select surf_pocket3, protein and id [670,695,791,376,488,669,626,124,377,506] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.678,0.278,0.702]
select surf_pocket4, protein and id [133,135,1058,112,114,1381,1371,107,115,116,117,1128,1097,1358,1129,1160] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.682]
select surf_pocket5, protein and id [228,231,1297,1293,1296,226,234,1605,1607,255,1606,235,239,309,251,253,254,211,212,233,215,217,1326,222,310,551,543,548,214] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.341]
select surf_pocket6, protein and id [4663,4666,4536,4659,4540,4543,4668,4671,4541,4542,4549,4101,4098,4099,4100,4156,4176,4095,4096,4075,4125,4140,4161,4109,4110,4111,4112] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.522,0.361]
select surf_pocket7, protein and id [6326,5163,6127,6152,6153,6164,5175,6261,6262,6263,6264,6285,6286,6287,5348,5162,5164,5169,5170,5173,5450,5446,6277,6278] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.596,0.278]
select surf_pocket8, protein and id [399,985,988,784,98,99,969,413,416,798,800,803,983,401,393,402,781,782,783] 
set surface_color,  pcol8, surf_pocket8 
   

deselect

orient
