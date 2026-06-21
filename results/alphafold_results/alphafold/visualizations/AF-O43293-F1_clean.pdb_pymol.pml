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

load "data/AF-O43293-F1_clean.pdb", protein
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

load "data/AF-O43293-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [444,445,1332,1333,192,326,327,1309,1310,1311,1173,208,181,172,177,185,187,178,165,166,169,189,193,346,1130,1131,1132,1450,1451,1452,1459,1149,479,2488,1444,1446,2450,475,477,478,480,2446,2451,2448,375,2486,517,1324,1325,1326,1327] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.416,0.278,0.702]
select surf_pocket2, protein and id [207,1303,1188,1189,153,205,208,150,325,326,327,1307,1309,1311,1173,156,159,161,163,151,803,805,1161,312,1302,628,764,780,783,1304,1305,758,1299,626,781] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.902,0.361,0.878]
select surf_pocket3, protein and id [2435,2437,2439,1813,1814,2438,2440,1432,1436,1440,1442,1517,1519,1523,1524,1492,1521,1529,1437,1527,1528,1465,1776,1779,1490,1493,1461,1452] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.702,0.278,0.380]
select surf_pocket4, protein and id [1062,1063,2138,2140,2142,2143,2144,2150,2155,1057,1058,1059,1065,1581,1584,1031,1586,1580,1582,1587,1593,1591] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.620,0.361]
select surf_pocket5, protein and id [818,1153,2378,2402,806] 
set surface_color,  pcol5, surf_pocket5 
   

deselect

orient
