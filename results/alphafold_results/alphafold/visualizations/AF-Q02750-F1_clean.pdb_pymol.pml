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

load "data/AF-Q02750-F1_clean.pdb", protein
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

load "data/AF-Q02750-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [1428,1430,1427,1429,1432,1433,553,554,570,1280,1282,1272,1281,1438,1446,1482,910,911,1488,1490,1512,1480,1416,1442,1444,1421,1426,729,730,800,801,533,921,924,815,925,706,727,811,404,423,437,405,416,417,395,397,434,1292,1296,1311,1312,1553,1333,1325,1323,401,403,396,399,1531,408,411,1539,1541,1542,974,1349,966,997,998,999,1348,949,965,945,389,390,391,393,394,436,1324] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.329,0.278,0.702]
select surf_pocket2, protein and id [1111,1112,1113,1114,2580,2566,2708,2705,2706,2707,1385,1395,1398,1109,1110,1139,1402,1134,1387,1135,2537,2535,157,130,131,155,774,785,1391,1393,788,1164] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.698,0.361,0.902]
select surf_pocket3, protein and id [1393,1394,1411,219,931,933,1363,215,183,1378,1381,182,184,181,210,939,830,519,245,246,935,938,1360,958,959,942] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.702,0.278,0.639]
select surf_pocket4, protein and id [832,835,294,299,286,290,737,739,740,816,818,823,284,675,839,712,714,677,701,703,707,708] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.545]
select surf_pocket5, protein and id [561,562,889,545,317,318,867,872,876,358,339,463,541,338,445,447,361,363] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.353,0.278]
select surf_pocket6, protein and id [784,786,110,111,769,123,1160,1185,1187,1184,100,102,106,107,1231,1229,1191] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.729,0.361]
select surf_pocket7, protein and id [44,263,261,72,230,255,735,254,252,253,742,762] 
set surface_color,  pcol7, surf_pocket7 
   

deselect

orient
