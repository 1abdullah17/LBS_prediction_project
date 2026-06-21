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

load "data/AF-Q13464-F1_clean.pdb", protein
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

load "data/AF-Q13464-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [853,854,1774,1770,1663,1675,1767,1765,679,681,688,689,675,677,727,724,669,674,2981,2983,2985,1691,1690,1325,1664,726,672,837,1290,1295,1272,1273] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.310,0.702]
select surf_pocket2, protein and id [975,977,1773,1788,1789,1787,1636,1790,1795,704,706,692,695,697,871,956,971,1915,1913,707,854,1652,1772,1774,1663,1666,1675,1651,679,688,689,727,710,711,712,1917,1918,1919,1920,1927] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.498,0.361,0.902]
select surf_pocket3, protein and id [9021,9635,9003,9004,9007,9018,9627,9812,9001,9636,9637,9638,9801,9648,9066,9646,9649,9022,9025,9878,9005,9006,9830,9203,9863,9828,9862,9843,9844,9845,9187,9188,9880,9881,9882] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.525,0.278,0.702]
select surf_pocket4, protein and id [9775,9468,9475,9656,9766,9776,9778,9781,9782,9783,9637,9638,9639,9653,9655,9792,9802,9803,9353,9788,9804] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.851,0.361,0.902]
select surf_pocket5, protein and id [1347,1349,1350,1352,1354,1388,1355,1698,2886,1313,2896,2898,2900,2902,2972,2974,1377] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.600]
select surf_pocket6, protein and id [3086,3090,3167,342,3158,314,341,3087,340,1204,1203,1207,298,299,300,302] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.596]
select surf_pocket7, protein and id [6662,6664,6660,2357,2355,2301,2318,2322,2327,2328,2302,2023,1996,1995,2324,2330] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.325]
select surf_pocket8, protein and id [1340,1658,1667,1954,2196,2238,1952,1940,1942,1339,1364,1666] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.478,0.361]
select surf_pocket9, protein and id [8345,8375,9257,8934,8937,8936,8952,8787] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.510,0.278]
select surf_pocket10, protein and id [717,718,720,877,3063,3061,903,875,699,701,708,709,1230,3066,863,1227,719] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.831,0.361]
select surf_pocket11, protein and id [1163,3252,582,3243,3245,3247,3239,3241,3242,560,224] 
set surface_color,  pcol11, surf_pocket11 
   

deselect

orient
