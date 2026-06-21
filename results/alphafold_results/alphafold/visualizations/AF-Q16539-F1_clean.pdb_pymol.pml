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

load "data/AF-Q16539-F1_clean.pdb", protein
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

load "data/AF-Q16539-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [874,876,1264,1344,1345,891,258,895,674,404,406,562,596,597,598,851,837,405,256,260,261,262,269,284,249,294,297,254,250,251,252,296,389,852,390,401,245,247,869,871,875,882,1341,1346,1347,1351,1353,1358,1352,1354,675,276,1375,559,560,586,1243,1252] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.302,0.278,0.702]
select surf_pocket2, protein and id [727,53,54,9,35,36,33,24,704,188,192,348,171,702,703,2741,709,715,723,728,2701,2702] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.631,0.361,0.902]
select surf_pocket3, protein and id [677,650,1329,678,679,683,685,2746,860,861,2748,696,2780,666,635,639,622,627,628,629,631,646] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.678,0.278,0.702]
select surf_pocket4, protein and id [1941,1943,1525,2014,1536,1537,2257,1519,1476,1470,1471,1473,1475,1472,1506,1805,2259] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.682]
select surf_pocket5, protein and id [1303,1007,1008,1285,926,956,958,927,1277,925,975,1006,1004,1005] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.341]
select surf_pocket6, protein and id [1070,1069,1034,1035,2413,1072,1095,653,660,661,1099,2459,1060,2432,2410,2411] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.522,0.361]
select surf_pocket7, protein and id [1575,2316,2318,2317,2326,1576,1144,1148,1116,1151,1154,2321,1113,1117,1119,2323] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.596,0.278]
select surf_pocket8, protein and id [1237,1239,1242,1244,909,910,911,942,933,1238,1449,1451,1439,1411,1413,1414,1229,1671,1672,1673,1718] 
set surface_color,  pcol8, surf_pocket8 
   

deselect

orient
