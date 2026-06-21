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

load "data/AF-P00533-F1_clean.pdb", protein
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

load "data/AF-P00533-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [5338,5202,6116,6117,5688,5692,5150,5151,5152,5154,5155,5157,5715,5718,5721,5156,5200,5556,5558,5675,5669,5670,5351,5354,5355,6197,6203,6204,6205,6088,5160,5161,5162,5203,5166,5167,5189,5178,5182,5494] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.333,0.702]
select surf_pocket2, protein and id [340,2094,2095,2108,2109,84,2097,324,2096,65,325,537,728,526,2098,2101,2102,712,713,716,1836,1839,525,64,2080,1813,1826,1828,336,337,727,2238] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.443,0.361,0.902]
select surf_pocket3, protein and id [110,364,369,247,365,384,215,244,124,213,214,121,238,165,210,237,239,241,242,276] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.463,0.278,0.702]
select surf_pocket4, protein and id [3359,3881,3892,3081,3083,3893,3882,3864,3867,3858,3328,3075,3322,3857] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.741,0.361,0.902]
select surf_pocket5, protein and id [2545,2574,2782,2803,2805,2806,2778,2587,2822,2823,2804,2546,2534,2536,2541,2552,2571,2742,2743,2575,2507,2518,2589,2591,2520,2583,2584,2585,2586,2588] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.694,0.278,0.702]
select surf_pocket6, protein and id [3913,3918,3904,3908,3910,3967,3969,3968,3970,3361,3365,3965,4143,3948,3950,3351,3353,3898,3879] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.757]
select surf_pocket7, protein and id [1804,1583,1208,1584,1592,1593,1521,1522,1523,1524,1585,989,1697,1002,1200,988,1202,1198,1005] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.475]
select surf_pocket8, protein and id [1950,1952,4450,4452,4454,4471,2010,1988,4425,4423,4392,4440,4411,4407,4409,4410,4303,4418] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.361,0.459]
select surf_pocket9, protein and id [1832,1835,1836,1839,1852,1792,1795,1797,1874,68,69,525,1781,1858,1861,1759,1761,1775,1786,1787] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.318,0.278]
select surf_pocket10, protein and id [3270,3439,3440,3249,3251,3226,3227,3215,3013,3285] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.561,0.361]
select surf_pocket11, protein and id [2233,1899,1902,2060,2077,2080,2082,2227,1827,1828,1924,2218,2225,1903,1913,1915,1893,1895,1897] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.553,0.278]
select surf_pocket12, protein and id [5435,5461,5181,5183,6228,5389,5437,5356,5372,5373,5462,5424] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.859,0.361]
select surf_pocket13, protein and id [2695,2981,3217,3218,2712,2547,2711,2544,2553,2488] 
set surface_color,  pcol13, surf_pocket13 
   

deselect

orient
